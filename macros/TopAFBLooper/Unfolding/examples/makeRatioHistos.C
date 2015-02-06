#include <iostream>
#include <fstream>
#include <cmath>
#include "AfbFinalUnfold.h"

#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TColor.h"
#include "TCut.h"
#include "TPaveText.h"

#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================

Int_t nVars = 12;

void makeRatioHistos()
{
  TH1::SetDefaultSumw2();

  setTDRStyle();
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout.precision(3);

  TFile* outfile1D = new TFile("ratios_1D.root", "RECREATE");

  for( int jVar=0; jVar<3; jVar++ ) {

	TString Var2D;

	switch( jVar ) {
	case 0:  Var2D = "mtt";  break;
	case 1:  Var2D = "ttrapidity2"; break;
	case 2:  Var2D = "ttpt"; break;
	}

	double scalettdil = 1.;
	double scalefake = 2.27055;
	double scalewjets = 1.;
	double scaleDYeemm = 1.46211;
	double scaleDYtautau = 1.17888;
	double scaletw = 1.;
	double scaleVV = 1.;

	// TString outfile2d_name = "ratio_histograms_"+Var2D+".root";
	// TString outfile1d_name = "ratio_histograms_1D.root";

	TString yaxisunit;
	if (Var2D == "mtt") yaxisunit = " (GeV/c^{2})";
	else if (Var2D == "ttrapidity2") yaxisunit = "";
	else if (Var2D == "ttpt") yaxisunit = " (GeV/c)";


	const int nBkg = 11;
	const int nSig = 3;
	TString path = "../";
	TString dataroot[nSig] = {"data_diel_baby.root", "data_dimu_baby.root", "data_mueg_baby.root"};
	TString bkgroot[nBkg];
	double bkgSF[nBkg];

	bkgroot[0] = "DY1to4Jeemm_baby.root";     bkgSF[0] = scaleDYeemm;
	bkgroot[1] = "DY1to4Jtautau_baby.root";   bkgSF[1] = scaleDYtautau;
	bkgroot[2] = "diboson_baby.root";		  bkgSF[2] = scaleVV;
	bkgroot[3] = "tW_lepdl_baby.root";	      bkgSF[3] = scaletw;
	bkgroot[4] = "tW_lepfake_baby.root";	  bkgSF[4] = scalefake;
	bkgroot[5] = "tW_lepsl_baby.root";		  bkgSF[5] = scalefake;
	bkgroot[6] = "triboson_baby.root";		  bkgSF[6] = scaleVV;
	bkgroot[7] = "ttV_baby.root";			  bkgSF[7] = scaleVV;
	bkgroot[8] = "ttfake_mcatnlo_baby.root";  bkgSF[8] = scalefake;
	bkgroot[9] = "ttsl_mcatnlo_baby.root";	  bkgSF[9] = scalefake;
	bkgroot[10] = "w1to4jets_baby.root";	  bkgSF[10] = scalefake;

	Float_t observable, observable_gen, tmass, ttmass, ttRapidity2;
	Float_t observableMinus, observableMinus_gen;
	Float_t obs2D, obs2D_gen;
	Double_t weight;
	Int_t Nsolns = 1;

	TFile* outfile2D = new TFile("ratios_"+Var2D+".root", "RECREATE");

	for (Int_t iVar = 0; iVar < nVars; iVar++)
	  {
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////// 1. Set up all our histograms //////////////////////////////////////////////

		//initialise 2D binning
		if (Var2D == "mtt") Initialize2DBinning(iVar);
		else if (Var2D == "ttrapidity2") Initialize2DBinningttrapidity2(iVar);
		else if (Var2D == "ttpt") Initialize2DBinningttpt(iVar);
		bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;

		cout << "Now making ratio histograms for " << acceptanceName << " vs " << Var2D << "..." << endl;

		//Do all our bin splitting
		int nbinsx_gen = -99;
		int nbinsx_reco = -99;
		// int nbinsunwrapped_gen = -99;
		// int nbinsunwrapped_reco = -99;

		if( iVar < 2 || iVar==9 ) nbinsx_gen = nbinsx2D*2;
		else nbinsx_gen = nbinsx2D;

		nbinsx_reco = nbinsx_gen*2;
		//nbinsx_reco = nbinsx_gen;

		double* genbins;
		double* recobins;

		genbins = new double[nbinsx_gen+1];
		recobins = new double[nbinsx_reco+1];

		//Make gen binning array
		if( iVar < 2 || iVar == 9 ) {
		  for( int i=0; i<=nbinsx2Dalt; i++ ) {
			genbins[i] = xbins2Dalt[i];
		  }
		}
		else {
		  for( int i=0; i<=nbins1D; i++ ) {
			genbins[i] = xbins2D[i];
		  }
		}

		//Make reco binning array
		for( int i=0; i<nbinsx_gen; i++ ) {
		  if( nbinsx_reco > nbinsx_gen ) {
			recobins[i*2] = genbins[i];
			recobins[i*2 +1] = ( genbins[i] + genbins[i+1] )/2.;
		  }
		  else recobins[i] = genbins[i];
		}
		recobins[nbinsx_reco] = genbins[nbinsx_gen];

		
		//Make histograms
		TH2D *h2Data = new TH2D ("Data_2d", "Data", nbinsx_reco, recobins, nbinsy2D, ybins2D);
		TH2D *h2Bkg = new TH2D ("Background_2d",  "Background", nbinsx_reco, recobins, nbinsy2D, ybins2D);
		TH2D *h2True = new TH2D ("true_2d", "Truth",    nbinsx_gen, genbins, nbinsy2D, ybins2D);
		TH2D *h2Meas = new TH2D ("meas_2d", "Measured", nbinsx_reco, recobins, nbinsy2D, ybins2D);

		TH1D *h1Data = new TH1D ("data_1d", "Data", nbinsx_reco, recobins);
		TH1D *h1Bkg  = new TH1D ("background_1d", "Background", nbinsx_reco, recobins);
		TH1D *h1True = new TH1D ("true_1d", "Truth", nbinsx_gen, genbins);
		TH1D *h1Meas = new TH1D ("meas_1d", "Measured", nbinsx_reco, recobins);


		/////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////// 2. Fill our histograms from the baby ntuples //////////////////

		// Set up chains
		TChain *ch_bkg[nBkg];
		TChain *ch_top = new TChain("tree");
		TChain *ch_data = new TChain("tree");

		// ch_data->Add(path + "data.root");
		for (int iSig = 0; iSig < nSig; ++iSig)
		  {
			ch_data->Add(path + dataroot[iSig]);
		  }

		ch_top->Add(path + "ttdl_mcatnlo_baby.root");

		for (int iBkg = 0; iBkg < nBkg; ++iBkg)
		  {
			ch_bkg[iBkg] = new TChain("tree");
			ch_bkg[iBkg]->Add(path + bkgroot[iBkg]);
		  }

		delete[] genbins;
		delete[] recobins;



		///// Load data from data chain, and fill hData //////////
		ch_data->SetBranchAddress(observablename,    &observable);
		if ( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
		ch_data->SetBranchAddress("weight", &weight);
		ch_data->SetBranchAddress("t_mass", &tmass);
		ch_data->SetBranchAddress("tt_mass", &ttmass);
		// ch_data->SetBranchAddress("ttRapidity2", &ttRapidity2);

		if (Var2D == "mtt")               ch_data->SetBranchAddress("tt_mass", &obs2D);
		else if (Var2D == "ttrapidity2")  ch_data->SetBranchAddress("ttRapidity2", &obs2D);
		else if (Var2D == "ttpt")         ch_data->SetBranchAddress("ttPt", &obs2D);

		for (Int_t i = 0; i < ch_data->GetEntries(); i++)
		  {
			ch_data->GetEntry(i);
			obs2D = fabs(obs2D);

			if ( tmass > 0 ) {
			  fillUnderOverFlow(h2Data, observable, obs2D, weight, Nsolns);
			  if (combineLepMinus) fillUnderOverFlow(h2Data, observableMinus, obs2D, weight, Nsolns);
			}
			if ( iVar < 2 || iVar == 9 || tmass > 0 ) {
			  fillUnderOverFlow(h1Data, observable, weight, Nsolns);
			  if (combineLepMinus) fillUnderOverFlow(h1Data, observableMinus, weight, Nsolns);
			}
		  }


		///// Load background MC from background chain, and fill h_bkg //////////
		for (int iBkg = 0; iBkg < nBkg; ++iBkg)
		  {
			ch_bkg[iBkg]->SetBranchAddress(observablename,    &observable);
			if ( combineLepMinus ) ch_bkg[iBkg]->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
			ch_bkg[iBkg]->SetBranchAddress("weight", &weight);
			ch_bkg[iBkg]->SetBranchAddress("t_mass", &tmass);
			ch_bkg[iBkg]->SetBranchAddress("tt_mass", &ttmass);
			// ch_bkg[iBkg]->SetBranchAddress("ttRapidity2", &ttRapidity2);

			if (Var2D == "mtt")              ch_bkg[iBkg]->SetBranchAddress("tt_mass", &obs2D);
			else if (Var2D == "ttrapidity2") ch_bkg[iBkg]->SetBranchAddress("ttRapidity2", &obs2D);
			else if (Var2D == "ttpt")        ch_bkg[iBkg]->SetBranchAddress("ttPt", &obs2D);

			for (Int_t i = 0; i < ch_bkg[iBkg]->GetEntries(); i++)
			  {
				ch_bkg[iBkg]->GetEntry(i);
				obs2D = fabs(obs2D);
				weight *= bkgSF[iBkg];

				if ( tmass > 0 ) {
				  fillUnderOverFlow(h2Bkg, observable, obs2D, weight, Nsolns);
				  if (combineLepMinus) fillUnderOverFlow(h2Bkg, observableMinus, obs2D, weight, Nsolns);
				}
				if ( iVar < 2 || iVar == 9 || tmass > 0 ) {
				  fillUnderOverFlow(h1Bkg, observable, weight, Nsolns);
				  if (combineLepMinus) fillUnderOverFlow(h1Bkg, observableMinus, weight, Nsolns);
				}
			  }
		  }


		///// Load true top MC from top chain, and fill h_true and hTrue_vs_Meas ///////////
		ch_top->SetBranchAddress(observablename,    &observable);
		ch_top->SetBranchAddress(observablename + "_gen", &observable_gen);
		if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
		if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
		ch_top->SetBranchAddress("weight", &weight);
		ch_top->SetBranchAddress("t_mass", &tmass);
		ch_top->SetBranchAddress("tt_mass", &ttmass);
		// ch_top->SetBranchAddress("ttRapidity2", &ttRapidity2);

		if (Var2D == "mtt") {
		  ch_top->SetBranchAddress("tt_mass", &obs2D);
		  ch_top->SetBranchAddress("tt_mass_gen", &obs2D_gen);
		}
		else if (Var2D == "ttrapidity2") {
		  ch_top->SetBranchAddress("ttRapidity2", &obs2D);
		  ch_top->SetBranchAddress("ttRapidity2_gen", &obs2D_gen);
		}
		else if (Var2D == "ttpt") {
		  ch_top->SetBranchAddress("ttPt", &obs2D);
		  ch_top->SetBranchAddress("ttPt_gen", &obs2D_gen);
		}

		for (Int_t i = 0; i < ch_top->GetEntries(); i++)
		  {
			ch_top->GetEntry(i);
			obs2D = fabs(obs2D);
			obs2D_gen = fabs(obs2D_gen);
			weight *= scalettdil;

			if ( tmass > 0 ) {
			  fillUnderOverFlow(h2Meas, observable, obs2D, weight, Nsolns);
			  fillUnderOverFlow(h2True, observable_gen, obs2D_gen, weight, Nsolns);
			  if ( combineLepMinus ) {
				fillUnderOverFlow(h2Meas, observableMinus, obs2D, weight, Nsolns);
				fillUnderOverFlow(h2True, observableMinus_gen, obs2D_gen, weight, Nsolns);
			  }
			}
			if ( iVar < 2 || iVar == 9 || tmass > 0 ) {
			  fillUnderOverFlow(h1Meas, observable, weight, Nsolns);
			  fillUnderOverFlow(h1True, observable_gen, weight, Nsolns);
			  if (combineLepMinus) {
				fillUnderOverFlow(h1Meas, observableMinus, weight, Nsolns);
				fillUnderOverFlow(h1True, observableMinus_gen, weight, Nsolns);
			  }
			}
		  }


		/////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////// 3. Calculate the ratios (data/MC) /////////////////////////////

		TH2D *h2ratio = (TH2D*)h2Data->Clone(acceptanceName+"_"+Var2D+"_ratio_2d");
		TH1D *h1ratio = (TH1D*)h1Data->Clone(acceptanceName+"_ratio_1d");

		h2ratio->Divide(h2Meas);
		h1ratio->Divide(h1Meas);

		//////////////////////////////////
		// Write ratio histograms to files

		outfile2D->cd();
		h2ratio->Write();

		if( jVar==1) {    // Only write the 1D histograms once
		  outfile1D->cd();
		  h1ratio->Write();
		}



		ch_data->Delete();

		ch_top->Delete();

		for (int iBkg = 0; iBkg < nBkg; ++iBkg)
		  {
			ch_bkg[iBkg]->Delete();
		  }

	  }  //End loop over asymmetry variables

	outfile2D->Close();

  }  //End loop over secondary variables

  outfile1D->Close();

}  //End main function

#ifndef __CINT__
int main ()
{
  makeRatioHistos();    // Main program when run stand-alone
  return 0;
}
#endif
