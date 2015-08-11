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

Int_t nVarsprim = 12;
Int_t nVars = nVarsprim+3;

void makeSolverEffHistos()
{
  TH1::SetDefaultSumw2();

  setTDRStyle();
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout.precision(3);


	double scalettdil = 1.;
	double scalefake = 2.18495;
	double scalewjets = 1.;
	double scaleDYeemm = 1.35973;
	double scaleDYtautau = 1.17793;
	double scaletw = 1.;
	double scaleVV = 1.;


	TString yaxisunit;

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

	Float_t observable, observable_gen, tmass, ttmass, ttRapidity2,ttPt;
	Float_t observableMinus, observableMinus_gen;
	Double_t weight;
	Int_t Nsolns = 1;


	TString obs2D[3] = {"tt_mass","ttRapidity2","ttPt"};


	for (Int_t iVar = 0; iVar < nVars; iVar++)
	  {
		///////////////////////////////////////////////////////////////////////////////////////////
		/////////////// 1. Set up all our histograms //////////////////////////////////////////////

		TCanvas *c1 = new TCanvas();
		//c1->Divide(2,2);

		if(iVar<nVarsprim) Initialize2DBinning(iVar);
		else if (iVar==nVarsprim) {Initialize2DBinning(6); ybins2D[0] = 340;}
		else if (iVar==nVarsprim+1) Initialize2DBinningttrapidity2(6);
		else if (iVar==nVarsprim+2) Initialize2DBinningttpt(6);

		if(iVar>=nVarsprim) observablename = obs2D[iVar-nVarsprim];
		if(iVar>=nVarsprim)  xaxislabel = yaxislabel;




		bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;

		cout << "Now making histograms for " << acceptanceName << "..." << endl;

		//Do all our bin splitting
		int nbinsx_gen = -99;
		int nbinsx_reco = -99;
		// int nbinsunwrapped_gen = -99;
		// int nbinsunwrapped_reco = -99;

		if( iVar < 2 || iVar==9 ) nbinsx_gen = nbinsx2D*2;
		else nbinsx_gen = nbinsx2D;

		if(iVar>=nVarsprim) nbinsx_gen = 6;

		//nbinsx_reco = nbinsx_gen*2;
		//nbinsx_reco = nbinsx_gen;

		double* genbins;
		//double* recobins;

		genbins = new double[nbinsx_gen+1];
		//recobins = new double[nbinsx_reco+1];

		//Make gen binning array
		if( iVar < 2 || iVar==9 ) {
		  for( int i=0; i<=nbinsx2Dalt; i++ ) {
			genbins[i] = xbins2Dalt[i];
		  }
		}
		else if(iVar<nVarsprim)  {
		  for( int i=0; i<=nbins1D; i++ ) {
			genbins[i] = xbins2D[i];
		  }
		}/*
		else {
		  for( int i=0; i<=3; i++ ) {
			genbins[i] = ybins2D[i];
		  }
		}*/
		else {
			//Make finer binning array for ttbar system vars
			for( int i=0; i<nbinsx_gen/2; i++ ) {
				genbins[i*2] = ybins2D[i];
				genbins[i*2 +1] = ( ybins2D[i] + ybins2D[i+1] )/2.;
			}
			genbins[nbinsx_gen] = ybins2D[nbinsx_gen/2];
		}

		
		//Make histograms

		TH1D *h1Data = new TH1D ("data_1d_"+observablename, "Data_"+observablename, nbinsx_gen, genbins);
		TH1D *h1Bkg  = new TH1D ("background_1d_"+observablename, "Background_"+observablename, nbinsx_gen, genbins);
		TH1D *h1Meas = new TH1D ("meas_1d_"+observablename, "Measured_"+observablename, nbinsx_gen, genbins);

		TH1D *h1sData = new TH1D ("data_1d_s_"+observablename, "Data_withsol_"+observablename, nbinsx_gen, genbins);
		TH1D *h1sBkg  = new TH1D ("background_1d_s_"+observablename, "Background_withsol_"+observablename, nbinsx_gen, genbins);
		TH1D *h1sMeas = new TH1D ("meas_1d_s_"+observablename, "Measured_withsol_"+observablename, nbinsx_gen, genbins);


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
		//delete[] recobins;



		///// Load data from data chain, and fill hData //////////
		ch_data->SetBranchAddress(observablename,    &observable);
		if ( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
		ch_data->SetBranchAddress("weight", &weight);
		ch_data->SetBranchAddress("t_mass", &tmass);
		ch_data->SetBranchAddress("tt_mass", &ttmass);
		ch_data->SetBranchAddress("ttRapidity2", &ttRapidity2);
		ch_data->SetBranchAddress("ttPt", &ttPt);

		for (Int_t i = 0; i < ch_data->GetEntries(); i++)
		  {
			ch_data->GetEntry(i);

			  fillUnderOverFlow(h1Data, observable, weight, Nsolns);
			  if (combineLepMinus) fillUnderOverFlow(h1Data, observableMinus, weight, Nsolns);

			if ( tmass > 0 ) {
			  fillUnderOverFlow(h1sData, observable, weight, Nsolns);
			  if (combineLepMinus) fillUnderOverFlow(h1sData, observableMinus, weight, Nsolns);
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
			ch_bkg[iBkg]->SetBranchAddress("ttRapidity2", &ttRapidity2);
			ch_bkg[iBkg]->SetBranchAddress("ttPt", &ttPt);

			for (Int_t i = 0; i < ch_bkg[iBkg]->GetEntries(); i++)
			  {
				ch_bkg[iBkg]->GetEntry(i);
				weight *= bkgSF[iBkg];

				  fillUnderOverFlow(h1Bkg, observable, weight, Nsolns);
				  if (combineLepMinus) fillUnderOverFlow(h1Bkg, observableMinus, weight, Nsolns);

				if ( tmass > 0 ) {
				  fillUnderOverFlow(h1sBkg, observable, weight, Nsolns);
				  if (combineLepMinus) fillUnderOverFlow(h1sBkg, observableMinus, weight, Nsolns);
				}
			  }
		  }


		///// Load true top MC from top chain, and fill h_true and hTrue_vs_Meas ///////////
		ch_top->SetBranchAddress(observablename + "_gen",    &observable);
		if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus);
		ch_top->SetBranchAddress("weight", &weight);
		ch_top->SetBranchAddress("t_mass", &tmass);
		ch_top->SetBranchAddress("tt_mass", &ttmass);
		ch_top->SetBranchAddress("ttRapidity2", &ttRapidity2);
		ch_top->SetBranchAddress("ttPt", &ttPt);

		for (Int_t i = 0; i < ch_top->GetEntries(); i++)
		  {
			ch_top->GetEntry(i);
			weight *= scalettdil;

			  fillUnderOverFlow(h1Meas, observable, weight, Nsolns);
			  if (combineLepMinus) {
				fillUnderOverFlow(h1Meas, observableMinus, weight, Nsolns);
			  }

			if ( tmass > 0 ) {
			  fillUnderOverFlow(h1sMeas, observable, weight, Nsolns);
			  if (combineLepMinus) {
				fillUnderOverFlow(h1sMeas, observableMinus, weight, Nsolns);
			  }
			}
		  }


		/////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////// 3. Calculate the ratios (data/MC) /////////////////////////////


		TH1D *h1Bkgratio = (TH1D*)h1sBkg->Clone(observablename+"_Bkgratio");
		h1Bkgratio->Divide(h1Bkg);
		h1Bkgratio->SetLineColor(kRed);

		TH1D *h1ttbarratio = (TH1D*)h1sMeas->Clone(observablename+"_ttbarratio");
		h1ttbarratio->Divide(h1Meas);
		h1ttbarratio->SetLineColor(kBlue+1);

		h1ttbarratio->SetMinimum(0.5);
		h1ttbarratio->SetMaximum(1.0);
		h1ttbarratio->GetYaxis()->SetTitle("Solver Efficiency");
		h1ttbarratio->GetXaxis()->SetTitle(xaxislabel);


		h1Meas->Add( h1Bkg, 1. );
		h1sMeas->Add( h1sBkg, 1. );

		TH1D *h1Measratio = (TH1D*)h1sMeas->Clone(observablename+"_Measratio");
		h1Measratio->Divide(h1Meas);
		h1Measratio->SetLineColor(kGreen+1);

		TH1D *h1Dataratio = (TH1D*)h1sData->Clone(observablename+"_Dataratio");
		h1Dataratio->Divide(h1Data);
		h1Dataratio->SetLineColor(kBlack);


		Float_t Afb, AfbErr;
		GetAfb(h1ttbarratio, Afb, AfbErr);
		cout<<observablename<<" "<<Afb<<" +/- "<<AfbErr<<endl;


		//c1->cd(1);
		h1ttbarratio->Draw("hist");
		//if( iVar < 2 || iVar==9 ) h1Measratio->Draw("histsame");
		//if( iVar < 2 || iVar==9 ) h1Dataratio->Draw("histsame");
		//if( iVar < 2 || iVar==9 )  h1Bkgratio->Draw("histsame");

		if(iVar<nVarsprim) c1->Print("SolverEfficiency_"+acceptanceName+".pdf");
		else c1->Print("SolverEfficiency_"+observablename+".pdf");


		ch_data->Delete();

		ch_top->Delete();

		for (int iBkg = 0; iBkg < nBkg; ++iBkg)
		  {
			ch_bkg[iBkg]->Delete();
		  }

	  }  //End loop over asymmetry variables



}  //End main function

#ifndef __CINT__
int main ()
{
  makeSolverEffHistos();    // Main program when run stand-alone
  return 0;
}
#endif
