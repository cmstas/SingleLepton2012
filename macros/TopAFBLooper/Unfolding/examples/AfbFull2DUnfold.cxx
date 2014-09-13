#include <iostream>
#include <fstream>
#include <cmath>
#include "AfbFinalUnfold.h"

#include "TRandom3.h"
#include "TH1.h"
#include "TH2.h"
#include "TStyle.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TMatrixD.h"
#include "TChain.h"
#include "TLegend.h"
#include "TColor.h"
#include "THStack.h"
#include "TCut.h"
#include "TMarker.h"
#include "TUnfoldSys.h"
#include "TPaveText.h"

#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================

//const Double_t _topScalingFactor=1.+(9824. - 10070.94)/9344.25;  //needs to be changed from preselection ratio to ratio for events with a ttbar solution


Int_t kterm = 3;  //note we used 4 here for ICHEP
Double_t tau = 3E-2;
Int_t nVars = 12;
Int_t includeSys = 0;
Int_t checkErrors = 0;
bool draw_truth_before_pT_reweighting = false;


void AfbUnfoldExample(TString Var2D = "mtt", double scalettdil = 1., double scalettotr = 1., double scalewjets = 1., double scaleDY = 1., double scaletw = 1., double scaleVV = 1.)
{
    TH1::SetDefaultSumw2();

    setTDRStyle();
    gStyle->SetOptFit();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    cout.precision(3);

    TString summary_name;
    if (Var2D == "mtt") summary_name = "summary_2Dunfolding";
    else if (Var2D == "ttrapidity2") summary_name = "summary_2Dunfolding_ttrapidity2";
    else if (Var2D == "ttpt") summary_name = "summary_2Dunfolding_ttpt";

    TString yaxisunit;
    if (Var2D == "mtt") yaxisunit = " (GeV/c^{2})";
    else if (Var2D == "ttrapidity2") yaxisunit = "";
    else if (Var2D == "ttpt") yaxisunit = " (GeV/c)";

    if (!(scalettotr == 1. && scalewjets == 1. && scaleDY == 1. && scaletw == 1. && scaleVV == 1.))  summary_name = summary_name + Form("_%i_%i_%i_%i_%i", int(10.*scalettotr + 0.5), int(10.*scalewjets + 0.5), int(10.*scaleDY + 0.5), int(10.*scaletw + 0.5), int(10.*scaleVV + 0.5));

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);

    ofstream myfile;
    myfile.open (summary_name + ".txt");
    cout.rdbuf(myfile.rdbuf());

    // OGU 130516: add second output txt file with format easier to be pasted into google docs
    ofstream second_output_file;
    second_output_file.open(summary_name + "_formated.txt");

    const int nBkg = 10;
	const int nSig = 3;
    TString path = "../";
    TString dataroot[nSig] = {"data_diel_baby.root", "data_dimu_baby.root", "data_mueg_baby.root"};
    // TString bkgroot[nBkg] = {"ttotr.root", "wjets.root", "DYee.root", "DYmm.root", "DYtautau.root", "tw.root", "VV.root"};
    TString bkgroot[nBkg] = {"DY1to4Jtot_baby.root", "diboson_baby.root", "tW_lepdl_baby.root", "tW_lepfake_baby.root", "tW_lepsl_baby.root", "triboson_baby.root", "ttV_baby.root", "ttfake_powheg_baby.root", "ttsl_powheg_baby.root", "w1to4jets_baby.root"};

    double bkgSF[nBkg] = {scalettotr, scalewjets, scaleDY, scaleDY, scaleDY, scaletw, scaleVV};

    Float_t observable, observable_gen, tmass, ttmass, ttRapidity2;
    Float_t observableMinus, observableMinus_gen;
    Float_t obs2D, obs2D_gen;
    Double_t weight;
    Int_t Nsolns = 1;
	Int_t channel = -99;

    for (Int_t iVar = 0; iVar < nVars; iVar++)
    {
	  ///////////////////////////////////////////////////////////////////////////////////////////
	  /////////////// 1. Set up all our histograms //////////////////////////////////////////////

        //initialise 2D binning
        if (Var2D == "mtt") Initialize2DBinning(iVar);
        else if (Var2D == "ttrapidity2") Initialize2DBinningttrapidity2(iVar);
        else if (Var2D == "ttpt") Initialize2DBinningttpt(iVar);
        bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;
		
		//Do all our bin splitting
		int nbinsx_gen = -99;
		int nbinsx_reco = -99;
		int nbinsx_reco_3ch = -99;
		int nbinsunwrapped_gen = -99;
		int nbinsunwrapped_reco = -99;
		int nbinsunwrapped_reco_3ch = -99;

		if( iVar < 2 || iVar==9 ) nbinsx_gen = nbinsx2D*2;
		else nbinsx_gen = nbinsx2D;

		nbinsx_reco = nbinsx_gen*2;
		//nbinsx_reco = nbinsx_gen;
		nbinsx_reco_3ch = nbinsx_reco*3;

		double* genbins;
		double* recobins;
		double* recobins_3ch;

		genbins = new double[nbinsx_gen+1];
		recobins = new double[nbinsx_reco+1];
		recobins_3ch = new double[nbinsx_reco_3ch+1];

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

		// Make reco binning array including the 3-channel split
		double recohist_width = fabs(recobins[nbinsx_reco] - recobins[0]);
		std::copy( recobins, recobins+nbinsx_reco+1, recobins_3ch );
		for( int i=nbinsx_reco+1; i<nbinsx_reco_3ch+1; i++ ) {
		  recobins_3ch[i] = recobins_3ch[i-nbinsx_reco] + recohist_width;
		}

		nbinsunwrapped_gen  = nbinsx_gen  * nbinsy2D;
		nbinsunwrapped_reco = nbinsx_reco * nbinsy2D;
		nbinsunwrapped_reco_3ch = nbinsx_reco_3ch * nbinsy2D;
		
		//Make histograms

		//Use proper 2D histograms, instead of the old 1D ones
        TH2D *hData_combined = new TH2D ("Data_combined", "Data combined", nbinsx_reco, recobins, nbinsy2D, ybins2D);
        TH2D *hData = new TH2D ("Data", "Data", nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
        TH2D *hBkg = new TH2D ("Background",  "Background",    nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
        TH2D *hBkg_combined = new TH2D ("Background_combined",  "Background combined", nbinsx_reco, recobins, nbinsy2D, ybins2D);
        TH2D *hData_unfolded = new TH2D ("Data_Unfold", "Data with background subtracted and unfolded", nbinsx_gen, genbins, nbinsy2D, ybins2D);
        TH2D *hTrue = new TH2D ("true", "Truth",    nbinsx_gen, genbins, nbinsy2D, ybins2D);
        TH2D *hTrue_split = new TH2D ("true_split", "Truth",    nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
        TH2D *hMeas = new TH2D ("meas", "Measured", nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
		TH2D *hPurity = new TH2D("purity", "Purity", nbinsx_gen, genbins, nbinsy2D, ybins2D);
		TH2D *hStability = new TH2D("stability", "Stability", nbinsx_gen, genbins, nbinsy2D, ybins2D);

		//Unwrapped histograms have n bins (where n = nx*ny), centered around the integers from 1 to n.
		TH1D *hData_unwrapped = new TH1D ("Data_BkgSub_Unwr", "Unwrapped data with background subtracted", nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5);
        TH1D *hBkg_unwrapped = new TH1D ("Background_Unwr",  "Background unwrapped",    nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5);
        TH1D *hData_unfolded_unwrapped = new TH1D ("Data_Unfold_Unwr", "Data unfolded and unwrapped", nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);
        TH1D *hTrue_unwrapped = new TH1D ("true_unwr", "Truth unwrapped",  nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);
        TH1D *hTrue_unwrapped_split = new TH1D ("true_unwr_split", "Truth unwrapped, split",  nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5);
        TH1D *hMeas_unwrapped = new TH1D ("meas_unwr", "Measured unwrapped", nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5);

		//Migration matrix, using the unwrapped binning on both axes
        TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5, nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);

        TH1D *hData_bkgSub;
		TH2D *hData_bkgSub_combined;
		// TH1D* hMeas_newErr;

		delete[] genbins;
		delete[] recobins_3ch;

		hTrue_split->RebinX(2);

        TMatrixD m_unfoldE(nbinsunwrapped_gen, nbinsunwrapped_gen);
        TMatrixD m_correctE(nbinsunwrapped_gen, nbinsunwrapped_gen);

		/////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////// 2. Fill our histograms from the baby ntuples //////////////////

		//cout << "Filling histo for " << acceptanceName << "." << endl;
		//cout << "X-axis: " << xaxislabel << endl;
		//cout << "Y-axis: " << yaxislabel << endl;

		// Set up chains
        TChain *ch_bkg[nBkg];
        TChain *ch_top = new TChain("tree");
        TChain *ch_data = new TChain("tree");

        // ch_data->Add(path + "data.root");
        for (int iSig = 0; iSig < nSig; ++iSig)
        {
            ch_data->Add(path + dataroot[iSig]);
        }

        ch_top->Add(path + "ttdl_mcatnlo_smallTree_baby.root");

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg] = new TChain("tree");
            ch_bkg[iBkg]->Add(path + bkgroot[iBkg]);
        }

		double offset = 0;
		double histmax = recobins[nbinsx_reco];
		double histmin = recobins[0];
		double hiBinCenter = hData_combined->GetXaxis()->GetBinCenter(nbinsx_reco);
		double loBinCenter = hData_combined->GetXaxis()->GetBinCenter(1);

		///// Load data from data chain, and fill hData //////////
        ch_data->SetBranchAddress(observablename,    &observable);
        if ( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        ch_data->SetBranchAddress("weight", &weight);
        // ch_data->SetBranchAddress("Nsolns", &Nsolns);
        ch_data->SetBranchAddress("t_mass", &tmass);
        ch_data->SetBranchAddress("tt_mass", &ttmass);
        ch_data->SetBranchAddress("ttRapidity2", &ttRapidity2);
		ch_data->SetBranchAddress("channel", &channel);

        if (Var2D == "mtt")               ch_data->SetBranchAddress("tt_mass", &obs2D);
        else if (Var2D == "ttrapidity2")  ch_data->SetBranchAddress("ttRapidity2", &obs2D);
        else if (Var2D == "ttpt")         ch_data->SetBranchAddress("ttPt", &obs2D);

        for (Int_t i = 0; i < ch_data->GetEntries(); i++)
        {
            ch_data->GetEntry(i);
            obs2D = fabs(obs2D);

			//Use an offset to sort events into 3 superbins: ee, mumu, emu
			offset = double(channel) * recohist_width;
			//Do the same thing as "fillUnderOverflow", except adapted for 3x1 histograms
			if( observable > histmax )        observable = hiBinCenter;
			else if( observable < histmin )   observable = loBinCenter;
			if( observableMinus > histmax )        observableMinus = hiBinCenter;
			else if( observableMinus < histmin )   observableMinus = loBinCenter;

            if ( tmass > 0 )
            {
			  fillUnderOverFlow(hData_combined, observable, obs2D, weight, Nsolns);
			  fillUnderOverFlow(hData, observable+offset, obs2D, weight, Nsolns);
			  if (combineLepMinus) {
				fillUnderOverFlow(hData_combined, observableMinus, obs2D, weight, Nsolns);
				fillUnderOverFlow(hData, observableMinus+offset, obs2D, weight, Nsolns);
			  }
            }
        }

		///// Load background MC from background chain, and fill h_bkg //////////
        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg]->SetBranchAddress(observablename,    &observable);
            if ( combineLepMinus ) ch_bkg[iBkg]->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
            ch_bkg[iBkg]->SetBranchAddress("weight", &weight);
            // ch_bkg[iBkg]->SetBranchAddress("Nsolns", &Nsolns);
            ch_bkg[iBkg]->SetBranchAddress("t_mass", &tmass);
			ch_bkg[iBkg]->SetBranchAddress("tt_mass", &ttmass);
			ch_bkg[iBkg]->SetBranchAddress("ttRapidity2", &ttRapidity2);
			ch_bkg[iBkg]->SetBranchAddress("channel", &channel);

            if (Var2D == "mtt")              ch_bkg[iBkg]->SetBranchAddress("tt_mass", &obs2D);
            else if (Var2D == "ttrapidity2") ch_bkg[iBkg]->SetBranchAddress("ttRapidity2", &obs2D);
            else if (Var2D == "ttpt")        ch_bkg[iBkg]->SetBranchAddress("ttPt", &obs2D);

            for (Int_t i = 0; i < ch_bkg[iBkg]->GetEntries(); i++)
            {
                ch_bkg[iBkg]->GetEntry(i);
                obs2D = fabs(obs2D);
                weight *= bkgSF[iBkg];

				offset = double(channel) * recohist_width;
				if( observable > histmax )        observable = hiBinCenter;
				else if( observable < histmin )   observable = loBinCenter;
				if( observableMinus > histmax )        observableMinus = hiBinCenter;
				else if( observableMinus < histmin )   observableMinus = loBinCenter;

                if ( tmass > 0 )
                {
				  fillUnderOverFlow(hBkg, observable+offset, obs2D, weight, Nsolns);
				  fillUnderOverFlow(hBkg_combined, observable, obs2D, weight, Nsolns);
				  if (combineLepMinus) {
					fillUnderOverFlow(hBkg, observableMinus+offset, obs2D, weight, Nsolns);
					fillUnderOverFlow(hBkg_combined, observableMinus, obs2D, weight, Nsolns);
				  }
                }
            }
        }

		///// Load true top MC from top chain, and fill h_true and hTrue_vs_Meas ///////////
        ch_top->SetBranchAddress(observablename,    &observable);
        ch_top->SetBranchAddress(observablename + "_gen", &observable_gen);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
        ch_top->SetBranchAddress("weight", &weight);
        // ch_top->SetBranchAddress("Nsolns", &Nsolns);
        ch_top->SetBranchAddress("t_mass", &tmass);
        ch_top->SetBranchAddress("tt_mass", &ttmass);
        ch_top->SetBranchAddress("ttRapidity2", &ttRapidity2);
		ch_top->SetBranchAddress("channel", &channel);

        if (Var2D == "mtt")
        {
            ch_top->SetBranchAddress("tt_mass", &obs2D);
            ch_top->SetBranchAddress("tt_mass_gen", &obs2D_gen);
        }
        else if (Var2D == "ttrapidity2")
        {
            ch_top->SetBranchAddress("ttRapidity2", &obs2D);
            ch_top->SetBranchAddress("ttRapidity2_gen", &obs2D_gen);
        }
        else if (Var2D == "ttpt")
        {
            ch_top->SetBranchAddress("ttPt", &obs2D);
            ch_top->SetBranchAddress("ttPt_gen", &obs2D_gen);
        }

		int measbin_3ch = -99;
		int measbin = -99;
		int genbin = -99;

        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
        {
            ch_top->GetEntry(i);
            obs2D = fabs(obs2D);
            obs2D_gen = fabs(obs2D_gen);
            weight *= scalettdil;

			offset = double(channel) * recohist_width;
			if( observable > histmax )        observable = hiBinCenter;
			else if( observable < histmin )   observable = loBinCenter;
			if( observableMinus > histmax )        observableMinus = hiBinCenter;
			else if( observableMinus < histmin )   observableMinus = loBinCenter;
			if( observable_gen > histmax )        observable_gen = hiBinCenter;
			else if( observable_gen < histmin )   observable_gen = loBinCenter;
			if( observableMinus_gen > histmax )        observableMinus_gen = hiBinCenter;
			else if( observableMinus_gen < histmin )   observableMinus_gen = loBinCenter;

            if ( tmass > 0 )
            {
			  genbin  = getUnwrappedBin(hTrue, observable_gen, obs2D_gen);
			  measbin = getUnwrappedBin(hData_combined, observable, obs2D);
			  measbin_3ch = measbin + (channel * nbinsunwrapped_reco);

			  // cout << "Channel: " << channel << ". Offset: " << offset << endl;

			  fillUnderOverFlow(hMeas, observable+offset, obs2D, weight, Nsolns);
			  fillUnderOverFlow(hTrue, observable_gen, obs2D_gen, weight, Nsolns);
			  fillUnderOverFlow(hTrue_split, observable_gen+offset, obs2D_gen, weight, Nsolns);
			  fillUnderOverFlow(hTrue_vs_Meas, measbin_3ch, genbin, weight, Nsolns);
			  if ( combineLepMinus )
                {
				  genbin  = getUnwrappedBin(hTrue, observableMinus_gen, obs2D_gen);
				  measbin = getUnwrappedBin(hData_combined, observableMinus, obs2D);
				  measbin_3ch = measbin + (channel * nbinsunwrapped_reco);

				  fillUnderOverFlow(hMeas, observableMinus+offset, obs2D, weight, Nsolns);
				  fillUnderOverFlow(hTrue, observableMinus_gen, obs2D_gen, weight, Nsolns);
				  fillUnderOverFlow(hTrue_split, observableMinus_gen+offset, obs2D_gen, weight, Nsolns);
				  fillUnderOverFlow(hTrue_vs_Meas, measbin_3ch, genbin, weight, Nsolns);
                }
            }
        }

		delete[] recobins;

		// Do the acceptance correction, by filling the migration matrix with events that have a gen-level value but no reco-level value
        TFile *file = new TFile("../denominator/acceptance/mcnlo/accept_" + acceptanceName + ".root");

		TH2D *acceptM[4];
		acceptM[0] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_diel"));
		acceptM[1] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_dimu"));
		acceptM[2] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_mueg"));
		acceptM[3] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_all" ));

		acceptM[3]->Scale(1.0 / acceptM[3]->Integral());

		for( int aChannel=0; aChannel<3; aChannel++ ) {
		  for( int yBin=1; yBin<=nbinsy2D; yBin++ ) {
			for( int xBin=1; xBin<=nbinsx_gen; xBin++ ) {

			  double acceptance = acceptM[aChannel]->GetBinContent( xBin, yBin );
			  double n_accepted = hTrue_split->GetBinContent( aChannel*nbinsx_gen + xBin, yBin );
			  genbin = (yBin-1)*nbinsx_gen + xBin;
			  hTrue_vs_Meas->Fill( -999999, double(genbin), n_accepted/acceptance - n_accepted );

			}
		  }
		}


		// Fill purity and stability plots
		if( hData->GetNbinsX() == nbinsx_gen ) {

		  for( int i=1; i<=nbinsx_gen; i++) {
			for( int j=1; j<=nbinsy2D; j++) {
			  int k = (j-1)*nbinsx_gen + i; //bin number in the unwrapped version
			  hPurity->SetBinContent( i, j, hTrue_vs_Meas->GetBinContent(k,k) / hMeas->GetBinContent(i,j) );
			  hStability->SetBinContent( i, j, hTrue_vs_Meas->GetBinContent(k,k) / hTrue->GetBinContent(i,j) );
			}
		  }

		}
		else if( hData->GetNbinsX() == nbinsx_reco_3ch ) {

		  for( int by=1; by<=nbinsy2D; by++ ) {
			for( int bx=1; bx<=nbinsx_gen; bx++ ) {

			  int outbin = (by-1)*nbinsx_gen + bx;
			  double numerator = 0;
			  double purity_denom = 0;
			  double stability_denom = hTrue->GetBinContent( bx, by );

			  for( int aChannel=0; aChannel<3; aChannel++ ) {
				numerator += hTrue_vs_Meas->GetBinContent((outbin*2)+(aChannel*nbinsunwrapped_reco), outbin);
				numerator += hTrue_vs_Meas->GetBinContent((outbin*2)+(aChannel*nbinsunwrapped_reco)-1, outbin);
				purity_denom += hMeas->GetBinContent( (bx*2)+(aChannel*nbinsx_reco), by );
				purity_denom += hMeas->GetBinContent( (bx*2)+(aChannel*nbinsx_reco)-1, by );
			  }

			  hPurity->SetBinContent( bx, by, numerator / purity_denom );
			  hStability->SetBinContent( bx, by, numerator / stability_denom );
			}
		  }
		}
		else {
		  cout << "\n***WARNING: Purity and stability plots are broken, because nbinsx_gen != nbinsx_reco!!!\n" << endl;
		}


		////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////// 3. Perform the unfolding procedure ////////////////////////////////////

		// First step: unwrap our TH2Ds into TH1Ds!
		unwrap2dhisto_3ch(hData, hData_unwrapped);
		unwrap2dhisto_3ch(hBkg,  hBkg_unwrapped);
		unwrap2dhisto_3ch(hTrue, hTrue_unwrapped);
		unwrap2dhisto_3ch(hTrue_split, hTrue_unwrapped_split);
		unwrap2dhisto_3ch(hMeas, hMeas_unwrapped);

        hData_bkgSub = (TH1D *) hData_unwrapped->Clone("data_bkgsub");
        hData_bkgSub->Add(hBkg_unwrapped, -1.0);

		hData_bkgSub_combined = (TH2D*) hData_combined->Clone("data_bkgsub_combined");
		hData_bkgSub_combined->Add(hBkg_combined, -1.0);
		

		/*
		hMeas_newErr = (TH1D *) hMeas_unwrapped->Clone();
		hMeas_newErr->Scale( hData_bkgSub->Integral() / hMeas_newErr->Integral() );
		for( int i=1; i<nbinsunwrapped_reco+1; i++) {
		  double n_sig = hMeas_unwrapped->GetBinContent(i);
		  double n_bkg = hBkg_unwrapped->GetBinContent(i);
		  double bkg_err = hBkg_unwrapped->GetBinError(i);
		  double mcerr = hMeas_newErr->GetBinError(i);
		  hMeas_newErr->SetBinError(i, sqrt(n_sig + n_bkg + bkg_err*bkg_err ) );
		  //hMeas_newErr->SetBinError(i, 2.0);
		  //hMeas_newErr->SetBinError(i, mcerr*10.);
		}
		*/


		// Now let's actually do the unfolding.
		TUnfoldSys unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone, TUnfold::kEConstraintArea);  //need to set reg mode "None" here if regularizing by hand
		unfold_TUnfold.SetInput(hData_bkgSub);
		scaleBias = hData_bkgSub->Integral() / hTrue_unwrapped->Integral();
		hTrue_unwrapped->Scale(scaleBias);
		//unfold_TUnfold.SetBias(hTrue_unwrapped);
		unfold_TUnfold.RegularizeBins2D(1,1,nbinsx_gen,nbinsx_gen,nbinsy2D,TUnfold::kRegModeCurvature);
		minimizeRhoAverage(&unfold_TUnfold, hData_bkgSub, -5.0, -0.5);
		unfold_TUnfold.GetOutput(hData_unfolded_unwrapped);
		tau = unfold_TUnfold.GetTau();

		TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
		for (Int_t cmi = 0; cmi < nbinsunwrapped_gen; cmi++)
		  {
			for (Int_t cmj = 0; cmj < nbinsunwrapped_gen; cmj++)
			  {
				m_unfoldE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
				m_correctE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
			  }
		  }
        

		// Generate a curve of rhoAvg vs log(tau)
		double ar_tau[100];
		double ar_rhoAvg[100];
		double tau_test = 0.0;
		double bestrhoavg = unfold_TUnfold.GetRhoAvg();

		TUnfoldSys unfold_getRhoAvg(hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone, TUnfold::kEConstraintArea);
		unfold_getRhoAvg.SetInput(hData_bkgSub);
		//unfold_getRhoAvg.SetBias(hTrue_unwrapped);
		unfold_getRhoAvg.RegularizeBins2D(1,1,nbinsx_gen,nbinsx_gen,nbinsy2D,TUnfold::kRegModeCurvature);

		for(int l=0; l<100; l++) {
		  tau_test = pow( 10, -5.5 + 0.05*l);
		  unfold_getRhoAvg.DoUnfold(tau_test, hData_bkgSub, scaleBias);
		  ar_tau[l] = tau_test;
		  ar_rhoAvg[l] = unfold_getRhoAvg.GetRhoAvg();
		}

		TGraph* gr_rhoAvg = new TGraph(100,ar_tau,ar_rhoAvg);
		TCanvas* c_rhoAvg = new TCanvas("c_rhoAvg","c_rhoAvg");
		c_rhoAvg->SetLogx();
		gr_rhoAvg->SetTitle("Global Correlation Coefficient;#tau;#rho_{avg}");
		gr_rhoAvg->SetLineColor(kRed);
		gr_rhoAvg->Draw("al");

		TMarker* m_rhoMin = new TMarker(tau,bestrhoavg,kCircle);
		m_rhoMin->Draw();
		c_rhoAvg->SaveAs("2D_minimizeRho_" + acceptanceName + ".pdf");
		

	  
		//Re-wrap 1D histograms into 2D
		rewrap1dhisto(hData_unfolded_unwrapped, hData_unfolded);
		// rewrap1dhisto(hData_bkgSub, hData_bkgSub_rewrapped);
		// rewrap1dhisto(hData_unwrapped, hData_unwrapped_rewrapped);


		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////// 4. Output a bunch of histograms and tables //////////////////////////////////
		float rmargin = gStyle->GetPadRightMargin();
		gStyle->SetPadRightMargin(0.17);

		// Data distribution (2D and unwrapped)
        TCanvas *c_data = new TCanvas("c_data", "c_data", 675, 600);
        gStyle->SetPalette(1);
        hData_bkgSub_combined->SetTitle("Data (background subtracted);"+xaxislabel+";"+yaxislabel);
        hData_bkgSub_combined->Draw("COLZ");
		hData_bkgSub_combined->GetZaxis()->SetMoreLogLabels();
        c_data->SetLogz();
        c_data->SaveAs("2D_data_" + acceptanceName +"_" + Var2D + ".pdf");
		hData_bkgSub->SetTitle("Unwrapped Data (background subtracted);Bin number;Entries per bin");
		hData_bkgSub->Draw();
        c_data->SaveAs("2D_dataunwrapped_" + acceptanceName +"_" + Var2D + ".pdf");

		// ttDil MC (2D and unwrapped)
        hTrue->SetTitle("True MC;"+xaxislabel+";"+yaxislabel);
		hTrue->GetZaxis()->SetMoreLogLabels();
		hTrue->Draw("COLZ");
        c_data->SaveAs("2D_true_" + acceptanceName +"_" + Var2D + ".pdf");
		hTrue_unwrapped->SetTitle("Unwrapped Truth;Bin number;Entries per bin");
		hTrue_unwrapped->Draw();
        c_data->SaveAs("2D_trueunwrapped_" + acceptanceName +"_" + Var2D + ".pdf");

		// Unfolded data (2D)
		hData_unfolded->SetTitle("Unfolded Data;"+xaxislabel+";"+yaxislabel);
		hData_unfolded->GetZaxis()->SetMoreLogLabels();
		hData_unfolded->Draw("COLZ");
        c_data->SaveAs("2D_unfolded_" + acceptanceName +"_" + Var2D + ".pdf");

		// Purity and stability (2D)
        TCanvas *c_purstab = new TCanvas("c_purstab", "c_purstab", 650, 600);
		gStyle->SetPadRightMargin(0.15);
        hPurity->SetTitle("Purity;"+xaxislabel+";"+yaxislabel);
        hStability->SetTitle("Stability;"+xaxislabel+";"+yaxislabel);
        hPurity->Draw("COLZ");
        c_purstab->SaveAs("2D_purity_" + acceptanceName +"_" + Var2D + ".pdf");
        hStability->Draw("COLZ");
        c_purstab->SaveAs("2D_stability_" + acceptanceName +"_" + Var2D + ".pdf");

		// Data_bkgsub, data unfolded, true top
		TCanvas *c1 = new TCanvas("c1","c1",1300,400);
		c1->Divide(3,1);
		c1->cd(1);
		hData_bkgSub_combined->Draw("COLZ");
		c1->cd(2);
		hData_unfolded->Draw("COLZ");
		c1->cd(3);
		hTrue->Draw("COLZ");
		c1->SaveAs("2D_data_comparison_"+acceptanceName+"_"+Var2D+".pdf");

		// Migration matrix
        TCanvas *c_resp = new TCanvas("c_resp", "c_resp", 1850, 600);
        TH2D *hResp = (TH2D*) hTrue_vs_Meas->Clone("response");
        gStyle->SetPalette(1);
		hResp->SetTitle("Migration matrix");
        hResp->GetXaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (reco)");
        hResp->GetYaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (gen)");
        hResp->Draw("COLZ");
        c_resp->SetLogz();
        c_resp->SaveAs("2D_Response_" + acceptanceName + "_" + Var2D + ".pdf");

		gStyle->SetPadRightMargin(rmargin);


		// A few lingering acceptance corrections ///////////////////////////////////////////////////
        for (Int_t x = 1; x <= nbinsx_gen; x++) {
		  for (Int_t y = 1; y<= nbinsy2D; y++) {
			if (acceptM[3]->GetBinContent(x,y) != 0) {
			  hTrue->SetBinContent(x,y, hTrue->GetBinContent(x,y) * 1.0 / acceptM[3]->GetBinContent(x,y));
			  hTrue->SetBinError  (x,y, hTrue->GetBinError(x,y)  * 1.0 / acceptM[3]->GetBinContent(x,y));
			}
		  }
		}

        TH2D *denomM_2d = (TH2D*) file->Get("denominator_" + acceptanceName + "_" + Var2D + "_all");


        //==================================================================
        //============== Print the asymmetry ===============================
		cout << "========= Variable: " << acceptanceName << "===================\n";
        Float_t Afb, AfbErr;

		cout << "Automatic tau = " << tau << endl;
		cout << "Minimum rhoAvg = " << bestrhoavg << endl;
		cout << "Bias scale = " << scaleBias << endl;

        GetAfb(hData_combined, Afb, AfbErr);
        cout << " Data: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hTrue, Afb, AfbErr);
        cout << " True Top: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hData_unfolded, Afb, AfbErr); // Formerly GetCorrectedAfb
        cout << " Unfolded: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " Unfolded: " << Afb << " +/-  " << AfbErr << endl;

        GetAfb(denomM_2d, Afb, AfbErr);
        cout << " True Top from acceptance denominator: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " True_Top_from_acceptance_denominator: " << Afb << " +/-  " << AfbErr << "\n";

        vector<double> afb_m;
        vector<double> afb_merr;
        vector<double> afb_m_denom;
        vector<double> afb_merr_denom;

		cout << "From unfolded data:" << endl;
        GetAvsY2d(hData_unfolded, afb_m, afb_merr, second_output_file);

        cout << " With corrected uncertainty: " << endl;  //this function fills the inclusive asymmetry at array index 0, and then the asym in each y bin
        GetCorrectedAfb2d(hData_unfolded, m_correctE, afb_m, afb_merr, second_output_file);

		cout << "From acceptance denominator:" << endl;
		GetAvsY2d(denomM_2d, afb_m_denom, afb_merr_denom, second_output_file);

		/*
        if (draw_truth_before_pT_reweighting)
        {
            GetAfb(denomM_nopTreweighting_0, AfbG[0], AfbErr);
            GetAfb(denomM_nopTreweighting_1, AfbG[1], AfbErr);
            GetAfb(denomM_nopTreweighting_2, AfbG[2], AfbErr);
        }
		*/

        TCanvas *c_afb = new TCanvas("c_afb", "c_afb", 500, 500);
        double ybinsForHisto[4] = {ybins2D[0], ybins2D[1], ybins2D[2], ybins2D[3]};
        if (Var2D == "mtt") ybinsForHisto[0] = 300.0;
        TH1D *hAfbVsMtt = new TH1D ("AfbVsMtt",  "AfbVsMtt",  3, ybinsForHisto);
        TH1D *hAfbVsMtt_statonly = new TH1D ("AfbVsMtt_statonly",  "AfbVsMtt_statonly",  3, ybinsForHisto);
        for (int nb = 0; nb < 3; nb++)
        {
            hAfbVsMtt->SetBinContent(nb + 1, afb_m[nb+1]);
            if (checkErrors)
            {
                if (includeSys)
                {
                    cout << "Difference between calculated and hard-coded stat errors: " << afb_merr[nb+1] - stat_corr[nb] << endl;
                }
                else
                {
                    cout << "Difference between calculated and hard-coded stat errors: " << afb_merr[nb+1] - stat_uncorr[nb] << endl;
                }
            }
            hAfbVsMtt->SetBinError(nb + 1,  sqrt( pow(afb_merr[nb+1], 2) + pow(syst_corr[nb], 2) ) );
            hAfbVsMtt_statonly->SetBinContent(nb + 1, afb_m[nb+1]);
            hAfbVsMtt_statonly->SetBinError(nb + 1, afb_merr[nb+1]);
        }

        //  GetAvsY(hTrue, m_unfoldE, afb_m, afb_merr);

        TH1D *hTop_AfbVsMtt = new TH1D ("Top_AfbVsMtt",  "Top_AfbVsMtt",  3, ybinsForHisto);
        for (int nb = 0; nb < 3; nb++)
        {
            hTop_AfbVsMtt->SetBinContent(nb + 1, afb_m_denom[nb]);
            hTop_AfbVsMtt->SetBinError(nb + 1, 0);
        }

        tdrStyle->SetErrorX(0.5);
        hAfbVsMtt->SetMinimum( hAfbVsMtt->GetMinimum() - 0.1 );
        hAfbVsMtt->SetMaximum( hAfbVsMtt->GetMaximum() + 0.1 );
        hAfbVsMtt->SetLineWidth( 2.0 );
        hAfbVsMtt->Draw("E");
        hAfbVsMtt_statonly->Draw("E1 same");
        hTop_AfbVsMtt->SetLineColor(TColor::GetColorDark(kRed));
        hTop_AfbVsMtt->SetMarkerColor(TColor::GetColorDark(kRed));
        hTop_AfbVsMtt->SetMarkerSize(0);
        hTop_AfbVsMtt->SetLineWidth( 2.0 );
        hAfbVsMtt->GetYaxis()->SetTitle(asymlabel);
        hAfbVsMtt->GetYaxis()->SetTitleOffset(1.2);
        hAfbVsMtt->GetXaxis()->SetTitle(yaxislabel + yaxisunit);
        hTop_AfbVsMtt->Draw("E same");

        TLegend* leg1 = new TLegend(0.6, 0.72, 0.9, 0.938, NULL, "brNDC");
        leg1->SetEntrySeparation(100);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.03);
        leg1->SetFillStyle(0);
        leg1->AddEntry(hAfbVsMtt, "data");
        leg1->AddEntry(hTop_AfbVsMtt,    "MC@NLO parton level");
        leg1->Draw();

        TPaveText *pt1 = new TPaveText(0.18, 0.88, 0.41, 0.91, "brNDC");
        pt1->SetName("pt1name");
        pt1->SetBorderSize(0);
        pt1->SetFillStyle(0);

        TText *blah;
        //blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
        blah = pt1->AddText("CMS, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
        blah->SetTextSize(0.032);
        blah->SetTextAlign(11);
        pt1->Draw();

        c_afb->SaveAs("2D_AfbVs" + Var2D + "_unfolded_" + acceptanceName + ".pdf");
        // c_afb->SaveAs("AfbVs" + Var2D + "_unfolded_" + acceptanceName + ".root");
        // c_afb->SaveAs("AfbVs" + Var2D + "_unfolded_" + acceptanceName + ".C");

		/*
        TCanvas *c_mttu = new TCanvas("c_mttu", "c_mttu", 500, 500);

        hData_unfolded->Scale(1. / hData_unfolded->Integral(), "width");
        hTrue->Scale(1. / hTrue->Integral(), "width");

        hData_unfolded->GetXaxis()->SetTitle(yaxislabel + " #times sign(" + xaxislabel + ")");
        hData_unfolded->GetYaxis()->SetTitle("d#sigma/d(" + yaxislabel + " #times sign(" + xaxislabel + "))");
        hData_unfolded->SetMinimum(0.0);
        hData_unfolded->SetMaximum( 2.0 * hData_unfolded->GetMaximum());
        hData_unfolded->SetMarkerStyle(23);
        hData_unfolded->SetMarkerSize(2.0);
        hData_unfolded->Draw("E");
        hData_unfolded->SetLineWidth(lineWidth);
        hTrue->SetLineWidth(lineWidth);
        hTrue->SetLineColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillColor(TColor::GetColorDark(kGreen));
        hTrue->Draw("hist same");
        hData_unfolded->Draw("E same");

        leg1 = new TLegend(0.6, 0.62, 0.9, 0.838, NULL, "brNDC");
        leg1->SetEntrySeparation(100);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetTextSize(0.03);
        leg1->SetFillStyle(0);
        leg1->AddEntry(hData_unfolded, "( Data-BG ) Unfolded");
        leg1->AddEntry(hTrue,    "MC@NLO parton level", "F");
        leg1->Draw();

        c_mttu->SaveAs(Var2D + "_2D_unfolded_" + acceptanceName + ".pdf");

        TFile *output = new TFile(Form("DataMC_%s.root", Var2D.Data()), "UPDATE");

        TH1D *hDataMCratio  = (TH1D *) hData_unfolded->Clone("hDataMCratio" + Var2D + acceptanceName);
        hDataMCratio->SetTitle("hDataMCratio" + Var2D + acceptanceName);
        hDataMCratio->Reset();
        hDataMCratio->Divide(hData_unfolded, hTrue, 1., 1.);
        hDataMCratio->Write();
		
        output->Close();
		*/
        ch_data->Delete();

        ch_top->Delete();

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
        {
            ch_bkg[iBkg]->Delete();
        }

    }
    myfile.close();
    second_output_file.close();
}

#ifndef __CINT__
int main ()
{
    AfbUnfoldExample();    // Main program when run stand-alone
    return 0;
}
#endif
