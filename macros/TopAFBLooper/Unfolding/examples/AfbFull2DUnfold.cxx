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


void AfbUnfoldExample(TString Var2D = "mtt", double scalettdil = 1., double scalefake = 2.18495, double scalewjets = 1., double scaleDYeemm = 1.35973, double scaleDYtautau = 1.17793, double scaletw = 1., double scaleVV = 1.)
{
  TH1::SetDefaultSumw2();

  setTDRStyle();
  gStyle->SetOptFit();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cout.precision(3);

  TString ChannelName[4] = {"diel", "dimu", "mueg", "all"};

  for( int iChan=0; iChan<4; iChan++ ) {

	TString channel_name;

	switch (iChan) {
	case 0:  channel_name = "ee";   break;
	case 1:  channel_name = "mm";   break;
	case 2:  channel_name = "em";   break;
	case 3:  channel_name = "all";  break;
	}

    TString summary_name;
    if (Var2D == "mtt") summary_name = "summary_2Dunfolding_" + channel_name ;
    else if (Var2D == "ttrapidity2") summary_name = "summary_2Dunfolding_ttrapidity2_" + channel_name;
    else if (Var2D == "ttpt") summary_name = "summary_2Dunfolding_ttpt_" + channel_name;

    TString yaxisunit;
    if (Var2D == "mtt") yaxisunit = " (GeV/c^{2})";
    else if (Var2D == "ttrapidity2") yaxisunit = "";
    else if (Var2D == "ttpt") yaxisunit = " (GeV/c)";

    // if (!(scalefake == 1. && scalewjets == 1. && scaleDY == 1. && scaletw == 1. && scaleVV == 1.))  summary_name = summary_name + Form("_%i_%i_%i_%i_%i", int(10.*scalefake + 0.5), int(10.*scalewjets + 0.5), int(10.*scaleDY + 0.5), int(10.*scaletw + 0.5), int(10.*scaleVV + 0.5));

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);

    ofstream myfile;
    myfile.open (summary_name + ".txt");
    cout.rdbuf(myfile.rdbuf());

    // OGU 130516: add second output txt file with format easier to be pasted into google docs
    ofstream second_output_file;
    second_output_file.open(summary_name + "_formated.txt");

    const int nBkg = 11;
	const int nSig = 3;
    TString path = "../";
    TString dataroot[nSig] = {"data_diel_baby.root", "data_dimu_baby.root", "data_mueg_baby.root"};
    TString bkgroot[nBkg];
	double bkgSF[nBkg];

	bkgroot[0] = "DY1to4Jeemm_baby.root";     	bkgSF[0] = scaleDYeemm;
	bkgroot[1] = "DY1to4Jtautau_baby.root";	 	bkgSF[1] = scaleDYtautau;
	bkgroot[2] = "diboson_baby.root";		 	bkgSF[2] = scaleVV;
	bkgroot[3] = "tW_lepdl_baby.root";		 	bkgSF[3] = scaletw;
	bkgroot[4] = "tW_lepfake_baby.root";		bkgSF[4] = scalefake;
	bkgroot[5] = "tW_lepsl_baby.root";		 	bkgSF[5] = scalefake;
	bkgroot[6] = "triboson_baby.root";		 	bkgSF[6] = scaleVV;
	bkgroot[7] = "ttV_baby.root";			 	bkgSF[7] = scaleVV;
	bkgroot[8] = "ttfake_mcatnlo_baby.root";	bkgSF[8] = scalefake;
	bkgroot[9] = "ttsl_mcatnlo_baby.root";	 	bkgSF[9] = scalefake;
	bkgroot[10] = "w1to4jets_baby.root";		bkgSF[10] = scalefake;

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
		int nbinsx_reco_split_ch = -99;
		int nbinsunwrapped_gen = -99;
		int nbinsunwrapped_reco = -99;

		if( iVar < 2 || iVar==9 ) nbinsx_gen = nbinsx2D*2;
		else nbinsx_gen = nbinsx2D;

		nbinsx_reco = nbinsx_gen*2;
		//nbinsx_reco = nbinsx_gen;
		nbinsx_reco_split_ch = nbinsx_reco*nSig;

		double* genbins;
		double* recobins;
		double* recobins_split_ch;

		genbins = new double[nbinsx_gen+1];
		recobins = new double[nbinsx_reco+1];
		recobins_split_ch = new double[nbinsx_reco_split_ch+1];

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

		// Make reco binning array including the channel split
		double recohist_width = fabs(recobins[nbinsx_reco] - recobins[0]);
		std::copy( recobins, recobins+nbinsx_reco+1, recobins_split_ch );
		for( int i=nbinsx_reco+1; i<nbinsx_reco_split_ch+1; i++ ) {
		  recobins_split_ch[i] = recobins_split_ch[i-nbinsx_reco] + recohist_width;
		}

		nbinsunwrapped_gen  = nbinsx_gen  * nbinsy2D;
		nbinsunwrapped_reco = nbinsx_reco * nbinsy2D;
		
		//Make histograms

		//Use 2D histograms to store the distributions
        TH2D *hData = new TH2D ("Data", "Data", nbinsx_reco, recobins, nbinsy2D, ybins2D);
        TH2D *hData_split = new TH2D ("Data_split", "Data_split", nbinsx_reco_split_ch, recobins_split_ch, nbinsy2D, ybins2D);
        TH2D *hBkg = new TH2D ("Background",  "Background", nbinsx_reco, recobins, nbinsy2D, ybins2D);
        TH2D *hBkg_split = new TH2D ("Background_split",  "Background_split", nbinsx_reco_split_ch, recobins_split_ch, nbinsy2D, ybins2D);
        TH2D *hData_unfolded = new TH2D ("Data_Unfold", "Data with background subtracted and unfolded", nbinsx_gen, genbins, nbinsy2D, ybins2D);
        TH2D *hTrue = new TH2D ("true", "Truth",    nbinsx_gen, genbins, nbinsy2D, ybins2D);
        TH2D *hTrue_split = new TH2D ("true_split", "Truth",    nbinsx_reco_split_ch, recobins_split_ch, nbinsy2D, ybins2D);
        TH2D *hMeas = new TH2D ("meas", "Measured", nbinsx_reco, recobins, nbinsy2D, ybins2D);
		TH2D *hData_bkgSub_rewrapped = new TH2D ("bkgsub", "bkgsub", nbinsx_reco, recobins, nbinsy2D, ybins2D);
		TH2D *hPurity = new TH2D("purity", "Purity", nbinsx_gen, genbins, nbinsy2D, ybins2D);
		TH2D *hStability = new TH2D("stability", "Stability", nbinsx_gen, genbins, nbinsy2D, ybins2D);

		//Unwrapped histograms have n bins (where n = nx*ny), centered around the integers from 1 to n.
		TH1D *hData_unwrapped = new TH1D ("Data_Unwr", "Unwrapped data with background subtracted", nbinsunwrapped_reco, 0.5, double(nbinsunwrapped_reco)+0.5);
        TH1D *hBkg_unwrapped = new TH1D ("Background_Unwr",  "Background unwrapped",    nbinsunwrapped_reco, 0.5, double(nbinsunwrapped_reco)+0.5);
        TH1D *hData_unfolded_unwrapped = new TH1D ("Data_Unfold_Unwr", "Data unfolded and unwrapped", nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);
        TH1D *hTrue_unwrapped = new TH1D ("true_unwr", "Truth unwrapped",  nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);
        TH1D *hMeas_unwrapped = new TH1D ("meas_unwr", "Measured unwrapped", nbinsunwrapped_reco, 0.5, double(nbinsunwrapped_reco)+0.5);

		//Migration matrix, using the unwrapped binning on both axes
        TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbinsunwrapped_reco, 0.5, double(nbinsunwrapped_reco)+0.5, nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);

        TH1D *hData_bkgSub;
        TH2D *hData_bkgSub_split;
		// TH1D* hMeas_newErr;

		hTrue_split->RebinX(2);

        TMatrixD m_smearingE(nbinsunwrapped_gen, nbinsunwrapped_gen);
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

        ch_top->Add(path + "ttdl_mcatnlo_baby.root");

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
		  {
            ch_bkg[iBkg] = new TChain("tree");
            ch_bkg[iBkg]->Add(path + bkgroot[iBkg]);
		  }

		double offset = 0;
		double histmax = recobins[nbinsx_reco];
		double histmin = recobins[0];
		double hiBinCenter = hData->GetXaxis()->GetBinCenter(nbinsx_reco);
		double loBinCenter = hData->GetXaxis()->GetBinCenter(1);

		delete[] genbins;
		delete[] recobins;
		delete[] recobins_split_ch;

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
			// If we're unfolding a single channel, skip events that aren't in that channel
			if( iChan < nSig && channel != iChan) continue;

			//Use an offset to sort events into 2 superbins: ee/mumu, emu
			//if( channel > 0 ) channel--;
			offset = double(channel) * recohist_width;
			//Do the same thing as "fillUnderOverflow", except adapted for 3x1 histograms
			if( observable > histmax )        observable = hiBinCenter;
			else if( observable < histmin )   observable = loBinCenter;
			if( observableMinus > histmax )        observableMinus = hiBinCenter;
			else if( observableMinus < histmin )   observableMinus = loBinCenter;

            if ( tmass > 0 )
			  {
				fillUnderOverFlow(hData, observable, obs2D, weight, Nsolns);
				fillUnderOverFlow(hData_split, observable+offset, obs2D, weight, Nsolns);
				if (combineLepMinus) {
				  fillUnderOverFlow(hData, observableMinus, obs2D, weight, Nsolns);
				  fillUnderOverFlow(hData_split, observableMinus+offset, obs2D, weight, Nsolns);
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
				if( iChan < nSig && channel != iChan) continue;

				offset = double(channel) * recohist_width;
				obs2D = fabs(obs2D);
                weight *= bkgSF[iBkg];

				if( observable > histmax )        observable = hiBinCenter;
				else if( observable < histmin )   observable = loBinCenter;
				if( observableMinus > histmax )        observableMinus = hiBinCenter;
				else if( observableMinus < histmin )   observableMinus = loBinCenter;

                if ( tmass > 0 )
				  {
					fillUnderOverFlow(hBkg, observable, obs2D, weight, Nsolns);
					fillUnderOverFlow(hBkg_split, observable+offset, obs2D, weight, Nsolns);
					if (combineLepMinus) fillUnderOverFlow(hBkg, observableMinus, obs2D, weight, Nsolns);
					if (combineLepMinus) fillUnderOverFlow(hBkg_split, observableMinus+offset, obs2D, weight, Nsolns);
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

		int measbin = -99;
		int genbin = -99;

        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
		  {
            ch_top->GetEntry(i);
			if( iChan < nSig && channel != iChan) continue;
            obs2D = fabs(obs2D);
            obs2D_gen = fabs(obs2D_gen);
            weight *= scalettdil;
			//if( channel>0 ) channel--;

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
				measbin = getUnwrappedBin(hData, observable, obs2D);

				// cout << "Channel: " << channel << ". Offset: " << offset << endl;

				fillUnderOverFlow(hMeas, observable, obs2D, weight, Nsolns);
				fillUnderOverFlow(hTrue, observable_gen, obs2D_gen, weight, Nsolns);
				fillUnderOverFlow(hTrue_split, observable_gen+offset, obs2D_gen, weight, Nsolns);
				fillUnderOverFlow(hTrue_vs_Meas, measbin, genbin, weight, Nsolns);
				if ( combineLepMinus )
				  {
					genbin  = getUnwrappedBin(hTrue, observableMinus_gen, obs2D_gen);
					measbin = getUnwrappedBin(hData, observableMinus, obs2D);

					fillUnderOverFlow(hMeas, observableMinus, obs2D, weight, Nsolns);
					fillUnderOverFlow(hTrue, observableMinus_gen, obs2D_gen, weight, Nsolns);
					fillUnderOverFlow(hTrue_split, observableMinus_gen+offset, obs2D_gen, weight, Nsolns);
					fillUnderOverFlow(hTrue_vs_Meas, measbin, genbin, weight, Nsolns);
				  }
			  }
		  }

		  hData_bkgSub_split = (TH2D *) hData_split->Clone("Data_BkgSub_split");
		  hData_bkgSub_split->Add(hBkg_split, -1.0);


		// Do the acceptance correction, by filling the migration matrix with events that have a gen-level value but no reco-level value
        TFile *file = new TFile("../denominator/acceptance/mcnlo/accept_" + acceptanceName + ".root");

		TH2D *acceptM[4];
		acceptM[0] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_diel"));
		acceptM[1] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_dimu"));
		acceptM[2] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_mueg"));
		acceptM[3] = (TH2D*)(file->Get("accept_" + acceptanceName + "_" + Var2D + "_all" ));

		TH2D *accNum[4];
		accNum[0] = (TH2D*)(file->Get("numerator_" + acceptanceName + "_" + Var2D +"_diel"));
		accNum[1] = (TH2D*)(file->Get("numerator_" + acceptanceName + "_" + Var2D +"_dimu"));
		accNum[2] = (TH2D*)(file->Get("numerator_" + acceptanceName + "_" + Var2D +"_mueg"));
		accNum[3] = (TH2D*)(file->Get("numerator_" + acceptanceName + "_" + Var2D +"_all"));

		TH2D *accDen[4];
		accDen[0] = (TH2D*)(file->Get("denominator_" + acceptanceName + "_" + Var2D +"_diel"));
		accDen[1] = (TH2D*)(file->Get("denominator_" + acceptanceName + "_" + Var2D +"_dimu"));
		accDen[2] = (TH2D*)(file->Get("denominator_" + acceptanceName + "_" + Var2D +"_mueg"));
		accDen[3] = (TH2D*)(file->Get("denominator_" + acceptanceName + "_" + Var2D +"_all"));

		//Tricks to make "channel 0" hold the combined same-flavor histograms
		//accNum[0]->Add( accNum[1] );
		//accDen[0]->Add( accDen[1] );
		//acceptM[0] = (TH2D*)(accNum[0]->Clone("accept_SF"));
		//acceptM[0]->Divide( accDen[0] );

		if( iChan==nSig ) {     // Combined channels...

		  double gen_integrals[4] = {0.}; 
		  double reco_integrals[4] = {0.};
  
		  //adjust A matrix to match the channel proportions in data

		  //Figure out the relative proportion of events in each channel
		  for( int aChannel=0; aChannel<nSig; aChannel++ ) {
		    gen_integrals[aChannel] = hTrue_split->Integral( aChannel*nbinsx_gen+1, (aChannel+1)*nbinsx_gen, 1, nbinsy2D );
		    reco_integrals[aChannel] = hData_bkgSub_split->Integral( aChannel*nbinsx_reco+1, (aChannel+1)*nbinsx_reco, 1, nbinsy2D );
		    gen_integrals[nSig] += gen_integrals[aChannel];
		    reco_integrals[nSig] += reco_integrals[aChannel];
		  }



		  TH2D* acceptNumcorrected = (TH2D*)(accNum[nSig]->Clone("acceptNumcorrected"));
		  TH2D* acceptDencorrected = (TH2D*)(accDen[nSig]->Clone("acceptDencorrected"));

		  double correction[nSig];
		  acceptNumcorrected->Reset();
		  acceptDencorrected->Reset();
		  //acceptNumcorrected->Print("all");

		  for( int aChannel=0; aChannel<nSig; aChannel++ ) {
			correction[aChannel] = (reco_integrals[aChannel] / reco_integrals[nSig]) / (gen_integrals[aChannel] / gen_integrals[nSig]);
			cout<<"Acceptance channel contribution factor for channel "<<ChannelName[aChannel]<<": "<<correction[aChannel]<<endl;
			acceptNumcorrected->Add(accNum[aChannel], correction[aChannel]);
			acceptDencorrected->Add(accDen[aChannel], correction[aChannel]);
			//acceptNumcorrected->Print("all");
		  }

		  TH2D* acceptMcorrected = (TH2D*)(acceptNumcorrected->Clone("acceptMcorrected"));
		  acceptMcorrected->Divide( acceptDencorrected );

		  // to be fully correct, we should really reweight the events by correction[aChannel] when we fill hTrue and hTrue_vs_Meas too
		  for( int yBin=1; yBin<=nbinsy2D; yBin++ ) {
			for( int xBin=1; xBin<=nbinsx_gen; xBin++ ) {
			  genbin = (yBin-1)*nbinsx_gen + xBin;
			  double acceptance = acceptMcorrected->GetBinContent( xBin, yBin );
			  double n_accepted = hTrue->GetBinContent( xBin, yBin );
			  double n_rejected = n_accepted/acceptance - n_accepted;
			  hTrue_vs_Meas->SetBinContent( 0, genbin, n_rejected );
			  double num_error = accNum[iChan]->GetBinError( xBin, yBin );
			  double den_error = accDen[iChan]->GetBinError( xBin, yBin );
			  double new_error = sqrt( den_error*den_error - num_error*num_error );
			  hTrue_vs_Meas->SetBinError(   0, genbin, new_error );
			}
		  }


/*
		  for( int aChannel=0; aChannel<2; aChannel++ ) {
			for( int yBin=1; yBin<=nbinsy2D; yBin++ ) {
			  for( int xBin=1; xBin<=nbinsx_gen; xBin++ ) {

				genbin = (yBin-1)*nbinsx_gen + xBin;
				double old_error = hTrue_vs_Meas->GetBinError( 0, genbin );

				//Calculate the number of rejected events, and fill it into the smearing matrix underflow
				double acceptance = acceptM[aChannel+1]->GetBinContent( xBin, yBin );
				double n_accepted = hTrue_split->GetBinContent( aChannel*nbinsx_gen + xBin, yBin );
				double n_rejected = n_accepted/acceptance - n_accepted;
				double correction = (reco_integrals[aChannel] / reco_integrals[2]) / (gen_integrals[aChannel] / gen_integrals[2]);
				hTrue_vs_Meas->Fill( -999999, double(genbin), n_rejected*correction );

				//Calculate the uncertainty on the rejected events
				double num_error = accNum[aChannel]->GetBinError( xBin, yBin );
				double den_error = accDen[aChannel]->GetBinError( xBin, yBin );
				double new_error = sqrt( old_error*old_error + den_error*den_error - num_error*num_error );
				hTrue_vs_Meas->SetBinError( 0, genbin, new_error );

			  }
			}
		  }
*/

		}
		else {     // Individual channels...
		  for( int yBin=1; yBin<=nbinsy2D; yBin++ ) {
			for( int xBin=1; xBin<=nbinsx_gen; xBin++ ) {
			  genbin = (yBin-1)*nbinsx_gen + xBin;
			  double acceptance = acceptM[iChan]->GetBinContent( xBin, yBin );
			  double n_accepted = hTrue->GetBinContent( xBin, yBin );
			  double n_rejected = n_accepted/acceptance - n_accepted;
			  hTrue_vs_Meas->SetBinContent( 0, genbin, n_rejected );
			  double num_error = accNum[iChan]->GetBinError( xBin, yBin );
			  double den_error = accDen[iChan]->GetBinError( xBin, yBin );
			  double new_error = sqrt( den_error*den_error - num_error*num_error );
			  hTrue_vs_Meas->SetBinError(   0, genbin, new_error );
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
		else if( hData->GetNbinsX() == 2*nbinsx_gen ) {

		  for( int by=1; by<=nbinsy2D; by++ ) {
			for( int bx=1; bx<=nbinsx_gen; bx++ ) {

			  int outbin = (by-1)*nbinsx_gen + bx;
			  double numerator = hTrue_vs_Meas->GetBinContent(outbin*2, outbin) + hTrue_vs_Meas->GetBinContent(outbin*2-1, outbin);
			  double purity_denom = hMeas->GetBinContent( bx*2, by ) + hMeas->GetBinContent( bx*2-1, by );
			  double stability_denom = hTrue->GetBinContent( bx, by );
 
			  hPurity->SetBinContent( bx, by, numerator / purity_denom );
			  hStability->SetBinContent( bx, by, numerator / stability_denom );
			}
		  }
		}
		else {
		  cout << "\n***WARNING: Purity and stability plots are broken, because nbinsx_gen != nbinsx_reco!!!\n" << endl;
		}

		///////////////////////////////////////////////////////////////////////
		// Set data-like stat errors on MC for optimizing tau

		unwrap2dhisto(hTrue, hTrue_unwrapped);
		unwrap2dhisto(hMeas, hMeas_unwrapped);
		unwrap2dhisto(hBkg,  hBkg_unwrapped);

		//Set data-like stat errors on MC for optimizing tau
		for( int i=1; i<=nbinsunwrapped_reco; i++) {
		  double n_sig = hMeas_unwrapped->GetBinContent(i);
		  double n_bkg = hBkg_unwrapped->GetBinContent(i);
		  double bkg_err = hBkg_unwrapped->GetBinError(i);
		  hMeas_unwrapped->SetBinError(i, sqrt(n_sig + n_bkg + bkg_err*bkg_err ) );
		}

		double tempScaleBias = hMeas_unwrapped->Integral() / hTrue_unwrapped->Integral();

		TUnfoldSys unfold_FindTau (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone, TUnfold::kEConstraintArea);
		unfold_FindTau.SetInput(hMeas_unwrapped);
		unfold_FindTau.RegularizeBins2D(1,1,nbinsx_gen,nbinsx_gen,nbinsy2D,TUnfold::kRegModeCurvature);
		minimizeRhoAverage(&unfold_FindTau, hMeas_unwrapped, -6.0, -1.0);
		tau = unfold_FindTau.GetTau();

		// Generate a curve of rhoAvg vs log(tau)
		double ar_tau[100];
		double ar_rhoAvg[100];
		double tau_test = 0.0;
		double bestrhoavg = unfold_FindTau.GetRhoAvg();

		for(int l=0; l<100; l++) {
		  tau_test = pow( 10, -6.0 + 0.05*l );
		  unfold_FindTau.DoUnfold(tau_test, hMeas_unwrapped, tempScaleBias);
		  ar_tau[l] = tau_test;
		  ar_rhoAvg[l] = unfold_FindTau.GetRhoAvg();
		}

		TGraph* gr_rhoAvg = new TGraph(100, ar_tau, ar_rhoAvg);
		TCanvas* c_rhoAvg = new TCanvas("c_rhoAvg","c_rhoAvg");
		c_rhoAvg->SetLogx();
		gr_rhoAvg->SetTitle("Global Correlation Coefficient;#tau;#rho_{avg}");
		gr_rhoAvg->SetLineColor(kRed);
		gr_rhoAvg->Draw("al");

		TMarker* m_rhoMin = new TMarker(tau,bestrhoavg,kCircle);
		m_rhoMin->Draw();
		c_rhoAvg->SaveAs("2D_" + acceptanceName + "_" + Var2D + "_" + channel_name + "_unfoldTests_minRho.pdf");


		////////////////////////////////////////////////////////////////////////////////////////////////
		//////////////////////// 3. Perform the unfolding procedure ////////////////////////////////////

		// First step: unwrap our TH2Ds into TH1Ds!
		unwrap2dhisto(hData, hData_unwrapped);
		unwrap2dhisto(hBkg,  hBkg_unwrapped);
		unwrap2dhisto(hTrue, hTrue_unwrapped);
		unwrap2dhisto(hMeas, hMeas_unwrapped);

        hData_bkgSub = (TH1D *) hData_unwrapped->Clone("data_bkgsub");
        hData_bkgSub->Add(hBkg_unwrapped, -1.0);

		// Now let's actually do the unfolding.
		TUnfoldSys unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone, TUnfold::kEConstraintArea);  //need to set reg mode "None" here if regularizing by hand
		unfold_TUnfold.SetInput(hData_bkgSub);
		scaleBias = (hData->Integral() - hBkg->Integral()) / hTrue_unwrapped->Integral();
		hTrue_unwrapped->Scale(scaleBias);
		//unfold_TUnfold.SetBias(hTrue_unwrapped);
		unfold_TUnfold.RegularizeBins2D(1,1,nbinsx_gen,nbinsx_gen,nbinsy2D,TUnfold::kRegModeCurvature);
		// minimizeRhoAverage(&unfold_TUnfold, hData_bkgSub, -5.0, -0.5);
		unfold_TUnfold.DoUnfold(tau, hData_bkgSub, scaleBias);
		unfold_TUnfold.GetOutput(hData_unfolded_unwrapped);
		// tau = unfold_TUnfold.GetTau();

		TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
		TH2D *ematrix_smearing = (TH2D*)ematrix->Clone("ematrix_smearing");
		unfold_TUnfold.GetEmatrixSysUncorr( ematrix_smearing, 0, false );

		for (Int_t cmi = 0; cmi < nbinsunwrapped_gen; cmi++)
		  {
			for (Int_t cmj = 0; cmj < nbinsunwrapped_gen; cmj++)
			  {
				m_smearingE(cmi, cmj) = ematrix_smearing->GetBinContent(cmi + 1, cmj + 1);
				m_unfoldE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
				m_correctE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
			  }
		  }


		//Re-wrap 1D histograms into 2D
		rewrap1dhisto(hData_unfolded_unwrapped, hData_unfolded);
		rewrap1dhisto(hData_bkgSub, hData_bkgSub_rewrapped);
		// rewrap1dhisto(hData_unwrapped, hData_unwrapped_rewrapped);


		/////////////////////////////////////////////////////////////////////////////////////////////
		/////////////// 4. Output a bunch of histograms and tables //////////////////////////////////
		float rmargin = gStyle->GetPadRightMargin();
		gStyle->SetPadRightMargin(0.17);

		// Data distribution (2D and unwrapped)
        TCanvas *c_data = new TCanvas("c_data", "c_data", 675, 600);
        gStyle->SetPalette(1);
        hData_bkgSub_rewrapped->SetTitle("Data (background subtracted);"+xaxislabel+";"+yaxislabel);
        hData_bkgSub_rewrapped->Draw("COLZ");
		hData_bkgSub_rewrapped->GetZaxis()->SetMoreLogLabels();
        c_data->SetLogz();
        c_data->SaveAs("2D_data_" + acceptanceName + "_" + Var2D + "_" + channel_name + ".pdf");
		hData_bkgSub->SetTitle("Unwrapped Data (background subtracted);Bin number;Entries per bin");
		hData_bkgSub->Draw();
        c_data->SaveAs("2D_dataunwrapped_" + acceptanceName +"_" + Var2D + "_" + channel_name + ".pdf");

		// ttDil MC (2D and unwrapped)
        hTrue->SetTitle("True MC;"+xaxislabel+";"+yaxislabel);
		hTrue->GetZaxis()->SetMoreLogLabels();
		hTrue->Draw("COLZ");
        c_data->SaveAs("2D_true_" + acceptanceName +"_" + Var2D + "_" + channel_name + ".pdf");
		hTrue_unwrapped->SetTitle("Unwrapped Truth;Bin number;Entries per bin");
		hTrue_unwrapped->Draw();
        c_data->SaveAs("2D_trueunwrapped_" + acceptanceName +"_" + Var2D + "_" + channel_name + ".pdf");

		// Unfolded data (2D)
		hData_unfolded->SetTitle("Unfolded Data;"+xaxislabel+";"+yaxislabel);
		hData_unfolded->GetZaxis()->SetMoreLogLabels();
		hData_unfolded->Draw("COLZ");
        c_data->SaveAs("2D_unfolded_" + acceptanceName +"_" + Var2D + "_" + channel_name + ".pdf");

		// Purity and stability (2D)
        TCanvas *c_purstab = new TCanvas("c_purstab", "c_purstab", 650, 600);
		gStyle->SetPadRightMargin(0.15);
        hPurity->SetTitle("Purity;"+xaxislabel+";"+yaxislabel);
        hStability->SetTitle("Stability;"+xaxislabel+";"+yaxislabel);
        hPurity->Draw("COLZ");
        c_purstab->SaveAs("2D_purity_" + acceptanceName +"_" + Var2D + "_" + channel_name + ".pdf");
        hStability->Draw("COLZ");
        c_purstab->SaveAs("2D_stability_" + acceptanceName +"_" + Var2D + "_" + channel_name + ".pdf");

		// Data_bkgsub, data unfolded, true top
		TCanvas *c1 = new TCanvas("c1","c1",1300,400);
		c1->Divide(3,1);
		c1->cd(1);
		hData_bkgSub_rewrapped->Draw("COLZ");
		c1->cd(2);
		hData_unfolded->Draw("COLZ");
		c1->cd(3);
		hTrue->Draw("COLZ");
		c1->SaveAs("2D_data_comparison_"+acceptanceName+"_"+Var2D+"_" + channel_name + ".pdf");

		// Migration matrix
        TCanvas *c_resp = new TCanvas("c_resp", "c_resp", 700, 600);
        TH2D *hResp = (TH2D*) hTrue_vs_Meas->Clone("response");
        gStyle->SetPalette(1);
		hResp->SetTitle("Migration matrix");
        hResp->GetXaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (reco)");
        hResp->GetYaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (gen)");
        hResp->Draw("COLZ");
        c_resp->SetLogz();
        c_resp->SaveAs("2D_Response_" + acceptanceName + "_" + Var2D + "_" + channel_name + ".pdf");

		gStyle->SetPadRightMargin(rmargin);


		// A few lingering acceptance corrections ///////////////////////////////////////////////////
		//acceptM[iChan]->Scale(1.0 / acceptM[iChan]->Integral());
        for (Int_t x = 1; x <= nbinsx_gen; x++) {
		  for (Int_t y = 1; y<= nbinsy2D; y++) {
			if (acceptM[iChan]->GetBinContent(x,y) != 0) {
			  hTrue->SetBinContent(x,y, hTrue->GetBinContent(x,y) * 1.0 / acceptM[iChan]->GetBinContent(x,y));
			  hTrue->SetBinError  (x,y, hTrue->GetBinError(x,y)  * 1.0 / acceptM[iChan]->GetBinContent(x,y));
			}
		  }
		}

        TH2D *denomM_2d = (TH2D*) file->Get("denominator_" + acceptanceName + "_" + Var2D + "_" + ChannelName[iChan]);


        //==================================================================
        //============== Print the asymmetry ===============================
		cout << "========= Variable: " << acceptanceName << "===================\n";
        Float_t Afb, AfbErr;

        vector<double> afb_m;
        vector<double> afb_merr;
        vector<double> afb_m_denom;
        vector<double> afb_merr_denom;

		cout << "Automatic tau = " << tau << endl;
		cout << "Minimum rhoAvg = " << bestrhoavg << endl;
		cout << "Bias scale = " << scaleBias << endl;

        GetAfb(hData, Afb, AfbErr);
        cout << " Data: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hTrue, Afb, AfbErr);
        cout << " True Top: " << Afb << " +/-  " << AfbErr << "\n";

		GetAfb(denomM_2d, Afb, AfbErr);
        cout << " True Top from acceptance denominator: " << Afb << " +/-  " << AfbErr << endl;
        second_output_file << acceptanceName << " " << observablename << " True_Top_from_acceptance_denominator: " << Afb << " +/-  " << AfbErr << "\n";

		GetAvsY2d(denomM_2d, afb_m_denom, afb_merr_denom, second_output_file);
		for( int i=0; i<nbinsy2D; i++ ) {
	        cout << " True Top from acceptance denominator bin" << i << ": " << afb_m_denom.at(i) << " +/-  " << afb_merr_denom.at(i) << endl;
	        second_output_file << acceptanceName << " " << observablename << " True_Top_from_acceptance_denominator_bin"<<i<<": "<< afb_m_denom.at(i) << " +/-  " << afb_merr_denom.at(i) << "\n";
	    }

		GetCorrectedAfb2d(hData_unfolded, m_smearingE, afb_m, afb_merr, second_output_file);
        cout << " Unfolded with smearing errors: " << afb_m.at(0) << " +/- " << afb_merr.at(0) << endl;

		GetCorrectedAfb2d(hData_unfolded, m_correctE, afb_m, afb_merr, second_output_file);
        cout << " Unfolded: " << afb_m.at(0) << " +/- " << afb_merr.at(0) << endl;
        second_output_file << acceptanceName << " " << observablename << " Unfolded: " << afb_m.at(0) << " +/- " << afb_merr.at(0) << endl;

		// Asymmetries in each y-bin
		for( int i=1; i<=nbinsy2D; i++ ) {
		  cout << Var2D << " bin" << i << ": " << afb_m.at(i) << " +/- " << afb_merr.at(i) << endl;
		  second_output_file << acceptanceName << " " << observablename << " " << Var2D << "bin" << i << ": " << afb_m.at(i) << " +/- " << afb_merr.at(i) << endl;
		}

		for( int j=0; j<nbinsy2D; j++ ) {
			for( int i=0; i<nbinsx_gen/2; i++ ) {
				cout << Var2D << " double differential bin (" << j + 1 << "," << i + 1 << "): " << afb_m.at(4 + j*nbinsx_gen/2 + i) << " +/- " << afb_merr.at(4 + j*nbinsx_gen/2 + i) << endl;
				second_output_file << acceptanceName << " " << observablename << " " << Var2D << "DDbin" << j + 1 << "x" << i + 1 << ": " << afb_m.at(4 + j*nbinsx_gen/2 + i) << " +/- " << afb_merr.at(4 + j*nbinsx_gen/2 + i) << endl;
			}
		}

		// Contents in each bin
		hData_unfolded->Scale( 1. / hData_unfolded->Integral() );

		for( int j=1; j<=nbinsy2D; j++ ) {
		  for( int i=1; i<=nbinsx_gen; i++ ) {
			cout << "Bin (" << j << "," << i << "): " << hData_unfolded->GetBinContent(i,j) << " +/- " << hData_unfolded->GetBinError(i,j) << endl;
			second_output_file << acceptanceName << " " << observablename << " bin" << j << "x" << i << ": " << hData_unfolded->GetBinContent(i,j) << " +/- " << hData_unfolded->GetBinError(i,j) << endl;
		  }
		}

		// Pairs: print out asym +/- err for each bin pair




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
		double minmin = min( hAfbVsMtt->GetMinimum(), hTop_AfbVsMtt->GetMinimum() );
		double maxmax = max( hAfbVsMtt->GetMaximum() + hAfbVsMtt->GetBinError(hAfbVsMtt->GetMaximumBin()), hTop_AfbVsMtt->GetMaximum() + hTop_AfbVsMtt->GetBinError(hTop_AfbVsMtt->GetMaximumBin()));
		double spread = maxmax - minmin;
		if( spread > 0.25 ) maxmax += 0.25*spread;
        hAfbVsMtt->SetMinimum( minmin - 0.1 );
        hAfbVsMtt->SetMaximum( maxmax + 0.1 );
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
        blah = pt1->AddText("CMS, 19.5 fb^{-1} at  #sqrt{s}=8 TeV");
        blah->SetTextSize(0.032);
        blah->SetTextAlign(11);
        pt1->Draw();

        c_afb->SaveAs("2D_AfbVs" + Var2D + "_unfolded_" + acceptanceName + "_" + channel_name + ".pdf");
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
  } // End loop over channels
}

#ifndef __CINT__
int main ()
{
  AfbUnfoldExample();    // Main program when run stand-alone
  return 0;
}
#endif
