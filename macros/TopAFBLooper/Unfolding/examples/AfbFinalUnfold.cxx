#include <iostream>
#include <fstream>

#include "TROOT.h"
#include "TRandom3.h"
#include "TH1D.h"
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
#include "TPaveText.h"
#include "TLatex.h"
#include "TMarker.h"
#include "TProfile.h"
#include "TF1.h"

#include "TUnfold.h"
#include "TUnfoldSys.h"

#include "AfbFinalUnfold.h"
#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================

Int_t kterm = 3; //for SVD
Double_t tau = 0.005; //for TUnfold - this is a more reasonable default (1E-4 gives very little regularisation)
Int_t nVars = 12;
Int_t includeSys = 0;
bool drawDiffs = true;
int lineWidthDiffs=drawDiffs?lineWidth*2/3:lineWidth;
bool checkErrors = false; //turn this on when making the final plots for the paper, to check the hard-coded systematics have been correctly entered
bool draw_truth_before_pT_reweighting = true; //turn this on when making the final plots for the paper (want to compare the data against the unweighted MC)
//bool drawTheory = true; //turn this on to show Bernreuther's predictions for AdeltaPhi and Ac1c2
bool maketauplots = false;


void AfbUnfoldExample(double scalettdil = 1., double scalefake = 2.18495, double scalewjets = 1., double scaleDYeemm = 1.35973, double scaleDYtautau = 1.17793, double scaletw = 1., double scaleVV = 1. )
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

    TString summary_name = "summary_1Dunfolding_" + channel_name;

    // if (!(scalefake == 1. && scalewjets == 1. && scaleDY == 1. && scaletw == 1. && scaleVV == 1.))  summary_name = Form("summary_1Dunfolding_%i_%i_%i_%i_%i", int(10.*scalefake + 0.5), int(10.*scalewjets + 0.5), int(10.*scaleDY + 0.5), int(10.*scaletw + 0.5), int(10.*scaleVV + 0.5));

    ofstream myfile;
    myfile.open (summary_name + ".txt");
	cout.rdbuf(myfile.rdbuf());

    // OGU 130516: add second output txt file with format easier to be pasted into google docs
    ofstream second_output_file;
    second_output_file.open(summary_name + "_formated.txt");

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);


    const int nBkg = 11;
    const int nSig = 3;
    TString path = "../";
    TString dataroot[nSig] = {"data_diel_baby.root", "data_dimu_baby.root", "data_mueg_baby.root"};
    TString bkgroot[nBkg];
	double bkgSF[nBkg];

	bkgroot[0] = "DY1to4Jeemm_baby.root";    	bkgSF[0] = scaleDYeemm;
	bkgroot[1] = "DY1to4Jtautau_baby.root";	 	bkgSF[1] = scaleDYtautau;
	bkgroot[2] = "diboson_baby.root";		 	bkgSF[2] = scaleVV;
	bkgroot[3] = "tW_lepdl_baby.root";		 	bkgSF[3] = scaletw;
	bkgroot[4] = "tW_lepfake_baby.root";	 	bkgSF[4] = scalefake;
	bkgroot[5] = "tW_lepsl_baby.root";		 	bkgSF[5] = scalefake;
	bkgroot[6] = "triboson_baby.root";		 	bkgSF[6] = scaleVV;
	bkgroot[7] = "ttV_baby.root";			 	bkgSF[7] = scaleVV;
	bkgroot[8] = "ttfake_mcatnlo_baby.root";	bkgSF[8] = scalefake;
	bkgroot[9] = "ttsl_mcatnlo_baby.root";	 	bkgSF[9] = scalefake;
	bkgroot[10] = "w1to4jets_baby.root";		bkgSF[10] = scalefake;

    Float_t observable, observable_gen, ttmass, ttRapidity2, tmass;
    Float_t observableMinus, observableMinus_gen;
    Double_t weight;
    Int_t Nsolns = 1;
	Int_t channel = -99;

    for (Int_t iVar = 0; iVar < nVars; iVar++)
	  {

        Initialize1DBinning(iVar);

		// Figure out all possible binning schemes based on the basic binning scheme
		int nbinsx_gen = -99;
		int nbinsx_reco = -99;
		int nbinsx_reco_split_ch = -99;

		if( iVar < 2 || iVar == 9 ) nbinsx_gen = nbinsx2Dalt;
		else nbinsx_gen = nbins1D;

		nbinsx_reco = nbinsx_gen*2;
		//nbinsx_reco = nbinsx_gen;
		nbinsx_reco_split_ch = nbinsx_reco*nSig;

		double* genbins;
		double* recobins;
		double* recobins_split_ch;

		genbins = new double[nbinsx_gen+1];
		recobins = new double[nbinsx_reco+1];
		recobins_split_ch = new double[nbinsx_reco_split_ch+1];

		// Make gen binning array
        if( iVar < 2 || iVar == 9 ) {
          for( int i=0; i<=nbinsx2Dalt; i++ ) {
            genbins[i] = xbins2Dalt[i];
          }
        }
        else {
          for( int i=0; i<=nbins1D; i++ ) {
            genbins[i] = xbins1D[i];
          }
        }

		// Make reco binning array
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
		std::copy( recobins, recobins+nbinsx_reco+1, recobins_split_ch );
		for( int i=nbinsx_reco+1; i<nbinsx_reco_split_ch+1; i++ ) {
		  recobins_split_ch[i] = recobins_split_ch[i-nbinsx_reco] + recohist_width;
		}


        bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;

		// Make our histograms!
        TH1D *hData = new TH1D ("Data", "Data with background subtracted",    nbinsx_reco, recobins);
        TH1D *hData_split = new TH1D ("Data_split", "Data split",    nbinsx_reco_split_ch, recobins_split_ch);
        //TH1D *hData_combined = new TH1D ("Data_combined", "Data with all channels combined", nbinsx_reco, recobins);
        TH1D *hBkg = new TH1D ("Background",  "Background",    nbinsx_reco, recobins);
        TH1D *hBkg_split = new TH1D ("Background_split",  "Background_split",    nbinsx_reco_split_ch, recobins_split_ch);
        //TH1D *hBkg = new TH1D ("Background",  "Background",    nbinsx_reco_split_ch, recobins_split_ch);
        TH1D *hData_unfolded = new TH1D ("Data_Unfold", "Data with background subtracted and unfolded", nbinsx_gen, genbins);

        TH1D *hTrue_split = new TH1D ("true_split", "Truth",    nbinsx_reco_split_ch, recobins_split_ch); //Rebinned below
		TH1D *hTrue = new TH1D ("true_combined", "Truth",    nbinsx_gen, genbins);
        TH1D *hMeas = new TH1D ("meas", "Measured", nbinsx_reco, recobins);
        //TH1D *hMeas = new TH1D ("meas", "Measured", nbinsx_reco_split_ch, recobins_split_ch);
        TH1D *denominatorM_nopTreweighting = new TH1D ("denominatorM_nopTreweighting", "denominatorM_nopTreweighting", nbinsx_gen, genbins);
		TH1D *hPurity = new TH1D("purity", "Purity", nbinsx_gen, genbins);
		TH1D *hStability = new TH1D("stability", "Stability", nbinsx_gen, genbins);

        TH2D *hTrue_vs_Meas = new TH2D("true_vs_meas", "True vs Measured", nbinsx_reco, recobins, nbinsx_gen, genbins);
        //TH2D *hTrue_vs_Meas = new TH2D("true_vs_meas", "True vs Measured", nbinsx_reco_split_ch, recobins_split_ch, nbinsx_gen, genbins);

        TH1D *hData_bkgSub;
        TH1D *hData_bkgSub_split;

		hTrue_split->Rebin(2); //This gives us the gen bin-widths, while preserving the split channel layout.

		int ntau = 1000;
		double taulogrange = 2.;
		
		TH1D *AFB_vs_tau = new TH1D("AFB_vs_tau", "AFB_vs_tau", ntau, -taulogrange, taulogrange);
		TH1D *AFB_vs_tau_u = new TH1D("AFB_vs_tau_u", "AFB_vs_tau_u", ntau, -taulogrange, taulogrange);
		TH1D *AFB_vs_tau_d = new TH1D("AFB_vs_tau_d", "AFB_vs_tau_d", ntau, -taulogrange, taulogrange);
		TH1D *AFB_vs_tau0_c = new TH1D("AFB_vs_tau0_c", "AFB_vs_tau0_c", ntau, -taulogrange, taulogrange);
		TH1D *AFB_vs_tau0_u = new TH1D("AFB_vs_tau0_u", "AFB_vs_tau0_u", ntau, -taulogrange, taulogrange);
		TH1D *AFB_vs_tau0_d = new TH1D("AFB_vs_tau0_d", "AFB_vs_tau0_d", ntau, -taulogrange, taulogrange);


		//delete[] genbins;
		//delete[] recobins;
		delete[] recobins_split_ch;

        TMatrixD m_smearingE (nbinsx_gen, nbinsx_gen);
        TMatrixD m_correctE(nbinsx_gen, nbinsx_gen);
        TMatrixD m_unfoldcorr (nbinsx_gen, nbinsx_gen);

        //  Now test with data and with BKG subtraction

        TChain *ch_bkg[nBkg];
        TChain *ch_top = new TChain("tree");
        TChain *ch_data = new TChain("tree");


		//ch_data->Add(path + "ttdl_mcatnlo_baby.root");
        //ch_data->Add(path + "data.root");
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


        ch_data->SetBranchAddress(observablename,    &observable);
        if ( combineLepMinus ) ch_data->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        ch_data->SetBranchAddress("weight", &weight);
        //ch_data->SetBranchAddress("Nsolns", &Nsolns);
        ch_data->SetBranchAddress("tt_mass", &ttmass);
        ch_data->SetBranchAddress("ttRapidity2", &ttRapidity2);
        ch_data->SetBranchAddress("t_mass", &tmass);
		ch_data->SetBranchAddress("channel", &channel);

		double offset = 0;
		double histmax = recobins[nbinsx_reco];
		double histmin = recobins[0];
		double hiBinCenter = hData->GetBinCenter(nbinsx_reco);
		double loBinCenter = hData->GetBinCenter(1);

        for (Int_t i = 0; i < ch_data->GetEntries(); i++)
		  {
            ch_data->GetEntry(i);

			// If we're unfolding a single channel, skip events that aren't in that channel
			if( iChan < nSig && channel != iChan) continue;

			// Calculate a correction to the asymmetry value, to put it in the correct bin in a 2x1 histogram.
			//Use an offset to sort events into 2 superbins: ee/mumu, emu
			//if(channel>0) channel--;
			offset = double(channel) * recohist_width;
			//Do the same thing as "fillUnderOverflow", except adapted for 2x1 histograms
			if( observable > histmax )        observable = hiBinCenter;
			else if( observable < histmin )   observable = loBinCenter;
			if( observableMinus > histmax )        observableMinus = hiBinCenter;
			else if( observableMinus < histmin )   observableMinus = loBinCenter;

            if ( iVar<2 || iVar==9 || tmass>0 )
			  {
                // leptonic asymmetries don't need valid top mass solution
                fillUnderOverFlow(hData, observable, weight, Nsolns);
                fillUnderOverFlow(hData_split, observable+offset, weight, Nsolns);
				if (combineLepMinus) {
				  fillUnderOverFlow(hData, observableMinus, weight, Nsolns);
				  fillUnderOverFlow(hData_split, observableMinus+offset, weight, Nsolns);
				}
			  }

		  }

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
		  {

            ch_bkg[iBkg]->SetBranchAddress(observablename,    &observable);
            if ( combineLepMinus ) ch_bkg[iBkg]->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
            ch_bkg[iBkg]->SetBranchAddress("weight", &weight);
            //ch_bkg[iBkg]->SetBranchAddress("Nsolns", &Nsolns);
            ch_bkg[iBkg]->SetBranchAddress("tt_mass", &ttmass);
            ch_bkg[iBkg]->SetBranchAddress("ttRapidity2", &ttRapidity2);
            ch_bkg[iBkg]->SetBranchAddress("t_mass", &tmass);
			ch_bkg[iBkg]->SetBranchAddress("channel", &channel);

            for (Int_t i = 0; i < ch_bkg[iBkg]->GetEntries(); i++)
			  {
                ch_bkg[iBkg]->GetEntry(i);
				if( iChan < nSig && channel != iChan) continue;

				offset = double(channel) * recohist_width;
                weight *= bkgSF[iBkg];

			    if( observable > histmax )        observable = hiBinCenter;
			    else if( observable < histmin )   observable = loBinCenter;
			    if( observableMinus > histmax )        observableMinus = hiBinCenter;
			    else if( observableMinus < histmin )   observableMinus = loBinCenter;

                if ( iVar<2 || iVar==9 || tmass > 0 )
				  {
					// leptonic asymmetries don't need valid top mass solution
                    fillUnderOverFlow(hBkg, observable, weight, Nsolns);
                    fillUnderOverFlow(hBkg_split, observable+offset, weight, Nsolns);
				    if (combineLepMinus) {
				      fillUnderOverFlow(hBkg, observableMinus, weight, Nsolns);
				      fillUnderOverFlow(hBkg_split, observableMinus+offset, weight, Nsolns);
				    }
				  }

			  }

		  }

        ch_top->SetBranchAddress(observablename,    &observable);
        ch_top->SetBranchAddress(observablename + "_gen", &observable_gen);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
        if ( combineLepMinus ) ch_top->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
        ch_top->SetBranchAddress("weight", &weight);
        //ch_top->SetBranchAddress("Nsolns", &Nsolns);
        ch_top->SetBranchAddress("tt_mass", &ttmass);
        ch_top->SetBranchAddress("ttRapidity2", &ttRapidity2);
        ch_top->SetBranchAddress("t_mass", &tmass);
		ch_top->SetBranchAddress("channel", &channel);

        for (Int_t i = 0; i < ch_top->GetEntries(); i++)
		  {
            ch_top->GetEntry(i);
			if( iChan < nSig && channel != iChan) continue;
			// Calculate a correction to the asymmetry value, to put it in the correct bin in our 3x1 histograms
			//if( channel > 0 ) channel --;
			offset = double(channel) * recohist_width;
			if( observable > histmax )        observable = hiBinCenter;
			else if( observable < histmin )   observable = loBinCenter;
			if( observableMinus > histmax )        observableMinus = hiBinCenter;
			else if( observableMinus < histmin )   observableMinus = loBinCenter;
			if( observable_gen > histmax )        observable_gen = hiBinCenter;
			else if( observable_gen < histmin )   observable_gen = loBinCenter;
			if( observableMinus_gen > histmax )        observableMinus_gen = hiBinCenter;
			else if( observableMinus_gen < histmin )   observableMinus_gen = loBinCenter;
            weight *= scalettdil;

            // if ( (acceptanceName == "lepChargeAsym") || (acceptanceName == "lepAzimAsym") || (acceptanceName == "lepAzimAsym2") )
            if ( iVar<2 || iVar==9 || tmass>0 )
			  {
                //response.Fill (observable, observable_gen, weight);
                fillUnderOverFlow(hMeas, observable, weight, Nsolns);
                fillUnderOverFlow(hTrue, observable_gen, weight, Nsolns);
                fillUnderOverFlow(hTrue_split, observable_gen+offset, weight, Nsolns);
                fillUnderOverFlow(hTrue_vs_Meas, observable, observable_gen, weight, Nsolns);
                if ( combineLepMinus )
				  {
                    //response.Fill (observableMinus, observableMinus_gen, weight);
                    fillUnderOverFlow(hMeas, observableMinus, weight, Nsolns);
                    fillUnderOverFlow(hTrue, observableMinus_gen, weight, Nsolns);
                    fillUnderOverFlow(hTrue_split, observableMinus_gen+offset, weight, Nsolns);
                    fillUnderOverFlow(hTrue_vs_Meas, observableMinus, observableMinus_gen, weight, Nsolns);
				  }
			  }

			// if(i % 10000 == 0) cout<<i<<" "<<ch_top->GetEntries()<<endl;
		  }

		delete[] recobins;

        hData_bkgSub = (TH1D *) hData->Clone("Data_BkgSub");
        hData_bkgSub->Add(hBkg, -1.0);

        hData_bkgSub_split = (TH1D *) hData_split->Clone("Data_BkgSub_split");
        hData_bkgSub_split->Add(hBkg_split, -1.0);

		// Do the acceptance correction, by filling the migration matrix with events that have a gen-level value but no reco-level value
        TFile *file = new TFile("../denominator/acceptance/mcnlo/accept_" + acceptanceName + ".root");

		TH1D *acceptM[4];
		acceptM[0] = (TH1D*)(file->Get("accept_" + acceptanceName + "_diel"));
		acceptM[1] = (TH1D*)(file->Get("accept_" + acceptanceName + "_dimu"));
		acceptM[2] = (TH1D*)(file->Get("accept_" + acceptanceName + "_mueg"));
		acceptM[3] = (TH1D*)(file->Get("accept_" + acceptanceName + "_all" ));

		TH1D *accNum[4];
		accNum[0] = (TH1D*)(file->Get("numerator_" + acceptanceName + "_diel"));
		accNum[1] = (TH1D*)(file->Get("numerator_" + acceptanceName + "_dimu"));
		accNum[2] = (TH1D*)(file->Get("numerator_" + acceptanceName + "_mueg"));
		accNum[3] = (TH1D*)(file->Get("numerator_" + acceptanceName + "_all"));

		TH1D *accDen[4];
		accDen[0] = (TH1D*)(file->Get("denominator_" + acceptanceName + "_diel"));
		accDen[1] = (TH1D*)(file->Get("denominator_" + acceptanceName + "_dimu"));
		accDen[2] = (TH1D*)(file->Get("denominator_" + acceptanceName + "_mueg"));
		accDen[3] = (TH1D*)(file->Get("denominator_" + acceptanceName + "_all"));

		//Tricks to make "channel 0" hold the combined same-flavor histograms
		//accNum[0]->Add( accNum[1] );
		//accDen[0]->Add( accDen[1] );
		//acceptM[0] = (TH1D*)(accNum[0]->Clone("accept_SF"));
		//acceptM[0]->Divide( accDen[0] );


		if( iChan==nSig ) {     // Combined channels...

		  double gen_integrals[4] = {0.};
		  double reco_integrals[4] = {0.};

		  //adjust A matrix to match the channel proportions in data

		  //Figure out the relative proportion of events in each channel
		  for( int aChannel=0; aChannel<nSig; aChannel++ ) {
		    gen_integrals[aChannel] = hTrue_split->Integral( aChannel*nbinsx_gen+1, (aChannel+1)*nbinsx_gen );
		    reco_integrals[aChannel] = hData_bkgSub_split->Integral( aChannel*nbinsx_reco+1, (aChannel+1)*nbinsx_reco );
		    gen_integrals[nSig] += gen_integrals[aChannel];
		    reco_integrals[nSig] += reco_integrals[aChannel];
		  }


		  TH1D* acceptNumcorrected = (TH1D*)(accNum[nSig]->Clone("acceptNumcorrected"));
		  TH1D* acceptDencorrected = (TH1D*)(accDen[nSig]->Clone("acceptDencorrected"));

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

		  TH1D* acceptMcorrected = (TH1D*)(acceptNumcorrected->Clone("acceptMcorrected"));
		  acceptMcorrected->Divide( acceptDencorrected );

		  // to be fully correct, we should really reweight the events by correction[aChannel] when we fill hTrue and hTrue_vs_Meas too
		  for( int acceptbin=1; acceptbin<=nbinsx_gen; acceptbin++ ) {
		    double acceptance = acceptMcorrected->GetBinContent(acceptbin);
		    double n_accepted = hTrue->GetBinContent(acceptbin);
		    double n_rejected = n_accepted/acceptance - n_accepted;
		    hTrue_vs_Meas->SetBinContent( 0, acceptbin, n_rejected );
		    double num_error = accNum[iChan]->GetBinError(acceptbin);
		    double den_error = accDen[iChan]->GetBinError(acceptbin);
		    double new_error = sqrt( den_error*den_error - num_error*num_error );
		    hTrue_vs_Meas->SetBinError( 0, acceptbin, new_error );
		  }

		/* //I think this method only corrects the channel fraction in the acceptance denominator (not the numerator)
		  for( int aChannel=0; aChannel<nSig; aChannel++ ) {
			for( int acceptbin=1; acceptbin<=nbinsx_gen; acceptbin++ ) {

			  double old_error = hTrue_vs_Meas->GetBinError( 0, acceptbin );

			  //Calculate the number of rejected events, and fill it into the smearing matrix underflow
			  double acceptance = acceptM[aChannel]->GetBinContent(acceptbin);
			  double n_accepted = hTrue_split->GetBinContent(aChannel*nbinsx_gen + acceptbin);
			  double n_rejected = n_accepted/acceptance - n_accepted;

			  hTrue_vs_Meas->Fill( -999999, hTrue->GetXaxis()->GetBinCenter(acceptbin), n_rejected*correction[aChannel] );

			  //Calculate the uncertainty on the rejected events
			  double num_error = accNum[aChannel]->GetBinError(acceptbin);
			  double den_error = accDen[aChannel]->GetBinError(acceptbin);
			  double new_error = sqrt( old_error*old_error + den_error*den_error - num_error*num_error  );
			  hTrue_vs_Meas->SetBinError( 0, acceptbin, new_error );

			}
		  }
		  */
		}
		else {     // Individual channels...
		  for( int acceptbin=1; acceptbin<=nbinsx_gen; acceptbin++ ) {
			double acceptance = acceptM[iChan]->GetBinContent(acceptbin);
			double n_accepted = hTrue->GetBinContent(acceptbin);
			double n_rejected = n_accepted/acceptance - n_accepted;
			hTrue_vs_Meas->SetBinContent( 0, acceptbin, n_rejected );
			double num_error = accNum[iChan]->GetBinError(acceptbin);
			double den_error = accDen[iChan]->GetBinError(acceptbin);
			double new_error = sqrt( den_error*den_error - num_error*num_error );
			hTrue_vs_Meas->SetBinError( 0, acceptbin, new_error );
		  }
		}


		// Fill purity and stability plots
        if( hMeas->GetNbinsX() == nbinsx_gen ) {
		  for( int i=1; i<=nbinsx_gen; i++ ) {
			hPurity->SetBinContent( i, hTrue_vs_Meas->GetBinContent(i,i) / hMeas->GetBinContent(i) );
			hStability->SetBinContent( i, hTrue_vs_Meas->GetBinContent(i,i) / hTrue->GetBinContent(i) );
		  }
        }
        else if( hMeas->GetNbinsX() ==  2*nbinsx_gen || hMeas->GetNbinsX() == 4*nbinsx_gen ) {
		  for( int i=1; i<=nbinsx_gen; i++ ) {
			hPurity->SetBinContent( i, (hTrue_vs_Meas->GetBinContent(2*i,i)+hTrue_vs_Meas->GetBinContent(2*i-1,i)) / (hMeas->GetBinContent(2*i)+hMeas->GetBinContent(2*i-1)) );
			hStability->SetBinContent( i, (hTrue_vs_Meas->GetBinContent(2*i,i)+hTrue_vs_Meas->GetBinContent(2*i-1,i)) / hTrue->GetBinContent(i) );
		  }
        }
		else {
		  cout << "\n***WARNING: Purity and stability plots are broken!!!\n" << endl;
		}




        scaleBias =  (hData->Integral() - hBkg->Integral()) / hMeas->Integral() ;
        hMeas->Scale(scaleBias);

        TCanvas *c_reco = new TCanvas("c_reco", "c_reco", 500, 500);

        hData->SetLineWidth(lineWidth + 2);

        hMeas->SetLineColor(TColor::GetColorDark(kGreen));
        hMeas->SetFillColor(TColor::GetColorDark(kGreen));
        hMeas->SetFillStyle(3353);

        hBkg->SetLineColor(kYellow);
        hBkg->SetFillColor(kYellow);


        THStack *hMC = new THStack("hMC", "Stacked Top+BG");

        hMC->Add(hBkg);
        hMC->Add(hMeas);

        hMC->SetMinimum(0.0);
        hMC->SetMaximum( 1.5 * hMC->GetMaximum());
        hMC->Draw("hist");
        hMC->GetXaxis()->SetTitle(xaxislabel);
        hMC->GetYaxis()->SetTitleOffset(1.3);
        hMC->GetYaxis()->SetTitle("Events/bin");

        hData->Draw("E same");

        TLegend *leg0 = new TLegend(0.58, 0.75, 0.9, 0.93, NULL, "brNDC");
        leg0->SetEntrySeparation(100);
        leg0->SetFillColor(0);
        leg0->SetLineColor(0);
        leg0->SetBorderSize(0);
        leg0->SetTextSize(0.03);
        leg0->SetFillStyle(0);
        leg0->AddEntry(hData, "Data");
        leg0->AddEntry(hMeas,  "MC@NLO reco level", "F");
        leg0->AddEntry(hBkg,  "Background", "F");
        leg0->Draw();
        c_reco->SaveAs("1D_Reco_" + acceptanceName + "_" + channel_name + ".pdf");

		/////////////////////////////////////////////////////
		// Set data-like stat errors on MC for optimizing tau
		for( int i=1; i<=nbinsx_reco; i++) {
		  double n_sig = hMeas->GetBinContent(i);
		  double n_bkg = hBkg->GetBinContent(i);
		  double bkg_err = hBkg->GetBinError(i);
		  hMeas->SetBinError(i, sqrt(n_sig + n_bkg + bkg_err*bkg_err ) );
		}

		double tempScaleBias = hMeas->Integral() / hTrue->Integral();

		TUnfoldSys unfold_FindTau (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintArea);
		unfold_FindTau.SetInput(hMeas);
		minimizeRhoAverage(&unfold_FindTau, hMeas, -5.0, -0.5);
		tau = unfold_FindTau.GetTau();

		// Generate a curve of rhoAvg vs tau
		double ar_tau[90];
		double ar_rhoAvg[90];
		double tau_test = 0.0;
		double bestrhoavg = unfold_FindTau.GetRhoAvg();

		for(int l=0; l<90; l++) {
		  tau_test = pow( 10, -5.0 + 0.05*l );
		  unfold_FindTau.DoUnfold(tau_test, hMeas, tempScaleBias);
		  ar_tau[l] = tau_test;
		  ar_rhoAvg[l] = unfold_FindTau.GetRhoAvg();
		}

		TGraph* gr_rhoAvg = new TGraph(90, ar_tau, ar_rhoAvg);
		TCanvas* c_rhoAvg = new TCanvas("c_rhoAvg","c_rhoAvg");
		c_rhoAvg->SetLogx();
		gr_rhoAvg->SetTitle("Global Correlation Coefficient;#tau;#rho_{avg}");
		gr_rhoAvg->SetLineColor(kRed);
		gr_rhoAvg->Draw("al");

		TMarker* m_rhoMin = new TMarker(tau,bestrhoavg,kCircle);
		m_rhoMin->Draw();
		c_rhoAvg->SaveAs("1D_" + acceptanceName + "_" + channel_name + "_minRho.pdf");

		// cout << "Optimal tau value: " << tau << endl;
		// cout << "Minimum rho average: " << bestrhoavg << endl;

		/////// Do the unfolding! /////////////////////////////////////////////////////

		TUnfoldSys unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeCurvature, TUnfold::kEConstraintArea); //Set kEConstraintNone instead of kEConstraintArea here to give bug-free covariance matrix (for us, None gives the same result as Area anyway)
		unfold_TUnfold.SetInput(hData_bkgSub);
		//unfold_TUnfold.SetBias(hTrue);  //doesn't make any difference, because if not set the bias distribution is
		//automatically determined from hTrue_vs_Meas, which gives exactly hTrue


		// minimizeRhoAverage(&unfold_TUnfold, hData_bkgSub, -5.0, 0.0);
		// tau = unfold_TUnfold.GetTau();
		if(maketauplots) {

	        Float_t Afb0, AfbErr0;
			unfold_TUnfold.DoUnfold(tau, hData_bkgSub, scaleBias);

			unfold_TUnfold.GetOutput(hData_unfolded);

			TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
			TH2D *ematrix_smearing = (TH2D*)ematrix->Clone("ematrix_smearing");
			unfold_TUnfold.GetEmatrixSysUncorr( ematrix_smearing, 0, false );
			for (Int_t cmi = 0; cmi < nbinsx_gen; cmi++)
			  {
				for (Int_t cmj = 0; cmj < nbinsx_gen; cmj++)
				  {
					m_smearingE(cmi, cmj) = ematrix_smearing->GetBinContent(cmi + 1, cmj + 1);
					m_correctE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
				  }
			  }
		  	GetCorrectedAfb(hData_unfolded, m_correctE, Afb0, AfbErr0);
		  	//cout << " Unfolded0: " << Afb0 << " +/-  " << AfbErr0 << "\n";


			for (int itau = 0; itau < ntau; ++itau)
			{
				double taufact = 10.;
				double power = double(2*itau+1)/double(ntau) - 1.;
				power *= taulogrange;

				double taumult = pow(taufact,power);
				unfold_TUnfold.DoUnfold(tau*taumult, hData_bkgSub, scaleBias);

				unfold_TUnfold.GetOutput(hData_unfolded);

		        Float_t Afb, AfbErr;

				ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
				ematrix_smearing = (TH2D*)ematrix->Clone("ematrix_smearing");
				unfold_TUnfold.GetEmatrixSysUncorr( ematrix_smearing, 0, false );
				for (Int_t cmi = 0; cmi < nbinsx_gen; cmi++)
				  {
					for (Int_t cmj = 0; cmj < nbinsx_gen; cmj++)
					  {
						m_smearingE(cmi, cmj) = ematrix_smearing->GetBinContent(cmi + 1, cmj + 1);
						m_correctE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
					  }
				  }
			  	GetCorrectedAfb(hData_unfolded, m_correctE, Afb, AfbErr);

		        //GetAfb(hData_unfolded, Afb, AfbErr);
		        //cout <<taumult<< " Unfolded: " << Afb << " +/-  " << AfbErr << "\n";

		        AFB_vs_tau->Fill(power, Afb);
		        AFB_vs_tau_u->Fill(power, Afb+AfbErr);
		        AFB_vs_tau_d->Fill(power, Afb-AfbErr);
		        AFB_vs_tau0_c->Fill(power, Afb0);
		        AFB_vs_tau0_u->Fill(power, Afb0+AfbErr0);
		        AFB_vs_tau0_d->Fill(power, Afb0-AfbErr0);
		    }

		}


		//scaleBias = 0.0; //set biasScale to 0 when using kRegModeSize, or to compare with unfoldingType == 1
		//do the unfolding with calculated bias scale (N_data/N_MC), and tau from ScanLcurve if doScanLCurve=true. Note that the results will only be the same as unfoldingType == 1 with scaleBias=0 and the same value of tau.
		unfold_TUnfold.DoUnfold(tau, hData_bkgSub, scaleBias);

		unfold_TUnfold.GetOutput(hData_unfolded);

		TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
		TH2D *ematrix_smearing = (TH2D*)ematrix->Clone("ematrix_smearing");
		unfold_TUnfold.GetEmatrixSysUncorr( ematrix_smearing, 0, true );
		TH2D *cmatrix = unfold_TUnfold.GetRhoIJ("cmatrix", "correlation matrix", 0, 0);
		for (Int_t cmi = 0; cmi < nbinsx_gen; cmi++)
		  {
			for (Int_t cmj = 0; cmj < nbinsx_gen; cmj++)
			  {
				m_smearingE(cmi, cmj) = ematrix_smearing->GetBinContent(cmi + 1, cmj + 1);
				m_correctE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
				m_unfoldcorr(cmi, cmj) = cmatrix->GetBinContent(cmi + 1, cmj + 1);
			  }
		  }

        //m_unfoldE.Print("f=%1.5g ");
        //m_unfoldcorr.Print("f=%1.5g ");


		/*
        //include the MC stat uncertainty in the stat error bars
        for (int i = 1; i < nbinsx_gen + 1; i++)
		{
            hData_unfolded->SetBinError(i,sqrt(m_smearingE(i-1, i-1)));
		}
		*/


		float rmargin = gStyle->GetPadRightMargin();
		gStyle->SetPadRightMargin(0.13);

        TCanvas *c_resp = new TCanvas("c_resp", "c_resp");
		c_resp->SetLogz();
		TH2D *hResp = (TH2D*) hTrue_vs_Meas->Clone("response");
        gStyle->SetPalette(1);
        hResp->GetXaxis()->SetTitle(xaxislabel);
        hResp->GetYaxis()->SetTitle(xaxislabel + "_{gen}");
        hResp->Draw("COLZ");
        //c_resp->SaveAs("1D_Response_" + acceptanceName + ".eps");
        c_resp->SaveAs("1D_Response_" + acceptanceName + "_" + channel_name + ".pdf");
        //c_resp->SaveAs("Response_" + acceptanceName + ".C");
        //c_resp->SaveAs("Response_" + acceptanceName + ".root");
		gStyle->SetPadRightMargin(rmargin);

        TCanvas *c_purstab = new TCanvas("c_purstab", "c_purstab");
        hPurity->SetTitle("Purity;"+xaxislabel+";"+yaxislabel);
        hStability->SetTitle("Stability;"+xaxislabel+";"+yaxislabel);
        hPurity->Draw();
        c_purstab->SaveAs("1D_purity_" + acceptanceName + "_" + channel_name + ".pdf");
        hStability->Draw();
        c_purstab->SaveAs("1D_stability_" + acceptanceName + "_" + channel_name + ".pdf");


		TH1D *denominatorM = (TH1D*)file->Get("denominator_" + acceptanceName + "_" + ChannelName[iChan] ); //Acceptance denominator

		TFile *file_nopTreweighting = new TFile("../denominator/acceptance/mcnlo_nopTreweighting/accept_" + acceptanceName + ".root");
		TH1D *denominatorM_nopTreweighting_raw = (TH1D *) file_nopTreweighting->Get("denominator_" + acceptanceName + "_" + ChannelName[iChan] );

		//acceptM[iChan]->Scale(1.0 / acceptM[iChan]->Integral());

		for (Int_t i = 1; i <= nbinsx_gen; i++)
		  {

			if (acceptM[iChan]->GetBinContent(i) != 0)
			  {
				hTrue->SetBinContent(i, hTrue->GetBinContent(i) * 1.0 / acceptM[iChan]->GetBinContent(i));
				hTrue->SetBinError (i, hTrue->GetBinError(i) * 1.0 / acceptM[iChan]->GetBinContent(i));
			  }

			denominatorM_nopTreweighting->SetBinContent(i, denominatorM_nopTreweighting_raw->GetBinContent(i));

		  }

		denominatorM_nopTreweighting->Scale(1. / denominatorM_nopTreweighting->Integral(), "width");



        TProfile *theoryProfileCorr = new TProfile("thprofilecorrelated", "correlated data from theory file", nbinsx_gen, genbins);
        TProfile *theoryProfileUnCorr = new TProfile("thprofileuncorrelated", "uncorrelated data from theory file", nbinsx_gen, genbins);
        delete[] genbins;

		TF1 *predict_corr = new TF1("prediction", "pol1", -1, 1);
		TF1 *predict_uncorr = new TF1("flat_line", "pol0", -1, 1);
		predict_corr->SetParameter( 0, 0.5 );
		predict_uncorr->SetParameter( 0, 0.5 );

		if( acceptanceName == "lepCosTheta" ) predict_corr->SetParameter( 1, 0.0015 );
		else if( observablename == "lep_cos_opening_angle" ) predict_corr->SetParameter( 1, 0.1085 );

        if (observablename == "lep_azimuthal_asymmetry2")
		  {

            Float_t dphi, v1, v2, v3;
            Int_t ncols, nlines;
            nlines = 0;
            FILE *fp = fopen("theory/lhc8-dphill-corr.dat", "r");
            while (1)
			  {
                ncols = fscanf(fp, "%f %f %f %f", &dphi, &v1, &v2, &v3);
                if (ncols < 0) break;
                if (nlines < 5) printf("dphi=%8f, v=%8f\n", dphi, v1);
                theoryProfileCorr->Fill(dphi, v1);
                nlines++;
			  }

            nlines = 0;
            fp = fopen("theory/lhc8-dphill-uncorr.dat", "r");
            while (1)
			  {
                ncols = fscanf(fp, "%f %f %f %f", &dphi, &v1, &v2, &v3);
                if (ncols < 0) break;
                if (nlines < 5) printf("dphi=%8f, v=%8f\n", dphi, v1);
                theoryProfileUnCorr->Fill(dphi, v1);
                nlines++;
			  }
		  }

        if (observablename == "top_spin_correlation")
		  {

            Float_t c1c2, v1, v2, v3;
            Int_t ncols, nlines;
            nlines = 0;
            FILE *fp = fopen("theory/Chelbin-8TeV.dat", "r");
            while (1)
			  {
                ncols = fscanf(fp, "%f %f %f %f", &c1c2, &v1, &v2, &v3);
                if (ncols < 0) break;
                if (nlines < 5) printf("c1c2=%8f, v=%8f\n", c1c2, v1);
                theoryProfileCorr->Fill(c1c2, v1);
                nlines++;
			  }

            nlines = 0;
            fp = fopen("theory/lhc7_uncorr_mu1m_cos1cos2.dat", "r");
            while (1)
			  {
                ncols = fscanf(fp, "%f %f", &c1c2, &v1);
                if (ncols < 0) break;
                if (nlines < 5) printf("c1c2=%8f, v=%8f\n", c1c2, v1);
                theoryProfileUnCorr->Fill(c1c2, v1);
                nlines++;
			  }

		  }



        //==================================================================
        //============== Print the asymetry ================================
        cout << "========= Variable: " << acceptanceName << "===================\n";

		cout << "Automated tau value: " << tau << endl;
		cout << "Minimum rhoAverage: " << bestrhoavg << endl;
		cout << "bias scale for TUnfold: " << scaleBias << endl;

        Float_t Afb, AfbErr;

        GetAfb(hData, Afb, AfbErr);
        cout << " Data: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hTrue, Afb, AfbErr);
        cout << " True Top: " << Afb << " +/-  " << AfbErr << "\n";

        GetCorrectedAfb(hData_unfolded, m_correctE, Afb, AfbErr);
        cout << " Unfolded: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " Unfolded: " << Afb << " +/-  " << AfbErr << endl;

        GetCorrectedAfb(hData_unfolded, m_smearingE, Afb, AfbErr);
        cout << " Unfolded with only smearing errors: " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(denominatorM, Afb, AfbErr);
        cout << " True Top from acceptance denominator: " << Afb << " +/-  " << AfbErr << "\n";
        second_output_file << acceptanceName << " " << observablename << " True_Top_from_acceptance_denominator: " << Afb << " +/-  " << AfbErr << "\n";

        //GetAfbBinByBin(hData_unfolded);

        //GetAfb(hData_unfolded, Afb, AfbErr);
        //cout<<" Unfolded (ignoring correlation): "<< Afb <<" +/-  "<< AfbErr<<"\n";


        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation")
		  {
            GetAfb_integratewidth( (TH1D *) theoryProfileCorr, Afb, AfbErr);
            cout << " Bernreuther correlated: " << Afb << " +/-  " << AfbErr << "\n";
            GetAfb_integratewidth( (TH1D *) theoryProfileUnCorr, Afb, AfbErr);
            cout << " Bernreuther uncorrelated: " << Afb << " +/-  " << AfbErr << "\n";
		  }

        vector<double> afb_bins;
        vector<double> afb_bins_err;
        GetCorrectedAfbBinByBin(hData_unfolded, m_correctE, afb_bins, afb_bins_err, second_output_file);

        //scale to total xsec with option "width",  so that differential xsec is plotted
        //hData_unfolded->Scale(xsection/hData_unfolded->Integral(),"width");
        //hTrue->Scale(xsection/hTrue->Integral(),"width");


        TH1D* hData_unfolded_clone = (TH1D *) hData_unfolded->Clone("Data_unfolded_clone");

        hData_unfolded->Scale(1. / hData_unfolded->Integral(), "width");
        hTrue->Scale(1. / hTrue->Integral(), "width");

        //calculate covariance matrix for normalised distribution
        for (int l = 0; l < nbinsx_gen; l++)
		  {
            for (int j = 0; j < nbinsx_gen; j++)
			  {
                m_correctE(l, j) = m_correctE(l, j) * (hData_unfolded->GetBinContent(l + 1) / hData_unfolded_clone->GetBinContent(l + 1)) * (hData_unfolded->GetBinContent(j + 1) / hData_unfolded_clone->GetBinContent(j + 1)); //this gives the covariance matrix for the bin values
                //m_correctE(l, j) = m_correctE(l, j) * (hData_unfolded->GetBinWidth(l + 1) * hData_unfolded->GetBinContent(l + 1) / hData_unfolded_clone->GetBinContent(l + 1)) * (hData_unfolded->GetBinWidth(j + 1) * hData_unfolded->GetBinContent(j + 1) / hData_unfolded_clone->GetBinContent(j + 1)); //this gives the covariance matrix for the integrated bin contents
                m_smearingE(l, j) = m_smearingE(l, j) * (hData_unfolded->GetBinContent(l + 1) / hData_unfolded_clone->GetBinContent(l + 1)) * (hData_unfolded->GetBinContent(j + 1) / hData_unfolded_clone->GetBinContent(j + 1)); //this gives the covariance matrix for the bin values
			  }
		  }


        for (int i = 1; i < nbinsx_gen + 1; i++)
		  {
            cout << i << " bin = " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
            second_output_file << acceptanceName << " " << observablename << " bin" << i << ": " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
            //second_output_file << acceptanceName << " " << observablename << " truthbin" << i << ": " << hTrue->GetBinContent(i) << " +/- " << hTrue->GetBinError(i) << endl;
		  }


        for (int i = 1; i < nbinsx_gen + 1; i++)
		  {
            //cout << i << " bin = " << hData_unfolded->GetBinContent(i) << " +/- " << hData_unfolded->GetBinError(i) << endl;
            second_output_file << acceptanceName << " " << observablename << " bwidth" << i << ": " << hData_unfolded->GetBinWidth(i) << " +/- " << 0 << endl;
		  }


		cout << "Statistical covariance matrix:" << endl;
        // m_correctE.Print("f=%1.5g ");
		char mystring[15];

		cout << "_____|";
		for( int col=0; col<nbinsx_gen; col++ ) cout << "_____" << col+1 << "_____|";
		cout << endl;

		for( int row=0; row<nbinsx_gen; row++ ) {
		  sprintf(mystring, "%4d | ", row+1);
		  cout << mystring;
		  for( int col=0; col<nbinsx_gen; col++ ){
			sprintf(mystring, "%1.6g  ", m_correctE(row, col) );
			cout << mystring;
		  }
		  cout << endl;
		}

		cout << "Statistical covariance matrix for MC stats:" << endl;

		cout << "_____|";
		for( int col=0; col<nbinsx_gen; col++ ) cout << "_____" << col+1 << "_____|";
		cout << endl;

		for( int row=0; row<nbinsx_gen; row++ ) {
		  sprintf(mystring, "%4d | ", row+1);
		  cout << mystring;
		  for( int col=0; col<nbinsx_gen; col++ ){
			sprintf(mystring, "%1.6g  ", m_smearingE(row, col) );
			cout << mystring;
		  }
		  cout << endl;
		}



        //confirm covariance matrix for normalised distribution is correct by re-calculating Afb
        GetCorrectedAfb_integratewidth_V(hData_unfolded, m_correctE, Afb, AfbErr); //uses covariance matrix for the bin values
        //GetCorrectedAfb_integratewidth(hData_unfolded, m_correctE, Afb, AfbErr); //uses covariance matrix for the integrated bin contents
        cout << " Unfolded_after_scaling: " << Afb << " +/-  " << AfbErr << "\n";

        TH1D *hData_unfolded_minussyst;
        TH1D *hData_unfolded_plussyst;
        hData_unfolded_minussyst = (TH1D *) hData_unfolded->Clone("Data_unfolded_minussyst");
        hData_unfolded_plussyst = (TH1D *) hData_unfolded->Clone("Data_unfolded_plussyst");

        for (Int_t i = 1; i <= nbinsx_gen; i++)
		  {
		  	syst_corr[i - 1] = sqrt(syst_corr[i - 1]*syst_corr[i - 1] + m_smearingE(i-1, i-1)); //include MC stat uncertainty in systematics

            hData_unfolded_minussyst->SetBinContent(i, hData_unfolded->GetBinContent(i) - syst_corr[i - 1] );  //hard-coded systs
            hData_unfolded_minussyst->SetBinError(i, 0);
            hData_unfolded_plussyst ->SetBinContent(i, 2 * syst_corr[i - 1]);  //hard-coded systs
            hData_unfolded_plussyst ->SetBinError(i, 0);
		  }

        THStack *hs = new THStack("hs_systband", "Systematic band");
        hData_unfolded_minussyst->SetLineColor(10);
        hData_unfolded_minussyst->SetFillColor(10);
        hData_unfolded_minussyst->SetFillStyle(0);
        hs->Add(hData_unfolded_minussyst);
        hData_unfolded_plussyst->SetFillStyle(3353);
        hData_unfolded_plussyst->SetLineColor(kWhite);
        hData_unfolded_plussyst->SetFillColor(15);
        hs->Add(hData_unfolded_plussyst);
        //hs->SetMinimum( 0 );
        Float_t hsmin = hData_unfolded->GetMinimum() - ( 0.2 * hData_unfolded->GetMaximum() ) > 0.10 ? hData_unfolded->GetMinimum() - ( 0.2 * hData_unfolded->GetMaximum() ) : 0;
        hs->SetMinimum( hsmin );
        //printf("min = %8f \n", hsmin);
        if( int(200.*hsmin) - 10*int(20.*hsmin) > 8 ) hs->SetMinimum( hsmin + 0.01 ); //avoid axis label very close to bottom

        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation") hs->SetMaximum(1.35 * hData_unfolded->GetMaximum());
        else hs->SetMaximum(1.3 * hData_unfolded->GetMaximum());


        TCanvas *c_test;
        if(drawDiffs) c_test = new TCanvas("c_final", "c_final", 500, 650);
        else c_test = new TCanvas("c_final", "c_final", 500, 500);

        float r = 0.31;
        float epsilon = 0.02;

        TPad *p1, *p2;
        TLine *line;
        if(!drawDiffs) {
            p1 = new TPad("p1", "dist", 0.0, 0.0, 1, 1.);
            p1->Draw();
            p1->cd();

        }
        if(drawDiffs) {
            p1 = new TPad("p1", "dist", 0.0, r - epsilon, 1, 1.);
            p2 = new TPad("p2", "diff", 0.0, 0.0, 1, r*(1. - epsilon));
            p1->Draw();
            p2->Draw();
            p1->cd();
        }

        theoryProfileCorr->SetLineColor(kBlue);
        theoryProfileCorr->SetLineWidth(lineWidthDiffs);
        theoryProfileCorr->SetMarkerStyle(1);

        theoryProfileUnCorr->SetLineColor(kBlue);
        theoryProfileUnCorr->SetLineWidth(lineWidthDiffs);
        theoryProfileUnCorr->SetLineStyle(2);
        theoryProfileUnCorr->SetMarkerStyle(1);

		predict_corr->SetLineColor(kBlue);
		predict_corr->SetLineWidth(lineWidthDiffs);
		predict_corr->SetLineStyle(1);

		predict_uncorr->SetLineColor(kBlue);
		predict_uncorr->SetLineWidth(lineWidthDiffs);
		predict_uncorr->SetLineStyle(2);

        hs->Draw();
        hs->GetXaxis()->SetTitle(xaxislabel);
        hs->GetYaxis()->SetTitle("1/#sigma d#sigma/d(" + xaxislabel + ")");
        hData_unfolded->GetXaxis()->SetTitle(xaxislabel);
        //hData_unfolded->GetYaxis()->SetTitle("1/#sigma d#sigma/d("+xaxislabel+")");
        //hData_unfolded->SetMinimum(0.0);
        //hData_unfolded->SetMaximum( 2.0* hData_unfolded->GetMaximum());
        hData_unfolded->SetMarkerStyle(23);
        hData_unfolded->SetMarkerSize(1);
        hData_unfolded->SetFillStyle(0);
        hData_unfolded->Draw("E same");
        hData_unfolded->SetLineWidth(lineWidthDiffs);
        denominatorM_nopTreweighting->SetLineWidth(lineWidthDiffs);
        denominatorM_nopTreweighting->SetLineColor(TColor::GetColorDark(kRed));
        denominatorM_nopTreweighting->SetFillStyle(0);
        hTrue->SetLineWidth(lineWidthDiffs);
        hTrue->SetLineColor(TColor::GetColorDark(kRed));
        //hTrue->SetFillColor(TColor::GetColorDark(kGreen));
        hTrue->SetFillStyle(0);
        if (!draw_truth_before_pT_reweighting) hTrue->Draw("hist same");
        else denominatorM_nopTreweighting->Draw("hist same");
        hData_unfolded->Draw("EP same");
        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation")
		  {
            theoryProfileUnCorr->Draw("hist same");
            theoryProfileCorr->Draw("hist same");
		  }
		else if(acceptanceName == "lepCosTheta" || observablename == "lep_cos_opening_angle") {
		  predict_corr->Draw("LSAME");
		  if( observablename == "lep_cos_opening_angle" ) predict_uncorr->Draw("LSAME");
		}

        //TLegend* leg1=new TLegend(0.55,0.62,0.9,0.838,NULL,"brNDC");
        TLegend *leg1 = new TLegend(0.58, 0.75, 0.9, 0.93, NULL, "brNDC");
        leg1->SetEntrySeparation(0.1);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.032);
        leg1->AddEntry(hData_unfolded, "( Data - BG ) unfolded");
        leg1->AddEntry(hData_unfolded_plussyst,    "Syst. uncertainty", "F");
        leg1->AddEntry(hTrue,    "MC@NLO parton level", "L");

        leg1->Draw();

        if (observablename == "lep_azimuthal_asymmetry2" || observablename == "top_spin_correlation")
		  {
            TLegend *leg2 = new TLegend(0.18, 0.745, 0.45, 0.88, NULL, "brNDC");
            leg2->SetEntrySeparation(0.5);
            leg2->SetFillColor(0);
            leg2->SetLineColor(0);
            leg2->SetBorderSize(0);
            leg2->SetFillStyle(0);
            leg2->SetTextSize(0.032);
            leg2->AddEntry(theoryProfileCorr,  "#splitline{W.Bernreuther & Z.G.Si}{(8 TeV SM, #mu=^{}m_{t})}", "L");
            leg2->AddEntry(theoryProfileUnCorr,  "#splitline{W.Bernreuther & Z.G.Si}{(8 TeV uncorrelated, #mu=^{}m_{t})}", "L");
            leg2->Draw();
		  }
		else if(acceptanceName == "lepCosTheta" || observablename == "lep_cos_opening_angle") {
		  TLegend *leg2 = new TLegend(0.18, 0.745, 0.45, 0.88, NULL, "brNDC");
		  leg2->SetEntrySeparation(0.5);
		  leg2->SetFillColor(0);
		  leg2->SetLineColor(0);
		  leg2->SetBorderSize(0);
		  leg2->SetFillStyle(0);
		  leg2->SetTextSize(0.032);
		  leg2->AddEntry(predict_corr,  "#splitline{W.Bernreuther & Z.G.Si}{(8 TeV SM, #mu=^{}m_{t})}", "L");
		  if( observablename == "lep_cos_opening_angle" ) leg2->AddEntry(predict_uncorr,  "#splitline{W.Bernreuther & Z.G.Si}{(8 TeV uncorrelated, #mu=^{}m_{t})}", "L");
		  leg2->Draw();
		}





        TPaveText *pt1 = new TPaveText(0.175, 0.885, 0.41, 0.91, "brNDC");
        pt1->SetName("pt1name");
        pt1->SetBorderSize(0);
        pt1->SetFillStyle(0);

        TText *blah;
        //blah = pt1->AddText("CMS Preliminary, 19.5 fb^{-1} at  #sqrt{s}=8 TeV");
        blah = pt1->AddText("CMS, 19.5 fb^{-1} at  #sqrt{s}=8 TeV");
        blah->SetTextSize(0.032);
        blah->SetTextAlign(11);
        pt1->Draw();

        if(drawDiffs) {

            //p1->SetTopMargin(0.02);
            p1->SetBottomMargin(epsilon);
            p2->SetTopMargin(0.0);
            p2->SetBottomMargin(0.29);
            p2->SetFillColor(0);
            p2->SetFillStyle(0);


			p2->cd();
			int tempnbins = hData_unfolded->GetNbinsX();
			TString s_hname = "diff_";
			TH1D *h_diff = (TH1D *) hData_unfolded->Clone(s_hname +  acceptanceName);

			TH1D *h_diff_minussyst;
			TH1D *h_diff_plussyst;
			h_diff_minussyst = (TH1D *) hData_unfolded->Clone("Data_unfolded_minussyst");
			h_diff_plussyst = (TH1D *) hData_unfolded->Clone("Data_unfolded_plussyst");


			h_diff->Reset();
			line = new TLine(hData_unfolded->GetXaxis()->GetBinLowEdge(1), 1.0, 
			hData_unfolded->GetXaxis()->GetBinUpEdge(tempnbins), 1.0);
			h_diff->TH1D::Sumw2();  

            for(int tempbin = 1; tempbin < hData_unfolded->GetNbinsX()+1; tempbin++) {
				double mc, mcerr;

				if (!draw_truth_before_pT_reweighting) mc = hTrue->GetBinContent(tempbin);
				else mc = denominatorM_nopTreweighting->GetBinContent(tempbin);
				mcerr = 0.;
				
				double data = hData_unfolded->GetBinContent(tempbin);
				double dataerr = hData_unfolded->GetBinError(tempbin);
				double err2 = pow(mcerr*data/mc/mc,2) + pow(dataerr/mc,2);
				if(mc < 1e-10) 
				continue;
				h_diff->SetBinContent(tempbin, (data)/mc);
				h_diff->SetBinError(tempbin, sqrt(err2));         
				//h_diff->GetXaxis()->SetBinLabel(tempbin, hData_unfolded->GetXaxis()->GetBinLabel(tempbin));


				h_diff_minussyst->SetBinContent(tempbin, (data)/mc
				                                    - syst_corr[tempbin - 1]/mc );  //now includes MC stats
				h_diff_minussyst->SetBinError(tempbin, 0);
				h_diff_plussyst ->SetBinContent(tempbin, 2 * syst_corr[tempbin - 1]/mc );  //now includes MC stats
				h_diff_plussyst ->SetBinError(tempbin, 0);
            }

            THStack *hsd = new THStack("hsd_systband", "Systematic band");
            h_diff_minussyst->SetLineColor(kWhite);
            h_diff_minussyst->SetFillColor(kWhite);
            h_diff_minussyst->SetFillStyle(0);
            hsd->Add(h_diff_minussyst);
            h_diff_plussyst->SetFillStyle(3353);
            h_diff_plussyst->SetLineColor(kWhite);
            h_diff_plussyst->SetFillColor(15);
            hsd->Add(h_diff_plussyst);

            hsd->Draw();

            hsd->SetMinimum(0.921);
            hsd->SetMaximum(1.079/1.05); //because THStack multiplies the max when gStyle->SetHistTopMargin(0.) is not set

            if(h_diff->GetBinContent(h_diff->GetMaximumBin()) > 1.079) {
            	hsd->SetMaximum(h_diff->GetBinContent(h_diff->GetMaximumBin())*1.01/1.05);
            	hsd->SetMinimum(2. - h_diff->GetBinContent(h_diff->GetMaximumBin())*1.01);
            }

            if(h_diff->GetBinContent(h_diff->GetMinimumBin()) < 0.921  &&  (2. - h_diff->GetBinContent(h_diff->GetMinimumBin())*0.99) > (h_diff->GetBinContent(h_diff->GetMaximumBin())*1.01)  ) {
            	hsd->SetMaximum( (2. - h_diff->GetBinContent(h_diff->GetMinimumBin())*0.99)/1.05);
            	hsd->SetMinimum(h_diff->GetBinContent(h_diff->GetMinimumBin())*0.99);
            }

            //hsd->GetXaxis()->SetTitle( hs->GetXaxis()->GetTitle() );
            hsd->GetXaxis()->SetTitle(xaxislabel);
            hsd->GetXaxis()->SetTitleSize( hs->GetXaxis()->GetTitleSize()*(1.-r)/r);
            hsd->GetXaxis()->SetLabelSize( hs->GetXaxis()->GetLabelSize()*(1.-r)/r);
            //hsd->GetXaxis()->SetLabelOffset(-0.88);

            hsd->GetYaxis()->SetTitle("Data/Simulation ");
            hsd->GetYaxis()->SetNdivisions(805);
            //hsd->GetYaxis()->SetTitleFont(hData_unfolded->GetYaxis()->GetTitleFont());
            hsd->GetYaxis()->SetTitleOffset(0.7);
            hsd->GetYaxis()->SetTitleSize(0.120);
            hsd->GetYaxis()->SetLabelSize( hs->GetYaxis()->GetLabelSize()*(1.-r)/r);
            //hsd->GetYaxis()->SetLabelFont(hData_unfolded->GetYaxis()->GetLabelFont());

            //hsd->SetMarkerSize(0.8);

            //hsd->Draw("same");

        	hs->GetXaxis()->SetLabelSize(0.);
        	hs->GetXaxis()->SetTitleSize(0.);

            h_diff->Draw("Pesames");
            line->Draw();
            c_test->Modified();
            c_test->Update();


            p1->cd();   
        }


		//c_test->SaveAs("1D_finalplot_unfolded_" + acceptanceName + ".eps");
		c_test->SaveAs("1D_finalplot_unfolded_" + acceptanceName + "_" + channel_name + ".pdf");
        //c_test->SaveAs("finalplot_unfolded_" + acceptanceName + ".C");
        //c_test->SaveAs("finalplot_unfolded_" + acceptanceName + ".root");

        if(maketauplots){
	        TCanvas *c_AFBvsTau = new TCanvas("c_AFBvsTau", "c_AFBvsTau");
	        AFB_vs_tau->SetMaximum( AFB_vs_tau0_u->GetMaximum() + 1.*(AFB_vs_tau0_u->GetMaximum()-AFB_vs_tau0_d->GetMinimum()) );
	        AFB_vs_tau->SetMinimum( AFB_vs_tau0_d->GetMinimum() - 1.*(AFB_vs_tau0_u->GetMaximum()-AFB_vs_tau0_d->GetMinimum()) );
	        AFB_vs_tau->SetTitle("tau dependence;log_{10}(#tau/#tau_{0});"+asymlabel);
	        AFB_vs_tau0_c->SetLineColor(kBlue);
	        AFB_vs_tau0_u->SetLineColor(kRed);
	        AFB_vs_tau0_d->SetLineColor(kRed);
	        AFB_vs_tau_u->SetLineStyle(2);
	        AFB_vs_tau_d->SetLineStyle(2);
	        AFB_vs_tau->Draw("hist");
	        AFB_vs_tau0_c->Draw("hist same");
	        AFB_vs_tau0_d->Draw("hist same");
	        AFB_vs_tau0_u->Draw("hist same");
	        AFB_vs_tau_d->Draw("hist same");
	        AFB_vs_tau_u->Draw("hist same");
	        c_AFBvsTau->SaveAs("1D_tau_" + acceptanceName + "_" + channel_name + ".pdf");
	    }


        ch_data->Delete();

        ch_top->Delete();

        for (int iBkg = 0; iBkg < nBkg; ++iBkg)
		  {
            ch_bkg[iBkg]->Delete();
		  }
		file->Close();
	  }

    myfile.close();
    second_output_file.close();
  } // End of loop over channels

}

#ifndef __CINT__
int main ()
{
  AfbUnfoldExample();    // Main program when run stand-alone
  return 0;
}
#endif
