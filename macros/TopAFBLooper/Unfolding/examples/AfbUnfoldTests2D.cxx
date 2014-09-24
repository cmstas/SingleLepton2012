
#include <iostream>
#include <vector>

#include "TRandom3.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TUnfold.h"
#include "TUnfoldSys.h"
#include "TMatrixD.h"
#include "TChain.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TMarker.h"
#include "TFitResult.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPaveText.h"

#include "AfbFinalUnfold.h"
#include "tdrstyle.C"

using std::cout;
using std::endl;


//==============================================================================
// Global definitions
//==============================================================================


TString Region = "";


Int_t kterm = 3;
Double_t tau = 0.003;
Int_t nPseudos = 1;   // Linearity tests can use 1. For pull width tests, normally set to 10k
Int_t includeSys = 0;

// int lineWidth = 5;
bool plot_inclusive_only = false;


TF1 *fx;

Double_t myfunction(Double_t * x, Double_t * par)
{
  Float_t xx = x[0];
  Float_t ac = (xmax + xmin) / 2.;
  Double_t fscaled = par[1] * ( 1. + sign(xx - ac) * par[0] * ( fx->Eval( fabs( (xx - ac) / (xmax - ac) ) ) ) ) ;
  return fscaled;
}

//TestType: "Pull" or "Linearity"
//slopeOption: 0 = continuous reweighting, 1 = 6-binned reweighting
void AfbUnfoldTests(Int_t iVar = 0, TString TestType = "Linearity", /*Int_t slopeOption = 0,*/ Int_t Nfunction = 0, TString Var2D = "mtt")
{
    TH1::SetDefaultSumw2();


    //    TF1 *fsin = new TF1("fsin","sin(TMath::Pi()*x)",-1,1);
    //    TF1 *fcos = new TF1("fcos","cos(TMath::Pi()*x)",-1,1);
    //    TF1 *fmx2 = new TF1("fmx2","1 - (2*x - 1)^2",-1,1);
    //    TF1 *fx = new TF1("fx","x",-1,1);
    //    TF1 *fx2 = new TF1("fx2","x^2",-1,1);
    //    TF1 *fx3 = new TF1("fx3","x^3",-1,1);
    //    TF1 *fxhalf = new TF1("fxhalf","x^0.5",-1,1);
    //    TF1 *fxquarter = new TF1("fxquarter","x^0.25",-1,1);
    //    TF1 *fexpx = new TF1("fexpx","exp(x)/exp(1)",-1,1);


    if (Nfunction == 0) fx = new TF1("fx", "x", -1, 1);
    if (Nfunction == 1) fx = new TF1("fx", "x^2", -1, 1);
    if (Nfunction == 2) fx = new TF1("fx", "x^3", -1, 1);
    if (Nfunction == 3) fx = new TF1("fx", "x^0.5", -1, 1);
    if (Nfunction == 4) fx = new TF1("fx", "x^0.25", -1, 1);
    if (Nfunction == 5) fx = new TF1("fx", "exp(x)/exp(1)", -1, 1);
    if (Nfunction == 6) fx = new TF1("fx", "1 - (2*x - 1)^2", -1, 1);
    if (Nfunction == 7) fx = new TF1("fx", "sin(TMath::Pi()*x)", -1, 1);
    if (Nfunction == 8) fx = new TF1("fx", "cos(TMath::Pi()*x)", -1, 1);
    if (Nfunction == 9) fx = new TF1("fx", "1/2", -1, 1);

    setTDRStyle();
    gStyle->SetOptTitle(0);
    gStyle->SetOptFit();
    gStyle->SetOptStat("emr");
    cout.precision(3);

    //Initialize1DBinning(iVar);
	if (Var2D == "mtt") Initialize2DBinning(iVar);
	else if (Var2D == "ttrapidity2") Initialize2DBinningttrapidity2(iVar);
	else if (Var2D == "ttpt") Initialize2DBinningttpt(iVar);

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


    bool combineLepMinus = acceptanceName == "lepCosTheta" ? true : false;
    // const int nDiff = nbins1D + nbins1D / 2;
    // const int nDiff = nbinsy2D;

    ofstream myfile;
    myfile.open (acceptanceName + "_summary_PEtest_2D.txt");
    cout.rdbuf(myfile.rdbuf());
    ofstream second_output_file;
    second_output_file.open(acceptanceName + "_summary_PEtest_2D_formated.txt");

    // TString FunctionName = formatFloat(Nfunction, "%6.0i");
	char tmpchar[10];
	sprintf( tmpchar, "%6.0i", Nfunction );
	TString FunctionName(Nfunction);
    FunctionName.ReplaceAll(" " , "" );
    if (Nfunction == 0) FunctionName = "0";
    TString file_name = "p1_f" + FunctionName + "_" + acceptanceName + "_2D";
    ofstream third_output_file;
    third_output_file.open(file_name + ".txt");

    TRandom3 *random = new TRandom3();
    random->SetSeed(5);


    double asym_centre = (xmax + xmin) / 2.;


    // TF1 *fx_scaled = new TF1("fx_scaled", myfunction, xmin, xmax, 2);
    //fx_scaled->SetParameters(0.3,30);
    //fx_scaled->Eval(0.5);
    //cout<<"***** "<<fx_scaled->Eval(0.5)<<endl;
    // double fscale = -9999;

	///////////////////////////////////////////////////////////////////////////////////////////
	/////////////// 1. Set up all our histograms //////////////////////////////////////////////

	// 2D histograms
    TH2D *hTrue_before = new TH2D ("trueBeforeScaling", "Truth",    nbinsx_gen, genbins, nbinsy2D, ybins2D);
    TH2D *hTrue_before_split = new TH2D ("trueBeforeScaling_split", "Truth",    nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
    TH2D *hMeas_before = new TH2D ("measBeforeScaling", "Measured", nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);

    TH2D *hTrue_after = new TH2D ("trueAfterScaling", "Truth",    nbinsx_gen, genbins, nbinsy2D, ybins2D);
    TH2D *hMeas_after = new TH2D ("measAfterScaling", "Measured", nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
    TH2D *hMeas_after_combined = new TH2D ("measAfterScaling_combined", "Measured", nbinsx_reco, recobins, nbinsy2D, ybins2D);

    TH2D *hSmeared = new TH2D ("smeared", "Smeared", nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
    TH2D *hUnfolded = new TH2D ("unfolded", "Unfolded", nbinsx_gen, genbins, nbinsy2D, ybins2D);

	TH2D *hBkg = new TH2D ("Background",  "Background",    nbinsx_reco_3ch, recobins_3ch, nbinsy2D, ybins2D);
	TH2D *hBkg_combined = new TH2D ("Background_combined",  "Background combined", nbinsx_reco, recobins, nbinsy2D, ybins2D);

	// 1D histograms to store the unwrapped distributions
    TH1D *hTrue_before_unwrapped = new TH1D ("trueBeforeScalingUnwr", "Truth Before Unwrapped",    nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5); //Bias distribution
    TH1D *hMeas_before_unwrapped = new TH1D ("measBeforeScalingUnwr", "Measured Unwrapped", nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5); //For stat-error calculations
    TH1D *hSmeared_unwrapped = new TH1D ("smearedUnwr", "Smeared Unwrapped", nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5); //Input to TUnfold
	TH1D *hUnfolded_unwrapped = new TH1D ("unfoldedUnwr", "Unfolded Unwrapped", nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5); //Output from TUnfold
	TH1D *hBkg_unwrapped = new TH1D ("Background_Unwr",  "Background unwrapped",    nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5); //For doing background subtraction

    TH1D *hTrue_after_unwrapped = new TH1D ("trueAfterScalingUnwr", "Truth Unwrapped",    nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);
    // TH1D *hMeas_after_unwrapped = new TH1D ("measAfterScalingUnwr", "Measured Unwrapped", nbinsunwrapped, 0.5, double(nbinsunwrapped)+0.5);


    double pullMax = 5;
    int pullBins = 50;
    if (TestType == "Linearity") pullMax = 100;
    if (TestType == "Linearity") pullBins = 1000;

    TH1D* AfbPull[nbinsy2D + nbinsunwrapped_gen + 1];

    for (int iD = 0; iD < nbinsy2D + nbinsunwrapped_gen + 1; ++iD)
    {
	  char hname[15], htitle[18];
	  sprintf( hname, "h_afbpull_%d", iD );
	  sprintf( htitle, "Pulls for Afb %d", iD );
	  AfbPull[iD] = new TH1D(hname, htitle, pullBins, -pullMax, pullMax);
    }

	TH2D *hTrue_vs_Meas = new TH2D ("true_vs_meas", "True vs Measured", nbinsunwrapped_reco_3ch, 0.5, double(nbinsunwrapped_reco_3ch)+0.5, nbinsunwrapped_gen, 0.5, double(nbinsunwrapped_gen)+0.5);

	hTrue_before_split->RebinX(2);

    TMatrixD m_unfoldE(nbinsunwrapped_gen, nbinsunwrapped_gen);

	/*
    TH1F *h_pulls[nbins1D];
    TH1F *h_resd[nbins1D];
    for (int i = 0; i < nbins1D; i++)
    {
        TString name = "h_pull_";
        name += i;
        h_pulls[i] = new TH1F(name, name, 50, -5.0, 5.0);
        name = "h_resd_";
        name += i;
        h_resd[i] = new TH1F(name, name, 20, -1, 1);
    }
	*/

	/////////////////////////////////////////////////////////////////////////////////////////
	///////////////////////// 2. Fill our histograms from the baby ntuples //////////////////

    Float_t observable, observable_gen, tmass, ttmass, ttRapidity2;
    Float_t observableMinus, observableMinus_gen;
	Float_t obs2D, obs2D_gen;
    Double_t weight;
    Int_t Nsolns = 1;
	Int_t channel = -99;

	double offset = 0;
	double histmax = recobins[nbinsx_reco];
	double histmin = recobins[0];
	double hiBinCenter = hMeas_after_combined->GetXaxis()->GetBinCenter(nbinsx_reco);
	double loBinCenter = hMeas_after_combined->GetXaxis()->GetBinCenter(1);

	delete[] genbins;
	delete[] recobins;
	delete[] recobins_3ch;

	// Load in data, just to get the integrals for scaling purposes //////////////////////////
	TChain *ch_data = new TChain("tree");
	ch_data->Add("../data_diel_baby.root");
	ch_data->Add("../data_dimu_baby.root");
	ch_data->Add("../data_mueg_baby.root");
	double integral_data = -999.9;
	if( iVar<2 || iVar==9 ) integral_data = double(ch_data->GetEntries());
	else integral_data = double(ch_data->GetEntries("t_mass > 0 "));

	// Background events /////////////////////////////////////////////////////
	TChain *ch_bkg = new TChain("tree");
	ch_bkg->Add("../DY1to4Jtot_baby.root");
	ch_bkg->Add("../diboson_baby.root");
	ch_bkg->Add("../tW_lepdl_baby.root");
	ch_bkg->Add("../tW_lepfake_baby.root");
	ch_bkg->Add("../tW_lepsl_baby.root");
	ch_bkg->Add("../triboson_baby.root");
	ch_bkg->Add("../ttV_baby.root");
	ch_bkg->Add("../ttfake_powheg_baby.root");
	ch_bkg->Add("../ttsl_powheg_baby.root");
	ch_bkg->Add("../w1to4jets_baby.root");

	ch_bkg->SetBranchAddress(observablename,    &observable);
	if ( combineLepMinus ) ch_bkg->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
	ch_bkg->SetBranchAddress("weight", &weight);
	// ch_bkg->SetBranchAddress("Nsolns", &Nsolns);
	ch_bkg->SetBranchAddress("t_mass", &tmass);
	ch_bkg->SetBranchAddress("tt_mass", &ttmass);
	ch_bkg->SetBranchAddress("ttRapidity2", &ttRapidity2);
	ch_bkg->SetBranchAddress("channel", &channel);

	if (Var2D == "mtt")              ch_bkg->SetBranchAddress("tt_mass", &obs2D);
	else if (Var2D == "ttrapidity2") ch_bkg->SetBranchAddress("ttRapidity2", &obs2D);
	else if (Var2D == "ttpt")        ch_bkg->SetBranchAddress("ttPt", &obs2D);

	for (Int_t i = 0; i < ch_bkg->GetEntries(); i++) {
	  ch_bkg->GetEntry(i);
	  obs2D = fabs(obs2D);
	  //weight *= bkgSF;

	  //Adjustments for 3-channel unfolding
	  offset = double(channel) * recohist_width;
	  if( observable > histmax )             observable = hiBinCenter;
	  else if( observable < histmin )        observable = loBinCenter;
	  if( observableMinus > histmax )        observableMinus = hiBinCenter;
	  else if( observableMinus < histmin )   observableMinus = loBinCenter;

	  if ( tmass > 0 ) {
		fillUnderOverFlow(hBkg, observable+offset, obs2D, weight, Nsolns);
		fillUnderOverFlow(hBkg_combined, observable, obs2D, weight, Nsolns);
		if (combineLepMinus) {
		  fillUnderOverFlow(hBkg, observableMinus+offset, obs2D, weight, Nsolns);
		  fillUnderOverFlow(hBkg_combined, observableMinus, obs2D, weight, Nsolns);
		}
	  }
	}

	// top MC events ///////////////////////////////////////////////////////////////
    TFile *ttfile = new TFile("../ttdl_mcatnlo_smallTree_baby.root");
    TTree *evtree = (TTree *) ttfile->Get("tree");
    Int_t entries = (Int_t)evtree->GetEntries();
    cout << "RESPONSE: Number of Entries: " << entries << endl;

    evtree->SetBranchAddress(observablename,    &observable);
    evtree->SetBranchAddress(observablename + "_gen", &observable_gen);
    if ( combineLepMinus ) evtree->SetBranchAddress("lepMinus_costheta_cms",    &observableMinus);
    if ( combineLepMinus ) evtree->SetBranchAddress("lepMinus_costheta_cms_gen",    &observableMinus_gen);
    evtree->SetBranchAddress("weight", &weight);
    // evtree->SetBranchAddress("Nsolns", &Nsolns);
    evtree->SetBranchAddress("t_mass", &tmass);
	evtree->SetBranchAddress("tt_mass", &ttmass);
	evtree->SetBranchAddress("ttRapidity2", &ttRapidity2);
	evtree->SetBranchAddress("channel", &channel);

	if (Var2D == "mtt")
	  {
		evtree->SetBranchAddress("tt_mass", &obs2D);
		evtree->SetBranchAddress("tt_mass_gen", &obs2D_gen);
	  }
	else if (Var2D == "ttrapidity2")
	  {
		evtree->SetBranchAddress("ttRapidity2", &obs2D);
		evtree->SetBranchAddress("ttRapidity2_gen", &obs2D_gen);
	  }
	else if (Var2D == "ttpt")
	  {
		evtree->SetBranchAddress("ttPt", &obs2D);
		evtree->SetBranchAddress("ttPt_gen", &obs2D_gen);
	  }


	// Set up stuff for the linearity and/or pull tests /////////////////////////////////
    Float_t slope = 0.0;
    const int Nlin = 7;
    //Float_t A_gen[Nlin], Aerr_gen[Nlin], A_unf[Nlin], Aerr_unf[Nlin], A_meas[Nlin], Aerr_meas[Nlin];
    //Float_t A_pull[Nlin], A_pullwidth[Nlin], Aerr_pull[Nlin], Aerr_pullwidth[Nlin];

    // Float_t  A_meas[Nlin], Aerr_meas[Nlin];
    vector<vector<Float_t> > A_gen, Aerr_gen, A_unf, Aerr_unf;
    vector<vector<Float_t> > A_pull, A_pullwidth, Aerr_pull, Aerr_pullwidth;

    A_gen.clear();
    Aerr_gen.clear();
    A_unf.clear();
    Aerr_unf.clear();
    A_pull.clear();
    A_pullwidth.clear();
    Aerr_pull.clear();
    Aerr_pullwidth.clear();


    // TH2D *hTrue_after_array[Nlin];
    // TH2D *hMeas_after_array[Nlin];

	// Reweight our events, and fill the histograms //////////////////////////////////

	//Begin loop over k
    for (int k = 0; k < Nlin; k++)
    {

        if ((TestType == "Pull") && (k == 1)) break;

        slope = -0.3 + 0.1 * k;
        //fx_scaled->SetParameters(slope,1.);

        cout << "slope =" << slope << "\n";

        hTrue_before->Reset();
		hTrue_before_split->Reset();
        hMeas_before->Reset();
        hTrue_after->Reset();
        hMeas_after->Reset();
		hMeas_after_combined->Reset();
        hTrue_vs_Meas->Reset();
        for (int iD = 0; iD < nbinsy2D + nbinsunwrapped_gen + 1; ++iD)
        {
            AfbPull[iD]->Reset();
        }

		int measbin_3ch = -99;
		int measbin = -99;
		int genbin = -99;

        for (Int_t i = 0; i < entries; i++)
        {
            evtree->GetEntry(i);
            double orig_weight = weight;

			obs2D = fabs(obs2D);
			obs2D_gen = fabs(obs2D_gen);

            /*if (slopeOption == 1)
            {
			  //fix the observable values to the bin centres so the acceptance function is unaffected by any reweighting
			  observable =  hEmpty->GetXaxis()->GetBinCenter( hEmpty->FindBin( observable, obs2D ) );
			  observable_gen =  hEmpty_gen->GetXaxis()->GetBinCenter( hEmpty_gen->FindBin( observable_gen, obs2D_gen ) );
			  obs2D =  hEmpty->GetYaxis()->GetBinCenter( hEmpty->FindBin( observable, obs2D ) );
			  obs2D_gen =  hEmpty_gen->GetYaxis()->GetBinCenter( hEmpty_gen->FindBin( observable_gen, obs2D_gen ) );
			  if ( combineLepMinus )
                {
				  observableMinus =  hEmpty->GetXaxis()->GetBinCenter( hEmpty->FindBin( observableMinus, obs2D ) );
				  observableMinus_gen =  hEmpty_gen->GetXaxis()->GetBinCenter( hEmpty_gen->FindBin( observableMinus_gen, obs2D_gen ) );
                }
			}*/

			offset = double(channel) * recohist_width;
			if( observable > histmax )        observable = hiBinCenter;
			else if( observable < histmin )   observable = loBinCenter;
			if( observableMinus > histmax )        observableMinus = hiBinCenter;
			else if( observableMinus < histmin )   observableMinus = loBinCenter;
			if( observable_gen > histmax )        observable_gen = hiBinCenter;
			else if( observable_gen < histmin )   observable_gen = loBinCenter;
			if( observableMinus_gen > histmax )        observableMinus_gen = hiBinCenter;
			else if( observableMinus_gen < histmin )   observableMinus_gen = loBinCenter;

            double xval = (observable_gen - asym_centre) / fabs(xmax - asym_centre);
            double xsign = sign(xval);
            //restrict range from -1 to +1
            if ( fabs(xval) > 1. ) xval = xsign;

            double xminusval = -9999;
            double xminussign = -9999;

            if ( combineLepMinus )
            {
                xminusval = (observableMinus_gen - asym_centre) / fabs(xmax - asym_centre);
                xminussign = sign(xminusval);
                //restrict range from -1 to +1
                if ( fabs(xminusval) > 1. ) xminusval = xminussign;
            }

            //if(i % 10000 == 0) cout<<i<<" "<<evtree->GetEntries()<<endl;

            if ( tmass > 0 )
            {
			  genbin  = getUnwrappedBin(hTrue_before, observable_gen, obs2D_gen);
			  measbin = getUnwrappedBin(hMeas_after_combined, observable, obs2D);
			  measbin_3ch = measbin + (channel * nbinsunwrapped_reco);

			  fillUnderOverFlow(hMeas_before, observable+offset, obs2D, weight, Nsolns);
			  fillUnderOverFlow(hTrue_before, observable_gen, obs2D_gen, weight, Nsolns);
			  fillUnderOverFlow(hTrue_before_split, observable_gen+offset, obs2D_gen, weight, Nsolns);
			  fillUnderOverFlow(hTrue_vs_Meas, measbin_3ch, genbin, weight, Nsolns);
			  if ( combineLepMinus )
              {
				genbin  = getUnwrappedBin(hTrue_before, observableMinus_gen, obs2D_gen);
				measbin = getUnwrappedBin(hMeas_after_combined, observableMinus, obs2D);
				measbin_3ch = measbin + (channel * nbinsunwrapped_reco);

				fillUnderOverFlow(hMeas_before, observableMinus+offset, obs2D, weight, Nsolns);
				fillUnderOverFlow(hTrue_before, observableMinus_gen, obs2D_gen, weight, Nsolns);
				fillUnderOverFlow(hTrue_before_split, observableMinus_gen+offset, obs2D_gen, weight, Nsolns);
				fillUnderOverFlow(hTrue_vs_Meas, measbin_3ch, genbin, weight, Nsolns);
              }
			  //if (TestType == "Linearity") weight = weight * fx_scaled->Eval(observable_gen); //this is very slow for some reason
			  if (TestType == "Linearity") weight = weight * (1.0 + slope * xsign * ( fx->Eval(fabs(xval)) ) );
			  fillUnderOverFlow(hMeas_after, observable+offset, obs2D, weight, Nsolns);
			  fillUnderOverFlow(hMeas_after_combined, observable, obs2D, weight, Nsolns);
			  fillUnderOverFlow(hTrue_after, observable_gen, obs2D_gen, weight, Nsolns);
			  if ( combineLepMinus )
              {
				//if (TestType == "Linearity") weight = orig_weight * fx_scaled->Eval(observableMinus_gen); //this is very slow for some reason
				if (TestType == "Linearity") weight = orig_weight * (1.0 + slope * xminussign * ( fx->Eval(fabs(xminusval)) ) );
				fillUnderOverFlow(hMeas_after, observableMinus+offset, obs2D, weight, Nsolns);
				fillUnderOverFlow(hMeas_after_combined, observableMinus, obs2D, weight, Nsolns);
				fillUnderOverFlow(hTrue_after, observableMinus_gen, obs2D_gen, weight, Nsolns);
              }
            }

        } // end of loop over ttdil entries

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
			  double n_accepted = hTrue_before_split->GetBinContent( aChannel*nbinsx_gen + xBin, yBin );
			  genbin = (yBin-1)*nbinsx_gen + xBin;
			  hTrue_vs_Meas->Fill( -999999, double(genbin), n_accepted/acceptance - n_accepted );

			}
		  }
		}

		// Normalize top MC to data
		double integral_top = hTrue_before->Integral();
		double integral_signal = integral_data - hBkg->Integral();

		hTrue_before->Scale(		 integral_signal / integral_top );
		hTrue_before_split->Scale(	 integral_signal / integral_top );
        hMeas_before->Scale(		 integral_signal / integral_top );
        hTrue_after->Scale(			 integral_signal / integral_top );
        hMeas_after->Scale(			 integral_signal / integral_top );
		hMeas_after_combined->Scale( integral_signal / integral_top );

		// Optimize tau for use in all subsequent unfoldings
		if (k == 0) {

		  cout << "Optimizing tau..." << endl;

		  unwrap2dhisto(hTrue_before, hTrue_before_unwrapped);
		  unwrap2dhisto_3ch(hMeas_before, hMeas_before_unwrapped);
		  unwrap2dhisto_3ch(hBkg,         hBkg_unwrapped);

		  //Set data-like stat errors on MC for optimizing tau
		  for( int i=1; i<=nbinsunwrapped_reco_3ch; i++) {
			double n_sig = hMeas_before_unwrapped->GetBinContent(i);
			double n_bkg = hBkg_unwrapped->GetBinContent(i);
			double bkg_err = hBkg_unwrapped->GetBinError(i);
			hMeas_before_unwrapped->SetBinError(i, sqrt(n_sig + n_bkg + bkg_err*bkg_err ) );
			//hMeas_before_unwrapped->SetBinError(i, 2.0);
		  }

		  scaleBias = hMeas_before_unwrapped->Integral() / hTrue_before_unwrapped->Integral();

		  TUnfoldSys unfold_FindTau (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone, TUnfold::kEConstraintArea);
		  unfold_FindTau.SetInput(hMeas_before_unwrapped);
		  unfold_FindTau.RegularizeBins2D(1,1,nbinsx_gen,nbinsx_gen,nbinsy2D,TUnfold::kRegModeCurvature);
		  minimizeRhoAverage(&unfold_FindTau, hMeas_before_unwrapped, -6.0, -1.0);
		  tau = unfold_FindTau.GetTau();

		  // Generate a curve of rhoAvg vs log(tau)
		  double ar_tau[100];
		  double ar_rhoAvg[100];
		  double tau_test = 0.0;
		  double bestrhoavg = unfold_FindTau.GetRhoAvg();

		  for(int l=0; l<100; l++) {
			tau_test = pow( 10, -6.0 + 0.05*l );
			unfold_FindTau.DoUnfold(tau_test, hMeas_before_unwrapped, scaleBias);
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
		  c_rhoAvg->SaveAs("2D_" + acceptanceName + "_unfoldTests_minRho.pdf");

		  cout << "Optimal tau value: " << tau << endl;
		  cout << "Minimum rho average: " << bestrhoavg << endl;
		} //End of tau optimization

		////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////// 3. Begin testing ///////////////////////////////////////////////////////

        //scale to keep total yield constant
        hTrue_after->Scale( hTrue_before->Integral() / hTrue_after->Integral() );
        hMeas_after->Scale( hMeas_before->Integral() / hMeas_after->Integral() );
		hMeas_after_combined->Scale( hMeas_before->Integral() / hMeas_after_combined->Integral() );

        // fscale = 0.65 * double(hTrue_before->Integral()) / double(nbinsx2D*nbinsy2D);

        // hTrue_after_array[k] = (TH2D *) hTrue_after->Clone();
        // hMeas_after_array[k] = (TH2D *) hMeas_after_combined->Clone();

        for (Int_t x = 1; x <= nbinsx_gen; x++) {
		  for (Int_t y = 1; y<= nbinsy2D; y++) {
			if (acceptM[3]->GetBinContent(x,y) != 0) {
			  hTrue_after->SetBinContent(x,y, hTrue_after->GetBinContent(x,y) * 1.0 / acceptM[3]->GetBinContent(x,y));
			  hTrue_after->SetBinError  (x,y, hTrue_after->GetBinError(x,y)  * 1.0 / acceptM[3]->GetBinContent(x,y));
			}
		  }
		}

        Float_t Afb, AfbErr;

        GetAfb(hTrue_after, Afb, AfbErr);
        // Float_t A_gen_k = Afb;
        // Float_t Aerr_gen_k = 0.0;
        cout << " True after re-weighting   : " << Afb << " +/-  " << AfbErr << "\n";

        GetAfb(hMeas_after_combined, Afb, AfbErr);
        // A_meas[k] = Afb;
        // Aerr_meas[k] = AfbErr;
        cout << " Measured after re-weighting   : " << Afb << " +/-  " << AfbErr << "\n";


        // Now do the pseudos

        // Float_t trialAsym = 0.0;
        vector <Float_t> SumAsym, SumErrAsym, SumTrueAsym, SumTrueErrAsym;
        SumAsym.clear();
        SumErrAsym.clear();
        SumTrueAsym.clear();
        SumTrueErrAsym.clear();

        for (int iD = 0; iD < nbinsy2D + nbinsunwrapped_gen + 1; ++iD)
        {
            SumAsym.push_back(0);
            SumErrAsym.push_back(0);
            SumTrueAsym.push_back(0);
            SumTrueErrAsym.push_back(0);
        }


        if (nPseudos > 0)
        {

		  ////////// Begin loop over all pseudoexperiments ////////////
		  for (int i = 0; i < nPseudos; i++)
            {

			  for (int x = 1; x <= hMeas_after->GetNbinsX(); x++)
                {
				  for (int y = 1; y < hMeas_after->GetNbinsY() + 1; y++)
					{
					  double fluct;
					  if (nPseudos > 1) fluct = random->Poisson(hMeas_after->GetBinContent(x,y));
					  else fluct = hMeas_after->GetBinContent(x,y);
					  hSmeared->SetBinError(x, y, sqrt(fluct));
					  hSmeared->SetBinContent(x, y, fluct);
					}
                }

			  unwrap2dhisto_3ch(hSmeared, hSmeared_unwrapped);
			  unwrap2dhisto(hTrue_before, hTrue_before_unwrapped);

			  // Unfold! /////////////////////////////
			  TUnfold unfold_TUnfold (hTrue_vs_Meas, TUnfold::kHistMapOutputVert, TUnfold::kRegModeNone, TUnfold::kEConstraintArea);
			  scaleBias =  hSmeared->Integral() / hMeas_before->Integral() ;
			  cout << "bias scale for TUnfold: " << scaleBias << endl;
			  unfold_TUnfold.SetBias(hTrue_before_unwrapped);  //doesn't make any difference, because if not set the bias distribution is automatically determined from hTrue_vs_Meas, which gives exactly hTrue
			  unfold_TUnfold.SetInput(hSmeared_unwrapped);
			  unfold_TUnfold.RegularizeBins2D(1,1,nbinsx_gen,nbinsx_gen,nbinsy2D,TUnfold::kRegModeCurvature);
			  unfold_TUnfold.DoUnfold(tau, hSmeared_unwrapped, scaleBias);
			  unfold_TUnfold.GetOutput(hUnfolded_unwrapped);

			  TH2D *ematrix = unfold_TUnfold.GetEmatrix("ematrix", "error matrix", 0, 0);
			  for (Int_t cmi = 0; cmi < nbinsunwrapped_gen; cmi++)
				{
				  for (Int_t cmj = 0; cmj < nbinsunwrapped_gen; cmj++)
					{
					  m_unfoldE(cmi, cmj) = ematrix->GetBinContent(cmi + 1, cmj + 1);
					}
				}


			  rewrap1dhisto(hUnfolded_unwrapped, hUnfolded);

			  GetAfb(hUnfolded, Afb, AfbErr);

			  //-------------------------------------------//
			  //   DAN'S DEBUGGING OUTPUT                  //
			  //-------------------------------------------//
			  // cout << "*******Unfolded result: Afb=" << Afb << ", err=" << AfbErr << endl;
			  // if( k == 3 ) {
			  // 	gStyle->SetPadRightMargin(0.17);

			  // 	TCanvas *c_data = new TCanvas("c_data", "c_data", 675, 600);
			  // 	gStyle->SetPalette(1);
			  // 	hSmeared->SetTitle("Smeared");
			  // 	hSmeared->GetZaxis()->SetMoreLogLabels();
			  // 	hSmeared->Draw("COLZ");
			  // 	c_data->SaveAs("DEBUG_2D_smeared_" + acceptanceName +"_" + Var2D + ".pdf");
			  // 	hMeas_after->SetTitle("Meas_after");
			  // 	hMeas_after->Draw("COLZ");
			  // 	c_data->SaveAs("DEBUG_2D_meas_" + acceptanceName +"_" + Var2D + ".pdf");

			  // 	TCanvas *c_resp = new TCanvas("c_resp", "c_resp", 1850, 600);
			  // 	TH2D *hResp = (TH2D*) hTrue_vs_Meas->Clone("response");
			  // 	gStyle->SetPalette(1);
			  // 	hResp->SetTitle("Migration matrix");
			  // 	hResp->GetXaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (reco)");
			  // 	hResp->GetYaxis()->SetTitle(yaxislabel + " and " + xaxislabel + ", unwrapped (gen)");
			  // 	hResp->Draw("COLZ");
			  // 	c_resp->SetLogz();
			  // 	c_resp->SaveAs("DEBUG_2D_Response_" + acceptanceName + "_" + Var2D + ".pdf");

			  // 	TCanvas *c_output = new TCanvas("c_output", "c_output", 675, 600);
			  // 	gStyle->SetPalette(1);
			  // 	// c_output->SetLogz();
			  // 	hUnfolded_unwrapped->Draw();
			  // 	c_output->SaveAs("DEBUG_2D_unfolded_unwrapped_" + acceptanceName + ".pdf");
			  // 	hUnfolded->SetTitle("Unfolded ttDil");
			  // 	hUnfolded->GetZaxis()->SetMoreLogLabels();
			  // 	hUnfolded->Draw("COLZ");
			  // 	c_output->SaveAs("DEBUG_2D_unfolded_" + acceptanceName + ".pdf");
			  // }
			  //-------------------------------------------//

			  vector<double> afb2D;
			  vector<double> afb2Derr;
			  afb2D.clear();
			  afb2Derr.clear();
			  GetCorrectedAfb2d(hUnfolded, m_unfoldE, afb2D, afb2Derr, second_output_file);

			  vector<double> afbtrue2D;
			  vector<double> afbtrue2Derr;
			  afbtrue2D.clear();
			  afbtrue2Derr.clear();
			  GetCorrectedAfb2d(hTrue_after, m_unfoldE, afbtrue2D, afbtrue2Derr, second_output_file);
			  //true errors are much smaller (from denominator)
			  afbtrue2Derr[0] = 0.0;
			  afbtrue2Derr[1] = 0.0;
			  afbtrue2Derr[2] = 0.0;
			  afbtrue2Derr[3] = 0.0;

			  
			  // AfbPull[0]->Fill( (Afb - A_gen_k)  / AfbErr );
			  // SumAsym[0] += Afb;
			  // SumErrAsym[0] += AfbErr;
			  // SumTrueAsym[0] += A_gen_k;
			  // SumTrueErrAsym[0] += Aerr_gen_k;

			  for (int iD = 0; iD < nbinsy2D + 1; ++iD)
                {
				  AfbPull[iD]->Fill( (afb2D[iD] - afbtrue2D[iD])  / afb2Derr[iD] );
				  SumAsym[iD] += afb2D[iD];
				  SumErrAsym[iD] += afb2Derr[iD];
				  SumTrueAsym[iD] += afbtrue2D[iD];
				  SumTrueErrAsym[iD] += afbtrue2Derr[iD];
                }

			  unwrap2dhisto(hTrue_after, hTrue_after_unwrapped);
			  unwrap2dhisto(hUnfolded, hUnfolded_unwrapped); //Get 1D version of unfolded results, including the latest error corrections and such

			  TH1D* hUnfolded_unwrapped_clone = (TH1D*) hUnfolded_unwrapped->Clone("hUnfolded_unwrapped_clone");
			  TH1D* hTrue_after_unwrapped_clone = (TH1D*) hTrue_after_unwrapped->Clone("hTrue_after_unwrapped_clone");

			  hUnfolded_unwrapped_clone->Scale(1. / hUnfolded_unwrapped_clone->Integral());
			  hTrue_after_unwrapped_clone->Scale(1. / hTrue_after_unwrapped_clone->Integral());

			  // For filling the bin-by-bin values
			  for (int iD = nbinsy2D + 1; iD < nbinsy2D+nbinsunwrapped_gen+1; ++iD)
                {
				  int binnumber = iD - nbinsy2D;
				  AfbPull[iD]->Fill( (hUnfolded_unwrapped_clone->GetBinContent(binnumber) - hTrue_after_unwrapped_clone->GetBinContent(binnumber)) / hUnfolded_unwrapped_clone->GetBinError(binnumber) );
				  SumAsym[iD] += hUnfolded_unwrapped_clone->GetBinContent(binnumber);
				  SumErrAsym[iD] += hUnfolded_unwrapped_clone->GetBinError(binnumber);
				  SumTrueAsym[iD] += hTrue_after_unwrapped_clone->GetBinContent(binnumber);
				  //SumTrueErrAsym[iD] + = hTrue_after_unwrapped_clone->GetBinError(binnumber);
				  SumTrueErrAsym[iD] += 0.0;
                }


			  /*
                for (int j = 0; j < nbinsunwrapped; j++)
                {
				double pull = (hUnfolded->GetBinContent(j + 1) - hTrue_after->GetBinContent(j + 1)) / hUnfolded->GetBinError(j + 1);
				h_pulls[j]->Fill(pull);
				double resd = (hUnfolded->GetBinContent(j + 1) - hTrue_after->GetBinContent(j + 1)) / hTrue_after->GetBinContent(j + 1);
				h_resd[j]->Fill(resd);
                }
			  */
            } // End loop over all pseudoexperiments

		  //cout << "Average Asymmetry =" << SumAsym / nPseudos << " +/-  " << SumErrAsym / (nPseudos) << "\n";

		  vector<Float_t> temp_A_pull;
		  vector<Float_t> temp_Aerr_pull;
		  vector<Float_t> temp_A_pullwidth;
		  vector<Float_t> temp_Aerr_pullwidth;
		  temp_A_pull.clear();
		  temp_Aerr_pull.clear();
		  temp_A_pullwidth.clear();
		  temp_Aerr_pullwidth.clear();

		  for (int iD = 0; iD < nbinsy2D + nbinsunwrapped_gen + 1; ++iD)
            {
			  temp_A_pull.push_back( AfbPull[iD]->GetMean() );
			  temp_Aerr_pull.push_back( AfbPull[iD]->GetMeanError() );
			  temp_A_pullwidth.push_back( AfbPull[iD]->GetRMS() );
			  temp_Aerr_pullwidth.push_back( AfbPull[iD]->GetRMSError() );
			  SumAsym[iD] /= nPseudos;
			  SumErrAsym[iD] /= (nPseudos * sqrt(nPseudos));
			  SumTrueAsym[iD] /= nPseudos;
			  SumTrueErrAsym[iD] /= (nPseudos * sqrt(nPseudos));
			  if (nPseudos == 1)
                {
				  SumErrAsym[iD] = 0.;
				  SumTrueErrAsym[iD] = 0.;
                }
            }

		  A_unf.push_back( SumAsym );
		  Aerr_unf.push_back( SumErrAsym );
		  A_gen.push_back( SumTrueAsym );
		  Aerr_gen.push_back( SumTrueErrAsym );

		  A_pull.push_back( temp_A_pull );
		  Aerr_pull.push_back( temp_Aerr_pull );
		  A_pullwidth.push_back( temp_A_pullwidth );
		  Aerr_pullwidth.push_back( temp_Aerr_pullwidth );

        }  // end of "if(npseudos>0)"

    } // end of loop over values of k




    if ((TestType == "Linearity"))
    {
        //vector<vector<Float_t>> transposed_A_gen, transposed_Aerr_gen, transposed_A_unf, transposed_Aerr_unf;
        //vector<vector<Float_t>> transposed_A_pull, transposed_A_pullwidth, transposed_Aerr_pull, transposed_Aerr_pullwidth;

        //vector<Float_t> transposed_A_gen, transposed_Aerr_gen, transposed_A_unf, transposed_Aerr_unf;
        //vector<Float_t> transposed_A_pull, transposed_A_pullwidth, transposed_Aerr_pull, transposed_Aerr_pullwidth;
        //
        //transposed_A_gen.clear();
        //transposed_Aerr_gen.clear();
        //transposed_A_unf.clear();
        //transposed_Aerr_unf.clear();
        //transposed_A_pull.clear();
        //transposed_A_pullwidth.clear();
        //transposed_Aerr_pull.clear();
        //transposed_Aerr_pullwidth.clear();

        //cout << " Got to line: " << __LINE__ << endl;

        //for (int iD = 0; iD < nDiff + 1; ++iD)
        //{

        /*
                vector<Float_t> tmp_A_gen, tmp_Aerr_gen, tmp_A_unf, tmp_Aerr_unf;
                vector<Float_t> tmp_A_pull, tmp_A_pullwidth, tmp_Aerr_pull, tmp_Aerr_pullwidth;

                tmp_A_gen.clear();
                tmp_Aerr_gen.clear();
                tmp_A_unf.clear();
                tmp_Aerr_unf.clear();
                tmp_A_pull.clear();
                tmp_A_pullwidth.clear();
                tmp_Aerr_pull.clear();
                tmp_Aerr_pullwidth.clear();
                */

        Float_t inclusive_A_gen[Nlin], inclusive_Aerr_gen[Nlin], inclusive_A_unf[Nlin], inclusive_Aerr_unf[Nlin];
        Float_t inclusive_A_pull[Nlin], inclusive_A_pullwidth[Nlin], inclusive_Aerr_pull[Nlin], inclusive_Aerr_pullwidth[Nlin];

        Float_t bin1_A_gen[Nlin], bin1_Aerr_gen[Nlin], bin1_A_unf[Nlin], bin1_Aerr_unf[Nlin];
        Float_t bin1_A_pull[Nlin], bin1_A_pullwidth[Nlin], bin1_Aerr_pull[Nlin], bin1_Aerr_pullwidth[Nlin];

        Float_t bin2_A_gen[Nlin], bin2_Aerr_gen[Nlin], bin2_A_unf[Nlin], bin2_Aerr_unf[Nlin];
        Float_t bin2_A_pull[Nlin], bin2_A_pullwidth[Nlin], bin2_Aerr_pull[Nlin], bin2_Aerr_pullwidth[Nlin];

        Float_t bin3_A_gen[Nlin], bin3_Aerr_gen[Nlin], bin3_A_unf[Nlin], bin3_Aerr_unf[Nlin];
        Float_t bin3_A_pull[Nlin], bin3_A_pullwidth[Nlin], bin3_Aerr_pull[Nlin], bin3_Aerr_pullwidth[Nlin];

        Float_t bin1_V_gen[Nlin], bin1_Verr_gen[Nlin], bin1_V_unf[Nlin], bin1_Verr_unf[Nlin];
        // Float_t bin1_V_pull[Nlin], bin1_V_pullwidth[Nlin], bin1_Verr_pull[Nlin], bin1_Verr_pullwidth[Nlin];

        Float_t bin2_V_gen[Nlin], bin2_Verr_gen[Nlin], bin2_V_unf[Nlin], bin2_Verr_unf[Nlin];
        // Float_t bin2_V_pull[Nlin], bin2_V_pullwidth[Nlin], bin2_Verr_pull[Nlin], bin2_Verr_pullwidth[Nlin];

        Float_t bin3_V_gen[Nlin], bin3_Verr_gen[Nlin], bin3_V_unf[Nlin], bin3_Verr_unf[Nlin];
        // Float_t bin3_V_pull[Nlin], bin3_V_pullwidth[Nlin], bin3_Verr_pull[Nlin], bin3_Verr_pullwidth[Nlin];

        Float_t bin4_V_gen[Nlin], bin4_Verr_gen[Nlin], bin4_V_unf[Nlin], bin4_Verr_unf[Nlin];
        // Float_t bin4_V_pull[Nlin], bin4_V_pullwidth[Nlin], bin4_Verr_pull[Nlin], bin4_Verr_pullwidth[Nlin];

        Float_t bin5_V_gen[Nlin], bin5_Verr_gen[Nlin], bin5_V_unf[Nlin], bin5_Verr_unf[Nlin];
        // Float_t bin5_V_pull[Nlin], bin5_V_pullwidth[Nlin], bin5_Verr_pull[Nlin], bin5_Verr_pullwidth[Nlin];

        Float_t bin6_V_gen[Nlin], bin6_Verr_gen[Nlin], bin6_V_unf[Nlin], bin6_Verr_unf[Nlin];
        // Float_t bin6_V_pull[Nlin], bin6_V_pullwidth[Nlin], bin6_Verr_pull[Nlin], bin6_Verr_pullwidth[Nlin];

		//Additional items necessitated by the switch from 1D to 2D
        Float_t bin7_V_gen[Nlin], bin7_Verr_gen[Nlin], bin7_V_unf[Nlin], bin7_Verr_unf[Nlin];
        Float_t bin8_V_gen[Nlin], bin8_Verr_gen[Nlin], bin8_V_unf[Nlin], bin8_Verr_unf[Nlin];
        Float_t bin9_V_gen[Nlin], bin9_Verr_gen[Nlin], bin9_V_unf[Nlin], bin9_Verr_unf[Nlin];
        Float_t bin10_V_gen[Nlin], bin10_Verr_gen[Nlin], bin10_V_unf[Nlin], bin10_Verr_unf[Nlin];
        Float_t bin11_V_gen[Nlin], bin11_Verr_gen[Nlin], bin11_V_unf[Nlin], bin11_Verr_unf[Nlin];
        Float_t bin12_V_gen[Nlin], bin12_Verr_gen[Nlin], bin12_V_unf[Nlin], bin12_Verr_unf[Nlin];
        Float_t bin13_V_gen[Nlin], bin13_Verr_gen[Nlin], bin13_V_unf[Nlin], bin13_Verr_unf[Nlin];
        Float_t bin14_V_gen[Nlin], bin14_Verr_gen[Nlin], bin14_V_unf[Nlin], bin14_Verr_unf[Nlin];
        Float_t bin15_V_gen[Nlin], bin15_Verr_gen[Nlin], bin15_V_unf[Nlin], bin15_Verr_unf[Nlin];
        Float_t bin16_V_gen[Nlin], bin16_Verr_gen[Nlin], bin16_V_unf[Nlin], bin16_Verr_unf[Nlin];
        Float_t bin17_V_gen[Nlin], bin17_Verr_gen[Nlin], bin17_V_unf[Nlin], bin17_Verr_unf[Nlin];
        Float_t bin18_V_gen[Nlin], bin18_Verr_gen[Nlin], bin18_V_unf[Nlin], bin18_Verr_unf[Nlin];


        for (int k = 0; k < Nlin; k++)
        {

            inclusive_A_gen[k] = A_gen[k][0];
            inclusive_Aerr_gen[k] = Aerr_gen[k][0];
            inclusive_A_unf[k] = A_unf[k][0];
            inclusive_Aerr_unf[k] = Aerr_unf[k][0];
            inclusive_A_pull[k] = A_pull[k][0];
            inclusive_A_pullwidth[k] = A_pullwidth[k][0];
            inclusive_Aerr_pull[k] = Aerr_pull[k][0];
            inclusive_Aerr_pullwidth[k] = Aerr_pullwidth[k][0];

            bin1_A_gen[k] = A_gen[k][1];
            bin1_Aerr_gen[k] = Aerr_gen[k][1];
            bin1_A_unf[k] = A_unf[k][1];
            bin1_Aerr_unf[k] = Aerr_unf[k][1];
            bin1_A_pull[k] = A_pull[k][1];
            bin1_A_pullwidth[k] = A_pullwidth[k][1];
            bin1_Aerr_pull[k] = Aerr_pull[k][1];
            bin1_Aerr_pullwidth[k] = Aerr_pullwidth[k][1];

            bin2_A_gen[k] = A_gen[k][2];
            bin2_Aerr_gen[k] = Aerr_gen[k][2];
            bin2_A_unf[k] = A_unf[k][2];
            bin2_Aerr_unf[k] = Aerr_unf[k][2];
            bin2_A_pull[k] = A_pull[k][2];
            bin2_A_pullwidth[k] = A_pullwidth[k][2];
            bin2_Aerr_pull[k] = Aerr_pull[k][2];
            bin2_Aerr_pullwidth[k] = Aerr_pullwidth[k][2];

            bin3_A_gen[k] = A_gen[k][3];
            bin3_Aerr_gen[k] = Aerr_gen[k][3];
            bin3_A_unf[k] = A_unf[k][3];
            bin3_Aerr_unf[k] = Aerr_unf[k][3];
            bin3_A_pull[k] = A_pull[k][3];
            bin3_A_pullwidth[k] = A_pullwidth[k][3];
            bin3_Aerr_pull[k] = Aerr_pull[k][3];
            bin3_Aerr_pullwidth[k] = Aerr_pullwidth[k][3];

			int binoffset = nbinsy2D;

            bin1_V_gen[k] = A_gen[k][binoffset+1];
            bin1_Verr_gen[k] = Aerr_gen[k][binoffset+1];
            bin1_V_unf[k] = A_unf[k][binoffset+1];
            bin1_Verr_unf[k] = Aerr_unf[k][binoffset+1];
            // bin1_V_pull[k] = A_pull[k][binoffset+1];
            // bin1_V_pullwidth[k] = A_pullwidth[k][binoffset+1];
            // bin1_Verr_pull[k] = Aerr_pull[k][binoffset+1];
            // bin1_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+1];

            bin2_V_gen[k] = A_gen[k][binoffset+2];
            bin2_Verr_gen[k] = Aerr_gen[k][binoffset+2];
            bin2_V_unf[k] = A_unf[k][binoffset+2];
            bin2_Verr_unf[k] = Aerr_unf[k][binoffset+2];
            // bin2_V_pull[k] = A_pull[k][binoffset+2];
            // bin2_V_pullwidth[k] = A_pullwidth[k][binoffset+2];
            // bin2_Verr_pull[k] = Aerr_pull[k][binoffset+2];
            // bin2_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+2];

            bin3_V_gen[k] = A_gen[k][binoffset+3];
            bin3_Verr_gen[k] = Aerr_gen[k][binoffset+3];
            bin3_V_unf[k] = A_unf[k][binoffset+3];
            bin3_Verr_unf[k] = Aerr_unf[k][binoffset+3];
            // bin3_V_pull[k] = A_pull[k][binoffset+3];
            // bin3_V_pullwidth[k] = A_pullwidth[k][binoffset+3];
            // bin3_Verr_pull[k] = Aerr_pull[k][binoffset+3];
            // bin3_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+3];

            bin4_V_gen[k] = A_gen[k][binoffset+4];
            bin4_Verr_gen[k] = Aerr_gen[k][binoffset+4];
            bin4_V_unf[k] = A_unf[k][binoffset+4];
            bin4_Verr_unf[k] = Aerr_unf[k][binoffset+4];
            // bin4_V_pull[k] = A_pull[k][binoffset+4];
            // bin4_V_pullwidth[k] = A_pullwidth[k][binoffset+4];
            // bin4_Verr_pull[k] = Aerr_pull[k][binoffset+4];
            // bin4_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+4];

            bin5_V_gen[k] = A_gen[k][binoffset+5];
            bin5_Verr_gen[k] = Aerr_gen[k][binoffset+5];
            bin5_V_unf[k] = A_unf[k][binoffset+5];
            bin5_Verr_unf[k] = Aerr_unf[k][binoffset+5];
            // bin5_V_pull[k] = A_pull[k][binoffset+5];
            // bin5_V_pullwidth[k] = A_pullwidth[k][binoffset+5];
            // bin5_Verr_pull[k] = Aerr_pull[k][binoffset+5];
            // bin5_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+5];

            bin6_V_gen[k] = A_gen[k][binoffset+6];
            bin6_Verr_gen[k] = Aerr_gen[k][binoffset+6];
            bin6_V_unf[k] = A_unf[k][binoffset+6];
            bin6_Verr_unf[k] = Aerr_unf[k][binoffset+6];
            // bin6_V_pull[k] = A_pull[k][binoffset+6];
            // bin6_V_pullwidth[k] = A_pullwidth[k][binoffset+6];
            // bin6_Verr_pull[k] = Aerr_pull[k][binoffset+6];
            // bin6_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+6];


			//Additional items necessitated by the switch from 1D to 2D
            bin7_V_gen[k] = A_gen[k][binoffset+7];
            bin7_Verr_gen[k] = Aerr_gen[k][binoffset+7];
            bin7_V_unf[k] = A_unf[k][binoffset+7];
            bin7_Verr_unf[k] = Aerr_unf[k][binoffset+7];
            // bin7_V_pull[k] = A_pull[k][binoffset+7];
            // bin7_V_pullwidth[k] = A_pullwidth[k][binoffset+7];
            // bin7_Verr_pull[k] = Aerr_pull[k][binoffset+7];
            // bin7_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+7];

            bin8_V_gen[k] = A_gen[k][binoffset+8];
            bin8_Verr_gen[k] = Aerr_gen[k][binoffset+8];
            bin8_V_unf[k] = A_unf[k][binoffset+8];
            bin8_Verr_unf[k] = Aerr_unf[k][binoffset+8];
            // bin8_V_pull[k] = A_pull[k][binoffset+8];
            // bin8_V_pullwidth[k] = A_pullwidth[k][binoffset+8];
            // bin8_Verr_pull[k] = Aerr_pull[k][binoffset+8];
            // bin8_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+8];

            bin9_V_gen[k] = A_gen[k][binoffset+9];
            bin9_Verr_gen[k] = Aerr_gen[k][binoffset+9];
            bin9_V_unf[k] = A_unf[k][binoffset+9];
            bin9_Verr_unf[k] = Aerr_unf[k][binoffset+9];
            // bin9_V_pull[k] = A_pull[k][binoffset+9];
            // bin9_V_pullwidth[k] = A_pullwidth[k][binoffset+9];
            // bin9_Verr_pull[k] = Aerr_pull[k][binoffset+9];
            // bin9_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+9];

            bin10_V_gen[k] = A_gen[k][binoffset+10];
            bin10_Verr_gen[k] = Aerr_gen[k][binoffset+10];
            bin10_V_unf[k] = A_unf[k][binoffset+10];
            bin10_Verr_unf[k] = Aerr_unf[k][binoffset+10];
            // bin10_V_pull[k] = A_pull[k][binoffset+10];
            // bin10_V_pullwidth[k] = A_pullwidth[k][binoffset+10];
            // bin10_Verr_pull[k] = Aerr_pull[k][binoffset+10];
            // bin10_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+10];

            bin11_V_gen[k] = A_gen[k][binoffset+11];
            bin11_Verr_gen[k] = Aerr_gen[k][binoffset+11];
            bin11_V_unf[k] = A_unf[k][binoffset+11];
            bin11_Verr_unf[k] = Aerr_unf[k][binoffset+11];
            // bin11_V_pull[k] = A_pull[k][binoffset+11];
            // bin11_V_pullwidth[k] = A_pullwidth[k][binoffset+11];
            // bin11_Verr_pull[k] = Aerr_pull[k][binoffset+11];
            // bin11_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+11];

            bin12_V_gen[k] = A_gen[k][binoffset+12];
            bin12_Verr_gen[k] = Aerr_gen[k][binoffset+12];
            bin12_V_unf[k] = A_unf[k][binoffset+12];
            bin12_Verr_unf[k] = Aerr_unf[k][binoffset+12];
            // bin12_V_pull[k] = A_pull[k][binoffset+12];
            // bin12_V_pullwidth[k] = A_pullwidth[k][binoffset+12];
            // bin12_Verr_pull[k] = Aerr_pull[k][binoffset+12];
            // bin12_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+12];

            bin13_V_gen[k] = A_gen[k][binoffset+13];
            bin13_Verr_gen[k] = Aerr_gen[k][binoffset+13];
            bin13_V_unf[k] = A_unf[k][binoffset+13];
            bin13_Verr_unf[k] = Aerr_unf[k][binoffset+13];
            // bin13_V_pull[k] = A_pull[k][binoffset+13];
            // bin13_V_pullwidth[k] = A_pullwidth[k][binoffset+13];
            // bin13_Verr_pull[k] = Aerr_pull[k][binoffset+13];
            // bin13_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+13];

            bin14_V_gen[k] = A_gen[k][binoffset+14];
            bin14_Verr_gen[k] = Aerr_gen[k][binoffset+14];
            bin14_V_unf[k] = A_unf[k][binoffset+14];
            bin14_Verr_unf[k] = Aerr_unf[k][binoffset+14];
            // bin14_V_pull[k] = A_pull[k][binoffset+14];
            // bin14_V_pullwidth[k] = A_pullwidth[k][binoffset+14];
            // bin14_Verr_pull[k] = Aerr_pull[k][binoffset+14];
            // bin14_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+14];

            bin15_V_gen[k] = A_gen[k][binoffset+15];
            bin15_Verr_gen[k] = Aerr_gen[k][binoffset+15];
            bin15_V_unf[k] = A_unf[k][binoffset+15];
            bin15_Verr_unf[k] = Aerr_unf[k][binoffset+15];
            // bin15_V_pull[k] = A_pull[k][binoffset+15];
            // bin15_V_pullwidth[k] = A_pullwidth[k][binoffset+15];
            // bin15_Verr_pull[k] = Aerr_pull[k][binoffset+15];
            // bin15_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+15];

            bin16_V_gen[k] = A_gen[k][binoffset+16];
            bin16_Verr_gen[k] = Aerr_gen[k][binoffset+16];
            bin16_V_unf[k] = A_unf[k][binoffset+16];
            bin16_Verr_unf[k] = Aerr_unf[k][binoffset+16];
            // bin16_V_pull[k] = A_pull[k][binoffset+16];
            // bin16_V_pullwidth[k] = A_pullwidth[k][binoffset+16];
            // bin16_Verr_pull[k] = Aerr_pull[k][binoffset+16];
            // bin16_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+16];

            bin17_V_gen[k] = A_gen[k][binoffset+17];
            bin17_Verr_gen[k] = Aerr_gen[k][binoffset+17];
            bin17_V_unf[k] = A_unf[k][binoffset+17];
            bin17_Verr_unf[k] = Aerr_unf[k][binoffset+17];
            // bin17_V_pull[k] = A_pull[k][binoffset+17];
            // bin17_V_pullwidth[k] = A_pullwidth[k][binoffset+17];
            // bin17_Verr_pull[k] = Aerr_pull[k][binoffset+17];
            // bin17_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+17];

            bin18_V_gen[k] = A_gen[k][binoffset+18];
            bin18_Verr_gen[k] = Aerr_gen[k][binoffset+18];
            bin18_V_unf[k] = A_unf[k][binoffset+18];
            bin18_Verr_unf[k] = Aerr_unf[k][binoffset+18];
            // bin18_V_pull[k] = A_pull[k][binoffset+18];
            // bin18_V_pullwidth[k] = A_pullwidth[k][binoffset+18];
            // bin18_Verr_pull[k] = Aerr_pull[k][binoffset+18];
            // bin18_Verr_pullwidth[k] = Aerr_pullwidth[k][binoffset+18];


            //cout << A_gen[k][0] << " " << Aerr_gen[k][0] << " " << A_unf[k][0] << " " << Aerr_unf[k][0] << " " << A_pull[k][0] << " " << A_pullwidth[k][0] << " " << Aerr_pull[k][0] << " " << Aerr_pullwidth[k][0] << endl;
            //cout << A_gen[k][1] << " " << Aerr_gen[k][1] << " " << A_unf[k][1] << " " << Aerr_unf[k][1] << " " << A_pull[k][1] << " " << A_pullwidth[k][1] << " " << Aerr_pull[k][1] << " " << Aerr_pullwidth[k][1] << endl;
            //cout << A_gen[k][2] << " " << Aerr_gen[k][2] << " " << A_unf[k][2] << " " << Aerr_unf[k][2] << " " << A_pull[k][2] << " " << A_pullwidth[k][2] << " " << Aerr_pull[k][2] << " " << Aerr_pullwidth[k][2] << endl;
            //cout << A_gen[k][3] << " " << Aerr_gen[k][3] << " " << A_unf[k][3] << " " << Aerr_unf[k][3] << " " << A_pull[k][3] << " " << A_pullwidth[k][3] << " " << Aerr_pull[k][3] << " " << Aerr_pullwidth[k][3] << endl;
            //cout << A_gen[k][4] << " " << Aerr_gen[k][4] << " " << A_unf[k][4] << " " << Aerr_unf[k][4] << " " << A_pull[k][4] << " " << A_pullwidth[k][4] << " " << Aerr_pull[k][4] << " " << Aerr_pullwidth[k][4] << endl;
            //cout << A_gen[k][5] << " " << Aerr_gen[k][5] << " " << A_unf[k][5] << " " << Aerr_unf[k][5] << " " << A_pull[k][5] << " " << A_pullwidth[k][5] << " " << Aerr_pull[k][5] << " " << Aerr_pullwidth[k][5] << endl;
            //cout << A_gen[k][6] << " " << Aerr_gen[k][6] << " " << A_unf[k][6] << " " << Aerr_unf[k][6] << " " << A_pull[k][6] << " " << A_pullwidth[k][6] << " " << Aerr_pull[k][6] << " " << Aerr_pullwidth[k][6] << endl;
            //cout << A_gen[k][7] << " " << Aerr_gen[k][7] << " " << A_unf[k][7] << " " << Aerr_unf[k][7] << " " << A_pull[k][7] << " " << A_pullwidth[k][7] << " " << Aerr_pull[k][7] << " " << Aerr_pullwidth[k][7] << endl;
            //cout << A_gen[k][8] << " " << Aerr_gen[k][8] << " " << A_unf[k][8] << " " << Aerr_unf[k][8] << " " << A_pull[k][8] << " " << A_pullwidth[k][8] << " " << Aerr_pull[k][8] << " " << Aerr_pullwidth[k][8] << endl;
            //cout << A_gen[k][9] << " " << Aerr_gen[k][9] << " " << A_unf[k][9] << " " << Aerr_unf[k][9] << " " << A_pull[k][9] << " " << A_pullwidth[k][9] << " " << Aerr_pull[k][9] << " " << Aerr_pullwidth[k][9] << endl;




            /*
                        tmp_A_gen.push_back( A_gen[k][iD] );
                        tmp_Aerr_gen.push_back( Aerr_gen[k][iD] );
                        tmp_A_unf.push_back( A_unf[k][iD] );
                        tmp_Aerr_unf.push_back( Aerr_unf[k][iD] );
                        tmp_A_pull.push_back( A_pull[k][iD] );
                        tmp_A_pullwidth.push_back( A_pullwidth[k][iD] );
                        tmp_Aerr_pull.push_back( Aerr_pull[k][iD] );
                        tmp_Aerr_pullwidth.push_back( Aerr_pullwidth[k][iD] );
            */

        }

        //transposed_A_gen.push_back( tmp_A_gen );
        //transposed_Aerr_gen.push_back( tmp_Aerr_gen );
        //transposed_A_unf.push_back( tmp_A_unf );
        //transposed_Aerr_unf.push_back( tmp_Aerr_unf );
        //transposed_A_pull.push_back( tmp_A_pull );
        //transposed_A_pullwidth.push_back( tmp_A_pullwidth );
        //transposed_Aerr_pull.push_back( tmp_Aerr_pull );
        //transposed_Aerr_pullwidth.push_back( tmp_Aerr_pullwidth );

        //}



        TGraphErrors *Asym2D_TrueUnf = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_unf, inclusive_Aerr_gen, inclusive_Aerr_unf);

        //TGraphErrors *Asym2D_TrueMeas = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_meas, inclusive_Aerr_gen, inclusive_Aerr_meas);

        TGraphErrors *Asym2D_PullWidth = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_pullwidth, inclusive_Aerr_gen, inclusive_Aerr_pullwidth);

        TGraphErrors *Asym2D_Pull = new TGraphErrors (Nlin, inclusive_A_gen, inclusive_A_pull, inclusive_Aerr_gen, inclusive_Aerr_pull);





        TGraphErrors *Asym2D_TrueUnfbin1 = new TGraphErrors (Nlin, bin1_A_gen, bin1_A_unf, bin1_Aerr_gen, bin1_Aerr_unf);
        TGraphErrors *Asym2D_PullWidthbin1 = new TGraphErrors (Nlin, bin1_A_gen, bin1_A_pullwidth, bin1_Aerr_gen, bin1_Aerr_pullwidth);
        TGraphErrors *Asym2D_Pullbin1 = new TGraphErrors (Nlin, bin1_A_gen, bin1_A_pull, bin1_Aerr_gen, bin1_Aerr_pull);

        TGraphErrors *Asym2D_TrueUnfbin2 = new TGraphErrors (Nlin, bin2_A_gen, bin2_A_unf, bin2_Aerr_gen, bin2_Aerr_unf);
        TGraphErrors *Asym2D_PullWidthbin2 = new TGraphErrors (Nlin, bin2_A_gen, bin2_A_pullwidth, bin2_Aerr_gen, bin2_Aerr_pullwidth);
        TGraphErrors *Asym2D_Pullbin2 = new TGraphErrors (Nlin, bin2_A_gen, bin2_A_pull, bin2_Aerr_gen, bin2_Aerr_pull);

        TGraphErrors *Asym2D_TrueUnfbin3 = new TGraphErrors (Nlin, bin3_A_gen, bin3_A_unf, bin3_Aerr_gen, bin3_Aerr_unf);
        TGraphErrors *Asym2D_PullWidthbin3 = new TGraphErrors (Nlin, bin3_A_gen, bin3_A_pullwidth, bin3_Aerr_gen, bin3_Aerr_pullwidth);
        TGraphErrors *Asym2D_Pullbin3 = new TGraphErrors (Nlin, bin3_A_gen, bin3_A_pull, bin3_Aerr_gen, bin3_Aerr_pull);





        TGraphErrors *Value_TrueUnfbin1 = new TGraphErrors (Nlin, bin1_V_gen, bin1_V_unf, bin1_Verr_gen, bin1_Verr_unf);
        // TGraphErrors *Value_PullWidthbin1 = new TGraphErrors (Nlin, bin1_V_gen, bin1_V_pullwidth, bin1_Verr_gen, bin1_Verr_pullwidth);
        // TGraphErrors *Value_Pullbin1 = new TGraphErrors (Nlin, bin1_V_gen, bin1_V_pull, bin1_Verr_gen, bin1_Verr_pull);

        TGraphErrors *Value_TrueUnfbin2 = new TGraphErrors (Nlin, bin2_V_gen, bin2_V_unf, bin2_Verr_gen, bin2_Verr_unf);
        // TGraphErrors *Value_PullWidthbin2 = new TGraphErrors (Nlin, bin2_V_gen, bin2_V_pullwidth, bin2_Verr_gen, bin2_Verr_pullwidth);
        // TGraphErrors *Value_Pullbin2 = new TGraphErrors (Nlin, bin2_V_gen, bin2_V_pull, bin2_Verr_gen, bin2_Verr_pull);

        TGraphErrors *Value_TrueUnfbin3 = new TGraphErrors (Nlin, bin3_V_gen, bin3_V_unf, bin3_Verr_gen, bin3_Verr_unf);
        // TGraphErrors *Value_PullWidthbin3 = new TGraphErrors (Nlin, bin3_V_gen, bin3_V_pullwidth, bin3_Verr_gen, bin3_Verr_pullwidth);
        // TGraphErrors *Value_Pullbin3 = new TGraphErrors (Nlin, bin3_V_gen, bin3_V_pull, bin3_Verr_gen, bin3_Verr_pull);

        TGraphErrors *Value_TrueUnfbin4 = new TGraphErrors (Nlin, bin4_V_gen, bin4_V_unf, bin4_Verr_gen, bin4_Verr_unf);
        // TGraphErrors *Value_PullWidthbin4 = new TGraphErrors (Nlin, bin4_V_gen, bin4_V_pullwidth, bin4_Verr_gen, bin4_Verr_pullwidth);
        // TGraphErrors *Value_Pullbin4 = new TGraphErrors (Nlin, bin4_V_gen, bin4_V_pull, bin4_Verr_gen, bin4_Verr_pull);

        TGraphErrors *Value_TrueUnfbin5 = new TGraphErrors (Nlin, bin5_V_gen, bin5_V_unf, bin5_Verr_gen, bin5_Verr_unf);
        // TGraphErrors *Value_PullWidthbin5 = new TGraphErrors (Nlin, bin5_V_gen, bin5_V_pullwidth, bin5_Verr_gen, bin5_Verr_pullwidth);
        // TGraphErrors *Value_Pullbin5 = new TGraphErrors (Nlin, bin5_V_gen, bin5_V_pull, bin5_Verr_gen, bin5_Verr_pull);

        TGraphErrors *Value_TrueUnfbin6 = new TGraphErrors (Nlin, bin6_V_gen, bin6_V_unf, bin6_Verr_gen, bin6_Verr_unf);
        // TGraphErrors *Value_PullWidthbin6 = new TGraphErrors (Nlin, bin6_V_gen, bin6_V_pullwidth, bin6_Verr_gen, bin6_Verr_pullwidth);
        // TGraphErrors *Value_Pullbin6 = new TGraphErrors (Nlin, bin6_V_gen, bin6_V_pull, bin6_Verr_gen, bin6_Verr_pull);

		//Additional plots necessitated by the move from 1D to 2D
		// At the time of last revision, the graphs "Value_PullWidthbin#" and "Value_Pullbin#" aren't used, so they're not expanded up to 18.
        TGraphErrors *Value_TrueUnfbin7 = new TGraphErrors (Nlin, bin7_V_gen, bin7_V_unf, bin7_Verr_gen, bin7_Verr_unf);
        TGraphErrors *Value_TrueUnfbin8 = new TGraphErrors (Nlin, bin8_V_gen, bin8_V_unf, bin8_Verr_gen, bin8_Verr_unf);
        TGraphErrors *Value_TrueUnfbin9 = new TGraphErrors (Nlin, bin9_V_gen, bin9_V_unf, bin9_Verr_gen, bin9_Verr_unf);
        TGraphErrors *Value_TrueUnfbin10 = new TGraphErrors (Nlin, bin10_V_gen, bin10_V_unf, bin10_Verr_gen, bin10_Verr_unf);
        TGraphErrors *Value_TrueUnfbin11 = new TGraphErrors (Nlin, bin11_V_gen, bin11_V_unf, bin11_Verr_gen, bin11_Verr_unf);
        TGraphErrors *Value_TrueUnfbin12 = new TGraphErrors (Nlin, bin12_V_gen, bin12_V_unf, bin12_Verr_gen, bin12_Verr_unf);
        TGraphErrors *Value_TrueUnfbin13 = new TGraphErrors (Nlin, bin13_V_gen, bin13_V_unf, bin13_Verr_gen, bin13_Verr_unf);
        TGraphErrors *Value_TrueUnfbin14 = new TGraphErrors (Nlin, bin14_V_gen, bin14_V_unf, bin14_Verr_gen, bin14_Verr_unf);
        TGraphErrors *Value_TrueUnfbin15 = new TGraphErrors (Nlin, bin15_V_gen, bin15_V_unf, bin15_Verr_gen, bin15_Verr_unf);
        TGraphErrors *Value_TrueUnfbin16 = new TGraphErrors (Nlin, bin16_V_gen, bin16_V_unf, bin16_Verr_gen, bin16_Verr_unf);
        TGraphErrors *Value_TrueUnfbin17 = new TGraphErrors (Nlin, bin17_V_gen, bin17_V_unf, bin17_Verr_gen, bin17_Verr_unf);
        TGraphErrors *Value_TrueUnfbin18 = new TGraphErrors (Nlin, bin18_V_gen, bin18_V_unf, bin18_Verr_gen, bin18_Verr_unf);



        TCanvas *c_ttbar = new TCanvas("c_ttbar", "c_ttbar", 500, 500);
        if (!plot_inclusive_only) c_ttbar->Divide(2, 2);
        if (!plot_inclusive_only) c_ttbar->cd(1);
        Asym2D_TrueUnf->SetTitle(asymlabel);
        Asym2D_TrueUnf->SetMarkerStyle(23);
        Asym2D_TrueUnf->SetMarkerColor(kBlack);
        Asym2D_TrueUnf->SetMarkerSize(0.6);
        Asym2D_TrueUnf->GetXaxis()->SetTitle(asymlabel + " inclusive (true)");
        Asym2D_TrueUnf->GetYaxis()->SetTitle(asymlabel + " inclusive (unfolded)");
        Asym2D_TrueUnf->Draw("AP same");
        //Asym2D_TrueUnf->Fit("pol1");
        TFitResultPtr r = Asym2D_TrueUnf->Fit("pol1", "S");
        Double_t par1   = r->Parameter(1);
        Double_t par0   = r->Parameter(0);
        Double_t par1err   = r->ParError(1);
        Double_t par0err   = r->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0: " << par0 << " +/- " << par0err << " p1: " << par1 << " +/- " << par1err <<  endl;


        if (!plot_inclusive_only)
        {
            c_ttbar->cd(2);
            Asym2D_TrueUnfbin1->SetTitle(asymlabel);
            Asym2D_TrueUnfbin1->SetMarkerStyle(23);
            Asym2D_TrueUnfbin1->SetMarkerColor(kBlue);
            Asym2D_TrueUnfbin1->SetMarkerSize(0.6);
            Asym2D_TrueUnfbin1->GetXaxis()->SetTitle(yaxislabel + " bin 1 (true)");
            Asym2D_TrueUnfbin1->GetYaxis()->SetTitle(yaxislabel + " bin 1 (unfolded)");
            Asym2D_TrueUnfbin1->Draw("AP same");
            //Asym2D_TrueUnfbin1->Fit("pol1");
            TFitResultPtr rbin1 = Asym2D_TrueUnfbin1->Fit("pol1", "S");
            par1   = rbin1->Parameter(1);
            par0   = rbin1->Parameter(0);
            par1err   = rbin1->ParError(1);
            par0err   = rbin1->ParError(0);
            third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0bin1: " << par0 << " +/- " << par0err << " p1bin1: " << par1 << " +/- " << par1err <<  endl;

            c_ttbar->cd(3);
            Asym2D_TrueUnfbin2->SetTitle(asymlabel);
            Asym2D_TrueUnfbin2->SetMarkerStyle(23);
            Asym2D_TrueUnfbin2->SetMarkerColor(kBlue);
            Asym2D_TrueUnfbin2->SetMarkerSize(0.6);
            Asym2D_TrueUnfbin2->GetXaxis()->SetTitle(yaxislabel + " bin 2 (true)");
            Asym2D_TrueUnfbin2->GetYaxis()->SetTitle(yaxislabel + " bin 2 (unfolded)");
            Asym2D_TrueUnfbin2->Draw("AP same");
            //Asym2D_TrueUnfbin2->Fit("pol1");
            TFitResultPtr rbin2 = Asym2D_TrueUnfbin2->Fit("pol1", "S");
            par1   = rbin2->Parameter(1);
            par0   = rbin2->Parameter(0);
            par1err   = rbin2->ParError(1);
            par0err   = rbin2->ParError(0);
            third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0bin2: " << par0 << " +/- " << par0err << " p1bin2: " << par1 << " +/- " << par1err <<  endl;

            c_ttbar->cd(4);
            Asym2D_TrueUnfbin3->SetTitle(asymlabel);
            Asym2D_TrueUnfbin3->SetMarkerStyle(23);
            Asym2D_TrueUnfbin3->SetMarkerColor(kBlue);
            Asym2D_TrueUnfbin3->SetMarkerSize(0.6);
            Asym2D_TrueUnfbin3->GetXaxis()->SetTitle(yaxislabel + " bin 3 (true)");
            Asym2D_TrueUnfbin3->GetYaxis()->SetTitle(yaxislabel + " bin 3 (unfolded)");
            Asym2D_TrueUnfbin3->Draw("AP same");
            //Asym2D_TrueUnfbin3->Fit("pol1");
            TFitResultPtr rbin3 = Asym2D_TrueUnfbin3->Fit("pol1", "S");
            par1   = rbin3->Parameter(1);
            par0   = rbin3->Parameter(0);
            par1err   = rbin3->ParError(1);
            par0err   = rbin3->ParError(0);
            third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0bin3: " << par0 << " +/- " << par0err << " p1bin3: " << par1 << " +/- " << par1err <<  endl;
        }

        c_ttbar->SaveAs("2D_" + acceptanceName + "_LinearityCheck.pdf");
        // c_ttbar->SaveAs(acceptanceName + "_LinearityCheck.C");












        TCanvas *c_LinearityBinByBin = new TCanvas("c_LinearityBinByBin", "c_LinearityBinByBin", 1500, 750);
        c_LinearityBinByBin->Divide(6, 3);

        c_LinearityBinByBin->cd(13);
        Value_TrueUnfbin1->SetTitle(asymlabel);
        Value_TrueUnfbin1->SetMarkerStyle(23);
        Value_TrueUnfbin1->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin1->SetMarkerSize(0.6);
        Value_TrueUnfbin1->GetXaxis()->SetTitle(asymlabel + " bin1 (true)");
        Value_TrueUnfbin1->GetYaxis()->SetTitle(asymlabel + " bin1 (unfolded)");
        Value_TrueUnfbin1->Draw("AP same");
        //Value_TrueUnfbin1->Fit("pol1");
        TFitResultPtr rVbin1 = Value_TrueUnfbin1->Fit("pol1", "S");
        par1   = rVbin1->Parameter(1);
        par0   = rVbin1->Parameter(0);
        par1err   = rVbin1->ParError(1);
        par0err   = rVbin1->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin1: " << par0 << " +/- " << par0err << " p1Vbin1: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(14);
        Value_TrueUnfbin2->SetTitle(asymlabel);
        Value_TrueUnfbin2->SetMarkerStyle(23);
        Value_TrueUnfbin2->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin2->SetMarkerSize(0.6);
        Value_TrueUnfbin2->GetXaxis()->SetTitle(asymlabel + " bin2 (true)");
        Value_TrueUnfbin2->GetYaxis()->SetTitle(asymlabel + " bin2 (unfolded)");
        Value_TrueUnfbin2->Draw("AP same");
        //Value_TrueUnfbin2->Fit("pol1");
        TFitResultPtr rVbin2 = Value_TrueUnfbin2->Fit("pol1", "S");
        par1   = rVbin2->Parameter(1);
        par0   = rVbin2->Parameter(0);
        par1err   = rVbin2->ParError(1);
        par0err   = rVbin2->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin2: " << par0 << " +/- " << par0err << " p1Vbin2: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(15);
        Value_TrueUnfbin3->SetTitle(asymlabel);
        Value_TrueUnfbin3->SetMarkerStyle(23);
        Value_TrueUnfbin3->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin3->SetMarkerSize(0.6);
        Value_TrueUnfbin3->GetXaxis()->SetTitle(asymlabel + " bin3 (true)");
        Value_TrueUnfbin3->GetYaxis()->SetTitle(asymlabel + " bin3 (unfolded)");
        Value_TrueUnfbin3->Draw("AP same");
        //Value_TrueUnfbin3->Fit("pol1");
        TFitResultPtr rVbin3 = Value_TrueUnfbin3->Fit("pol1", "S");
        par1   = rVbin3->Parameter(1);
        par0   = rVbin3->Parameter(0);
        par1err   = rVbin3->ParError(1);
        par0err   = rVbin3->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin3: " << par0 << " +/- " << par0err << " p1Vbin3: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(16);
        Value_TrueUnfbin4->SetTitle(asymlabel);
        Value_TrueUnfbin4->SetMarkerStyle(23);
        Value_TrueUnfbin4->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin4->SetMarkerSize(0.6);
        Value_TrueUnfbin4->GetXaxis()->SetTitle(asymlabel + " bin4 (true)");
        Value_TrueUnfbin4->GetYaxis()->SetTitle(asymlabel + " bin4 (unfolded)");
        Value_TrueUnfbin4->Draw("AP same");
        //Value_TrueUnfbin4->Fit("pol1");
        TFitResultPtr rVbin4 = Value_TrueUnfbin4->Fit("pol1", "S");
        par1   = rVbin4->Parameter(1);
        par0   = rVbin4->Parameter(0);
        par1err   = rVbin4->ParError(1);
        par0err   = rVbin4->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin4: " << par0 << " +/- " << par0err << " p1Vbin4: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(17);
        Value_TrueUnfbin5->SetTitle(asymlabel);
        Value_TrueUnfbin5->SetMarkerStyle(23);
        Value_TrueUnfbin5->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin5->SetMarkerSize(0.6);
        Value_TrueUnfbin5->GetXaxis()->SetTitle(asymlabel + " bin5 (true)");
        Value_TrueUnfbin5->GetYaxis()->SetTitle(asymlabel + " bin5 (unfolded)");
        Value_TrueUnfbin5->Draw("AP same");
        //Value_TrueUnfbin5->Fit("pol1");
        TFitResultPtr rVbin5 = Value_TrueUnfbin5->Fit("pol1", "S");
        par1   = rVbin5->Parameter(1);
        par0   = rVbin5->Parameter(0);
        par1err   = rVbin5->ParError(1);
        par0err   = rVbin5->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin5: " << par0 << " +/- " << par0err << " p1Vbin5: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(18);
        Value_TrueUnfbin6->SetTitle(asymlabel);
        Value_TrueUnfbin6->SetMarkerStyle(23);
        Value_TrueUnfbin6->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin6->SetMarkerSize(0.6);
        Value_TrueUnfbin6->GetXaxis()->SetTitle(asymlabel + " bin6 (true)");
        Value_TrueUnfbin6->GetYaxis()->SetTitle(asymlabel + " bin6 (unfolded)");
        Value_TrueUnfbin6->Draw("AP same");
        //Value_TrueUnfbin6->Fit("pol1");
        TFitResultPtr rVbin6 = Value_TrueUnfbin6->Fit("pol1", "S");
        par1   = rVbin6->Parameter(1);
        par0   = rVbin6->Parameter(0);
        par1err   = rVbin6->ParError(1);
        par0err   = rVbin6->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin6: " << par0 << " +/- " << par0err << " p1Vbin6: " << par1 << " +/- " << par1err <<  endl;

		// Additions due to the change from 1D to 2D
        c_LinearityBinByBin->cd(7);
        Value_TrueUnfbin7->SetTitle(asymlabel);
        Value_TrueUnfbin7->SetMarkerStyle(23);
        Value_TrueUnfbin7->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin7->SetMarkerSize(0.6);
        Value_TrueUnfbin7->GetXaxis()->SetTitle(asymlabel + " bin7 (true)");
        Value_TrueUnfbin7->GetYaxis()->SetTitle(asymlabel + " bin7 (unfolded)");
        Value_TrueUnfbin7->Draw("AP same");
        //Value_TrueUnfbin7->Fit("pol1");
        TFitResultPtr rVbin7 = Value_TrueUnfbin7->Fit("pol1", "S");
        par1   = rVbin7->Parameter(1);
        par0   = rVbin7->Parameter(0);
        par1err   = rVbin7->ParError(1);
        par0err   = rVbin7->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin7: " << par0 << " +/- " << par0err << " p1Vbin7: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(8);
        Value_TrueUnfbin8->SetTitle(asymlabel);
        Value_TrueUnfbin8->SetMarkerStyle(23);
        Value_TrueUnfbin8->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin8->SetMarkerSize(0.6);
        Value_TrueUnfbin8->GetXaxis()->SetTitle(asymlabel + " bin8 (true)");
        Value_TrueUnfbin8->GetYaxis()->SetTitle(asymlabel + " bin8 (unfolded)");
        Value_TrueUnfbin8->Draw("AP same");
        //Value_TrueUnfbin8->Fit("pol1");
        TFitResultPtr rVbin8 = Value_TrueUnfbin8->Fit("pol1", "S");
        par1   = rVbin8->Parameter(1);
        par0   = rVbin8->Parameter(0);
        par1err   = rVbin8->ParError(1);
        par0err   = rVbin8->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin8: " << par0 << " +/- " << par0err << " p1Vbin8: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(9);
        Value_TrueUnfbin9->SetTitle(asymlabel);
        Value_TrueUnfbin9->SetMarkerStyle(23);
        Value_TrueUnfbin9->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin9->SetMarkerSize(0.6);
        Value_TrueUnfbin9->GetXaxis()->SetTitle(asymlabel + " bin9 (true)");
        Value_TrueUnfbin9->GetYaxis()->SetTitle(asymlabel + " bin9 (unfolded)");
        Value_TrueUnfbin9->Draw("AP same");
        //Value_TrueUnfbin9->Fit("pol1");
        TFitResultPtr rVbin9 = Value_TrueUnfbin9->Fit("pol1", "S");
        par1   = rVbin9->Parameter(1);
        par0   = rVbin9->Parameter(0);
        par1err   = rVbin9->ParError(1);
        par0err   = rVbin9->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin9: " << par0 << " +/- " << par0err << " p1Vbin9: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(10);
        Value_TrueUnfbin10->SetTitle(asymlabel);
        Value_TrueUnfbin10->SetMarkerStyle(23);
        Value_TrueUnfbin10->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin10->SetMarkerSize(0.6);
        Value_TrueUnfbin10->GetXaxis()->SetTitle(asymlabel + " bin10 (true)");
        Value_TrueUnfbin10->GetYaxis()->SetTitle(asymlabel + " bin10 (unfolded)");
        Value_TrueUnfbin10->Draw("AP same");
        //Value_TrueUnfbin10->Fit("pol1");
        TFitResultPtr rVbin10 = Value_TrueUnfbin10->Fit("pol1", "S");
        par1   = rVbin10->Parameter(1);
        par0   = rVbin10->Parameter(0);
        par1err   = rVbin10->ParError(1);
        par0err   = rVbin10->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin10: " << par0 << " +/- " << par0err << " p1Vbin10: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(11);
        Value_TrueUnfbin11->SetTitle(asymlabel);
        Value_TrueUnfbin11->SetMarkerStyle(23);
        Value_TrueUnfbin11->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin11->SetMarkerSize(0.6);
        Value_TrueUnfbin11->GetXaxis()->SetTitle(asymlabel + " bin11 (true)");
        Value_TrueUnfbin11->GetYaxis()->SetTitle(asymlabel + " bin11 (unfolded)");
        Value_TrueUnfbin11->Draw("AP same");
        //Value_TrueUnfbin11->Fit("pol1");
        TFitResultPtr rVbin11 = Value_TrueUnfbin11->Fit("pol1", "S");
        par1   = rVbin11->Parameter(1);
        par0   = rVbin11->Parameter(0);
        par1err   = rVbin11->ParError(1);
        par0err   = rVbin11->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin11: " << par0 << " +/- " << par0err << " p1Vbin11: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(12);
        Value_TrueUnfbin12->SetTitle(asymlabel);
        Value_TrueUnfbin12->SetMarkerStyle(23);
        Value_TrueUnfbin12->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin12->SetMarkerSize(0.6);
        Value_TrueUnfbin12->GetXaxis()->SetTitle(asymlabel + " bin12 (true)");
        Value_TrueUnfbin12->GetYaxis()->SetTitle(asymlabel + " bin12 (unfolded)");
        Value_TrueUnfbin12->Draw("AP same");
        //Value_TrueUnfbin12->Fit("pol1");
        TFitResultPtr rVbin12 = Value_TrueUnfbin12->Fit("pol1", "S");
        par1   = rVbin12->Parameter(1);
        par0   = rVbin12->Parameter(0);
        par1err   = rVbin12->ParError(1);
        par0err   = rVbin12->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin12: " << par0 << " +/- " << par0err << " p1Vbin12: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(1);
        Value_TrueUnfbin13->SetTitle(asymlabel);
        Value_TrueUnfbin13->SetMarkerStyle(23);
        Value_TrueUnfbin13->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin13->SetMarkerSize(0.6);
        Value_TrueUnfbin13->GetXaxis()->SetTitle(asymlabel + " bin13 (true)");
        Value_TrueUnfbin13->GetYaxis()->SetTitle(asymlabel + " bin13 (unfolded)");
        Value_TrueUnfbin13->Draw("AP same");
        //Value_TrueUnfbin13->Fit("pol1");
        TFitResultPtr rVbin13 = Value_TrueUnfbin13->Fit("pol1", "S");
        par1   = rVbin13->Parameter(1);
        par0   = rVbin13->Parameter(0);
        par1err   = rVbin13->ParError(1);
        par0err   = rVbin13->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin13: " << par0 << " +/- " << par0err << " p1Vbin13: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(2);
        Value_TrueUnfbin14->SetTitle(asymlabel);
        Value_TrueUnfbin14->SetMarkerStyle(23);
        Value_TrueUnfbin14->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin14->SetMarkerSize(0.6);
        Value_TrueUnfbin14->GetXaxis()->SetTitle(asymlabel + " bin14 (true)");
        Value_TrueUnfbin14->GetYaxis()->SetTitle(asymlabel + " bin14 (unfolded)");
        Value_TrueUnfbin14->Draw("AP same");
        //Value_TrueUnfbin14->Fit("pol1");
        TFitResultPtr rVbin14 = Value_TrueUnfbin14->Fit("pol1", "S");
        par1   = rVbin14->Parameter(1);
        par0   = rVbin14->Parameter(0);
        par1err   = rVbin14->ParError(1);
        par0err   = rVbin14->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin14: " << par0 << " +/- " << par0err << " p1Vbin14: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(3);
        Value_TrueUnfbin15->SetTitle(asymlabel);
        Value_TrueUnfbin15->SetMarkerStyle(23);
        Value_TrueUnfbin15->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin15->SetMarkerSize(0.6);
        Value_TrueUnfbin15->GetXaxis()->SetTitle(asymlabel + " bin15 (true)");
        Value_TrueUnfbin15->GetYaxis()->SetTitle(asymlabel + " bin15 (unfolded)");
        Value_TrueUnfbin15->Draw("AP same");
        //Value_TrueUnfbin15->Fit("pol1");
        TFitResultPtr rVbin15 = Value_TrueUnfbin15->Fit("pol1", "S");
        par1   = rVbin15->Parameter(1);
        par0   = rVbin15->Parameter(0);
        par1err   = rVbin15->ParError(1);
        par0err   = rVbin15->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin15: " << par0 << " +/- " << par0err << " p1Vbin15: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(4);
        Value_TrueUnfbin16->SetTitle(asymlabel);
        Value_TrueUnfbin16->SetMarkerStyle(23);
        Value_TrueUnfbin16->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin16->SetMarkerSize(0.6);
        Value_TrueUnfbin16->GetXaxis()->SetTitle(asymlabel + " bin16 (true)");
        Value_TrueUnfbin16->GetYaxis()->SetTitle(asymlabel + " bin16 (unfolded)");
        Value_TrueUnfbin16->Draw("AP same");
        //Value_TrueUnfbin16->Fit("pol1");
        TFitResultPtr rVbin16 = Value_TrueUnfbin16->Fit("pol1", "S");
        par1   = rVbin16->Parameter(1);
        par0   = rVbin16->Parameter(0);
        par1err   = rVbin16->ParError(1);
        par0err   = rVbin16->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin16: " << par0 << " +/- " << par0err << " p1Vbin16: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(5);
        Value_TrueUnfbin17->SetTitle(asymlabel);
        Value_TrueUnfbin17->SetMarkerStyle(23);
        Value_TrueUnfbin17->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin17->SetMarkerSize(0.6);
        Value_TrueUnfbin17->GetXaxis()->SetTitle(asymlabel + " bin17 (true)");
        Value_TrueUnfbin17->GetYaxis()->SetTitle(asymlabel + " bin17 (unfolded)");
        Value_TrueUnfbin17->Draw("AP same");
        //Value_TrueUnfbin17->Fit("pol1");
        TFitResultPtr rVbin17 = Value_TrueUnfbin17->Fit("pol1", "S");
        par1   = rVbin17->Parameter(1);
        par0   = rVbin17->Parameter(0);
        par1err   = rVbin17->ParError(1);
        par0err   = rVbin17->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin17: " << par0 << " +/- " << par0err << " p1Vbin17: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->cd(6);
        Value_TrueUnfbin18->SetTitle(asymlabel);
        Value_TrueUnfbin18->SetMarkerStyle(23);
        Value_TrueUnfbin18->SetMarkerColor(kGreen - 1);
        Value_TrueUnfbin18->SetMarkerSize(0.6);
        Value_TrueUnfbin18->GetXaxis()->SetTitle(asymlabel + " bin18 (true)");
        Value_TrueUnfbin18->GetYaxis()->SetTitle(asymlabel + " bin18 (unfolded)");
        Value_TrueUnfbin18->Draw("AP same");
        //Value_TrueUnfbin18->Fit("pol1");
        TFitResultPtr rVbin18 = Value_TrueUnfbin18->Fit("pol1", "S");
        par1   = rVbin18->Parameter(1);
        par0   = rVbin18->Parameter(0);
        par1err   = rVbin18->ParError(1);
        par0err   = rVbin18->ParError(0);
        third_output_file << acceptanceName << " Linearity, f " << Nfunction << " p0Vbin18: " << par0 << " +/- " << par0err << " p1Vbin18: " << par1 << " +/- " << par1err <<  endl;

        c_LinearityBinByBin->SaveAs("2D_" + acceptanceName + "_LinearityCheck_binbybin.pdf");
        // c_LinearityBinByBin->SaveAs(acceptanceName + "_LinearityCheck_binbybin.C");













        TCanvas *c_Pull_lin = new TCanvas("c_Pull_lin", "c_Pull_lin", 500, 500);
        if (!plot_inclusive_only) c_Pull_lin->Divide(2, 2);
        if (!plot_inclusive_only) c_Pull_lin->cd(1);
        Asym2D_Pull->SetTitle(asymlabel);
        Asym2D_Pull->SetMarkerStyle(23);
        Asym2D_Pull->SetMarkerColor(kBlack);
        Asym2D_Pull->SetMarkerSize(0.6);
        Asym2D_Pull->GetXaxis()->SetTitle(asymlabel + " inclusive (true)");
        Asym2D_Pull->GetYaxis()->SetTitle(asymlabel + " inclusive pull");
        Asym2D_Pull->Draw("AP same");
        //Asym2D_Pull->Fit("pol1");
        TFitResultPtr rp = Asym2D_Pull->Fit("pol1", "S");
        par1   = rp->Parameter(1);
        par0   = rp->Parameter(0);
        par1err   = rp->ParError(1);
        par0err   = rp->ParError(0);
        third_output_file << acceptanceName << " Pull, f " << Nfunction << " p0: " << par0 << " +/- " << par0err << " p1: " << par1 << " +/- " << par1err <<  endl;

        if (!plot_inclusive_only)
        {
            c_Pull_lin->cd(2);
            Asym2D_Pullbin1->SetTitle(asymlabel);
            Asym2D_Pullbin1->SetMarkerStyle(23);
            Asym2D_Pullbin1->SetMarkerColor(kBlue);
            Asym2D_Pullbin1->SetMarkerSize(0.6);
            Asym2D_Pullbin1->GetXaxis()->SetTitle(yaxislabel + " bin 1 (true)");
            Asym2D_Pullbin1->GetYaxis()->SetTitle(yaxislabel + " bin 1 pull");
            Asym2D_Pullbin1->Draw("AP same");
            //Asym2D_Pullbin1->Fit("pol1");
            TFitResultPtr rpbin1 = Asym2D_Pullbin1->Fit("pol1", "S");
            par1   = rpbin1->Parameter(1);
            par0   = rpbin1->Parameter(0);
            par1err   = rpbin1->ParError(1);
            par0err   = rpbin1->ParError(0);
            third_output_file << acceptanceName << " Pull, f " << Nfunction << " p0bin1: " << par0 << " +/- " << par0err << " p1bin1: " << par1 << " +/- " << par1err <<  endl;

            c_Pull_lin->cd(3);
            Asym2D_Pullbin2->SetTitle(asymlabel);
            Asym2D_Pullbin2->SetMarkerStyle(23);
            Asym2D_Pullbin2->SetMarkerColor(kBlue);
            Asym2D_Pullbin2->SetMarkerSize(0.6);
            Asym2D_Pullbin2->GetXaxis()->SetTitle(yaxislabel + " bin 2 (true)");
            Asym2D_Pullbin2->GetYaxis()->SetTitle(yaxislabel + " bin 2 pull");
            Asym2D_Pullbin2->Draw("AP same");
            //Asym2D_Pullbin2->Fit("pol1");
            TFitResultPtr rpbin2 = Asym2D_Pullbin2->Fit("pol1", "S");
            par1   = rpbin2->Parameter(1);
            par0   = rpbin2->Parameter(0);
            par1err   = rpbin2->ParError(1);
            par0err   = rpbin2->ParError(0);
            third_output_file << acceptanceName << " Pull, f " << Nfunction << " p0bin2: " << par0 << " +/- " << par0err << " p1bin2: " << par1 << " +/- " << par1err <<  endl;

            c_Pull_lin->cd(4);
            Asym2D_Pullbin3->SetTitle(asymlabel);
            Asym2D_Pullbin3->SetMarkerStyle(23);
            Asym2D_Pullbin3->SetMarkerColor(kBlue);
            Asym2D_Pullbin3->SetMarkerSize(0.6);
            Asym2D_Pullbin3->GetXaxis()->SetTitle(yaxislabel + " bin 3 (true)");
            Asym2D_Pullbin3->GetYaxis()->SetTitle(yaxislabel + " bin 3 pull");
            Asym2D_Pullbin3->Draw("AP same");
            //Asym2D_Pullbin3->Fit("pol1");
            TFitResultPtr rpbin3 = Asym2D_Pullbin3->Fit("pol1", "S");
            par1   = rpbin3->Parameter(1);
            par0   = rpbin3->Parameter(0);
            par1err   = rpbin3->ParError(1);
            par0err   = rpbin3->ParError(0);
            third_output_file << acceptanceName << " Pull, f " << Nfunction << " p0bin3: " << par0 << " +/- " << par0err << " p1bin3: " << par1 << " +/- " << par1err <<  endl;
        }

        c_Pull_lin->SaveAs("2D_" + acceptanceName + "_LinearityCheck_Pull.pdf");
        // c_Pull_lin->SaveAs(acceptanceName + "_LinearityCheck_Pull.C");




        TCanvas *c_PullWidth_lin = new TCanvas("c_PullWidth_lin", "c_PullWidth_lin", 500, 500);
        if (!plot_inclusive_only) c_PullWidth_lin->Divide(2, 2);
        if (!plot_inclusive_only) c_PullWidth_lin->cd(1);
        Asym2D_PullWidth->SetTitle(asymlabel);
        Asym2D_PullWidth->SetMarkerStyle(23);
        Asym2D_PullWidth->SetMarkerColor(kBlack);
        Asym2D_PullWidth->SetMarkerSize(0.6);
        Asym2D_PullWidth->GetXaxis()->SetTitle(Var2D + " inclusive (true)");
        Asym2D_PullWidth->GetYaxis()->SetTitle(Var2D + " inclusive pull width");
        Asym2D_PullWidth->Draw("AP same");
        //Asym2D_PullWidth->Fit("pol1");
        TFitResultPtr rpw = Asym2D_PullWidth->Fit("pol1", "S");
        par1   = rpw->Parameter(1);
        par0   = rpw->Parameter(0);
        par1err   = rpw->ParError(1);
        par0err   = rpw->ParError(0);
        third_output_file << acceptanceName << " PullWidth, f " << Nfunction << " p0: " << par0 << " +/- " << par0err << " p1: " << par1 << " +/- " << par1err <<  endl;

        if (!plot_inclusive_only)
        {
            c_PullWidth_lin->cd(2);
            Asym2D_PullWidthbin1->SetTitle(asymlabel);
            Asym2D_PullWidthbin1->SetMarkerStyle(23);
            Asym2D_PullWidthbin1->SetMarkerColor(kBlue);
            Asym2D_PullWidthbin1->SetMarkerSize(0.6);
            Asym2D_PullWidthbin1->GetXaxis()->SetTitle(yaxislabel + " bin 1 (true)");
            Asym2D_PullWidthbin1->GetYaxis()->SetTitle(yaxislabel + " bin 1 pull width");
            Asym2D_PullWidthbin1->Draw("AP same");
            //Asym2D_PullWidthbin1->Fit("pol1");
            TFitResultPtr rpwbin1 = Asym2D_PullWidthbin1->Fit("pol1", "S");
            par1   = rpwbin1->Parameter(1);
            par0   = rpwbin1->Parameter(0);
            par1err   = rpwbin1->ParError(1);
            par0err   = rpwbin1->ParError(0);
            third_output_file << acceptanceName << " PullWidth, f " << Nfunction << " p0bin1: " << par0 << " +/- " << par0err << " p1bin1: " << par1 << " +/- " << par1err <<  endl;

            c_PullWidth_lin->cd(3);
            Asym2D_PullWidthbin2->SetTitle(asymlabel);
            Asym2D_PullWidthbin2->SetMarkerStyle(23);
            Asym2D_PullWidthbin2->SetMarkerColor(kBlue);
            Asym2D_PullWidthbin2->SetMarkerSize(0.6);
            Asym2D_PullWidthbin2->GetXaxis()->SetTitle(yaxislabel + " bin 2 (true)");
            Asym2D_PullWidthbin2->GetYaxis()->SetTitle(yaxislabel + " bin 2 pull width");
            Asym2D_PullWidthbin2->Draw("AP same");
            //Asym2D_PullWidthbin2->Fit("pol1");
            TFitResultPtr rpwbin2 = Asym2D_PullWidthbin2->Fit("pol1", "S");
            par1   = rpwbin2->Parameter(1);
            par0   = rpwbin2->Parameter(0);
            par1err   = rpwbin2->ParError(1);
            par0err   = rpwbin2->ParError(0);
            third_output_file << acceptanceName << " PullWidth, f " << Nfunction << " p0bin2: " << par0 << " +/- " << par0err << " p1bin2: " << par1 << " +/- " << par1err <<  endl;

            c_PullWidth_lin->cd(4);
            Asym2D_PullWidthbin3->SetTitle(asymlabel);
            Asym2D_PullWidthbin3->SetMarkerStyle(23);
            Asym2D_PullWidthbin3->SetMarkerColor(kBlue);
            Asym2D_PullWidthbin3->SetMarkerSize(0.6);
            Asym2D_PullWidthbin3->GetXaxis()->SetTitle(yaxislabel + " bin 3 (true)");
            Asym2D_PullWidthbin3->GetYaxis()->SetTitle(yaxislabel + " bin 3 pull width");
            Asym2D_PullWidthbin3->Draw("AP same");
            //Asym2D_PullWidthbin3->Fit("pol1");
            TFitResultPtr rpwbin3 = Asym2D_PullWidthbin3->Fit("pol1", "S");
            par1   = rpwbin3->Parameter(1);
            par0   = rpwbin3->Parameter(0);
            par1err   = rpwbin3->ParError(1);
            par0err   = rpwbin3->ParError(0);
            third_output_file << acceptanceName << " PullWidth, f " << Nfunction << " p0bin3: " << par0 << " +/- " << par0err << " p1bin3: " << par1 << " +/- " << par1err <<  endl;
        }

        c_PullWidth_lin->SaveAs("2D_" + acceptanceName + "_LinearityCheck_PullWidth.pdf");
        // c_PullWidth_lin->SaveAs(acceptanceName + "_LinearityCheck_PullWidth.C");

		//These plots don't make sense in 2D! Maybe someday I'll find an equivalent, but for now, they get commented out.
		/*
		TH2D* hMeas_before_combined = (TH2D*)( hMeas_after_array[3]->Clone("meas_before_combined") ); //A proxy for the non-weighted

        gStyle->SetOptStat(0);
        TCanvas *c_asymdist_lin = new TCanvas("c_asymdist_lin", "c_asymdist_lin", 1000, 500);
        c_asymdist_lin->Divide(4, 2);
        c_asymdist_lin->cd(1);
        hTrue_before->SetLineColor(TColor::GetColorDark(kRed));
        hTrue_before->SetLineWidth(1);
        hTrue_before->SetMinimum(0);
        hTrue_before->SetMaximum(1.3 * hTrue_before->GetMaximum());
        hTrue_before->SetFillStyle(0);
        hTrue_before->GetXaxis()->SetTitle(xaxislabel);
        hTrue_before->GetYaxis()->SetTitle("Number of events");
        hTrue_before->Draw("hist");
        hMeas_before_combined->SetLineColor(TColor::GetColorDark(kBlue));
        hMeas_before_combined->SetLineWidth(1);
        hMeas_before_combined->SetFillStyle(0);
        hMeas_before_combined->GetXaxis()->SetTitle(xaxislabel);
        hMeas_before_combined->GetYaxis()->SetTitle("Number of events");
        hMeas_before_combined->Draw("hist same");

        TLegend *leg1 = new TLegend(0.70, 0.76, 0.9, 0.93, NULL, "brNDC");
        leg1->SetEntrySeparation(0.1);
        leg1->SetFillColor(0);
        leg1->SetLineColor(0);
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->SetTextSize(0.03);
        leg1->AddEntry(hTrue_before,    "gen",  "L");
        leg1->AddEntry(hMeas_before_combined,    "reco", "L");
        leg1->Draw();



        for (int k = 0; k < Nlin; ++k)
        {
            slope = -0.3 + 0.1 * k;
            fx_scaled->SetParameters(slope, fscale);
            if (fabs(slope) < 0.001) slope = 0;
            // TString slope_temp = formatFloat(slope, "%6.1f");
			sprintf( tmpchar, "%6.1f", slope );
			TString slope_temp(tmpchar);
            slope_temp.ReplaceAll(" " , "" );
            c_asymdist_lin->cd(k + 2);
            hTrue_after_array[k]->SetLineColor(TColor::GetColorDark(kRed));
            hTrue_after_array[k]->SetLineWidth(1);
            hTrue_after_array[k]->SetMinimum(0);
            hTrue_after_array[k]->SetMaximum(1.3 * hTrue_after_array[k]->GetMaximum());
            hTrue_after_array[k]->SetFillStyle(0);
            hTrue_after_array[k]->GetXaxis()->SetTitle(xaxislabel + ", slope = " + slope_temp);
            hTrue_after_array[k]->GetYaxis()->SetTitle("Number of events");
            hTrue_after_array[k]->Draw("hist");
            hMeas_after_array[k]->SetLineColor(TColor::GetColorDark(kBlue));
            hMeas_after_array[k]->SetLineWidth(1);
            hMeas_after_array[k]->SetFillStyle(0);
            hMeas_after_array[k]->GetXaxis()->SetTitle(xaxislabel + ", slope = " + slope_temp);
            hMeas_after_array[k]->GetYaxis()->SetTitle("Number of events");
            hMeas_after_array[k]->Draw("hist same");
            fx_scaled->SetLineColor(TColor::GetColorDark(kGreen));
            fx_scaled->DrawCopy("LSAME");

            TPaveText *pt1 = new TPaveText(0.20, 0.80, 0.45, 0.93, "brNDC");
            pt1->SetName("pt1name");
            pt1->SetBorderSize(0);
            pt1->SetFillStyle(0);

            TText *blah;

            Float_t Afb, AfbErr;
            GetAfb(hTrue_after_array[k], Afb, AfbErr);

            // TString Asym1_temp = formatFloat(Afb, "%6.2f");
			sprintf( tmpchar, "%6.2f", Afb );
			TString Asym1_temp(tmpchar);
            Asym1_temp.ReplaceAll(" " , "" );
            Asym1_temp = TString(" Asym (gen): ") +  Asym1_temp;
            blah = pt1->AddText(Asym1_temp.Data());
            blah->SetTextSize(0.03);
            blah->SetTextAlign(11);
            blah->SetTextColor(TColor::GetColorDark(kRed));

            GetAfb(hMeas_after_array[k], Afb, AfbErr);

            // TString Asym2_temp = formatFloat(Afb, "%6.2f");
			sprintf( tmpchar, "%6.2f", Afb );
			TString Asym2_temp(tmpchar);
            Asym2_temp.ReplaceAll(" " , "" );
            Asym2_temp = TString(" Asym (reco): ") +  Asym2_temp;
            blah = pt1->AddText(Asym2_temp.Data());
            blah->SetTextSize(0.03);
            blah->SetTextAlign(11);
            blah->SetTextColor(TColor::GetColorDark(kBlue));

            pt1->Draw();

            TLegend *leg2 = new TLegend(0.70, 0.84, 0.9, 0.93, NULL, "brNDC");
            leg2->SetEntrySeparation(0.1);
            leg2->SetFillColor(0);
            leg2->SetLineColor(0);
            leg2->SetBorderSize(0);
            leg2->SetFillStyle(0);
            leg2->SetTextSize(0.03);
            leg2->AddEntry(fx_scaled,    "weight function",  "L");
            leg2->Draw();

        }

        c_asymdist_lin->SaveAs("2D_" + acceptanceName + "_LinearityCheck_AsymDists.pdf");
        // c_asymdist_lin->SaveAs(acceptanceName + "_LinearityCheck_AsymDists.C");
		*/

    } // end "if TestType==Linearity"
    else
    { // presumably TestType==Pull


        TCanvas *c_pull = new TCanvas("c_pull", "c_pull", 500, 500);
        if (!plot_inclusive_only) c_pull->Divide(2, 2);
        if (!plot_inclusive_only) c_pull->cd(1);
        AfbPull[0]->SetMarkerStyle(23);
        AfbPull[0]->SetMarkerColor(kBlack);
        AfbPull[0]->SetMarkerSize(0.6);
        AfbPull[0]->GetXaxis()->SetTitle(asymlabel + " inclusive pull");
        AfbPull[0]->GetYaxis()->SetTitle("Number of PEs / 0.2");
        AfbPull[0]->Fit("gaus");
        AfbPull[0]->Draw();

        if (!plot_inclusive_only)
        {
            c_pull->cd(2);
            AfbPull[1]->SetMarkerStyle(23);
            AfbPull[1]->SetMarkerColor(kBlue);
            AfbPull[1]->SetMarkerSize(0.6);
            AfbPull[1]->GetXaxis()->SetTitle(yaxislabel + " bin 1 pull");
            AfbPull[1]->GetYaxis()->SetTitle("Number of PEs / 0.2");
            AfbPull[1]->Fit("gaus");
            AfbPull[1]->Draw();

            c_pull->cd(3);
            AfbPull[2]->SetMarkerStyle(23);
            AfbPull[2]->SetMarkerColor(kBlue);
            AfbPull[2]->SetMarkerSize(0.6);
            AfbPull[2]->GetXaxis()->SetTitle(yaxislabel + " bin 2 pull");
            AfbPull[2]->GetYaxis()->SetTitle("Number of PEs / 0.2");
            AfbPull[2]->Fit("gaus");
            AfbPull[2]->Draw();

            c_pull->cd(4);
            AfbPull[3]->SetMarkerStyle(23);
            AfbPull[3]->SetMarkerColor(kBlue);
            AfbPull[3]->SetMarkerSize(0.6);
            AfbPull[3]->GetXaxis()->SetTitle(yaxislabel + " bin 3 pull");
            AfbPull[3]->GetYaxis()->SetTitle("Number of PEs / 0.2");
            AfbPull[3]->Fit("gaus");
            AfbPull[3]->Draw();
        }

        c_pull->SaveAs("2D_" + acceptanceName + "_Pull.pdf");
        // c_pull->SaveAs(acceptanceName + "_Pull.C");


		
        // TFile *plots = new TFile(acceptanceName + "_plots.root", "RECREATE");
		/*
        for (int i = 0; i < nbins2D; i++)
        {
            h_pulls[i] ->Write();
            h_resd[i] ->Write();
        }
		*/
        AfbPull[0]->Write();
        AfbPull[1]->Write();
        AfbPull[2]->Write();
        AfbPull[3]->Write();

    }

    //myfile.close();
    //second_output_file.close();
    third_output_file.close();

	//delete[] AfbPull;
}

#ifndef __CINT__
int main ()
{
    AfbUnfoldTests();    // Main program when run stand-alone
    return 0;
}
#endif
