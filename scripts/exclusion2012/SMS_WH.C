//#include "Utils/SMS_utils.C"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>
#include "TCanvas.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>

#include "limits_TChiWH_MET100.C"
#include "limits_TChiWH_MET125.C"
#include "limits_TChiWH_MET150.C"
#include "limits_TChiWH_MET175.C"

using namespace std;

void SMS_WH(char* sample = "TChiWH" , bool print = true){

  gStyle->SetPaintTextFormat(".2f");

  //--------------------------------------------------
  // input parameters
  //--------------------------------------------------

  const float lumi      = 19600; 
  char* suffix          = (char*) "";
  char* denomhistoname  = (char*) "masses";
  char* filename        = (char*) "";
  char* denomname       = (char*) "";
  char* label           = (char*)"";
  bool  doISRWeight     = true;

  bool doBDT            = false;
  char* BDTchar         = (char*)"";

  cout << "Doing cut-based analysis" << endl;

  //--------------------------------------------------
  // set up input files
  //--------------------------------------------------

  if( TString(sample).Contains("TChiWH") ){    
    filename  = (char*) "/tas/olivito/data/whmet/output_V00-02-36_2012/Minibabies_V00-04-05/TChiWH_whmet.root";
    //    filename  = (char*) "/tas01/disk03/olivito/whmet/output_V00-02-36_2012/Minibabies_V00-03-16/TChiWH_whmet.root";
    // filename  = (char*) "/tas01/disk03/olivito/whmet/output_V00-02-36_2012_2jskim/Minibabies_V00-04-01_BDT_V00-00-01/TChiWH_whmet.root";
    //denomname = (char*) "/tas/cms2/stop/cms2V05-03-26_stoplooperV00-02-24/T2tt_mad/minibaby_V00-03-06/Skim_4jets_MET100_MT120/myMassDB_T2tt_combined_25GeVbins.root";
    denomname = (char*) "/tas/olivito/data/whmet/output_V00-02-36_2012/Minibabies_V00-04-05/myMassDB_TChiWH.root";

    label     = (char*)"pp #rightarrow #tilde{#chi}_{1}^{#pm}#tilde{#chi}_{2}^{0} #rightarrow W(l#nu)H(b#bar{b})#tilde{#chi}_{1}^{0}#tilde{#chi}_{1}^{0}";
  }

  else{
    cout << "ERROR! unrecognized sample " << sample << ", quitting!!!" << endl;
    exit(0);
  }

  //--------------------------------------------------
  // set up logfile
  //--------------------------------------------------

  char* logfilename = Form("logfiles/%s%s.log",sample,suffix);
  ofstream* logfile = new ofstream();
  logfile->open(logfilename,ios::trunc);

  cout << "----------------------------------------------------" << endl;
  cout << "Sample            " << sample          << endl;
  cout << "logfile           " << logfilename     << endl;
  cout << "Using file        " << filename        << endl;
  cout << "Using denominator " << denomname       << endl;
  cout << "Denom histo       " << denomhistoname  << endl;
  cout << "Using lumi        " << lumi            << " pb-1" << endl;
  cout << "----------------------------------------------------" << endl;

  *logfile << "----------------------------------------------------" << endl;
  *logfile << "Sample            " << sample          << endl;
  *logfile << "logfile           " << logfilename     << endl;
  *logfile << "Using file        " << filename        << endl;
  *logfile << "Using denominator " << denomname       << endl;
  *logfile << "Denom histo       " << denomhistoname  << endl;
  *logfile << "Using lumi        " << lumi            << " pb-1" << endl;
  *logfile << "----------------------------------------------------" << endl;

  TFile* fdenom = TFile::Open(denomname);
  TH2F*  hdenom = (TH2F*) fdenom->Get(denomhistoname);

  //--------------------------------------------------
  // read in TChain
  //--------------------------------------------------

  TChain *ch = new TChain("t");
  ch->Add(filename);

  //--------------------------------------------------
  // read in reference cross section
  //--------------------------------------------------

  TFile *xsecfile = TFile::Open("c1n2_xsec.root");
  TH1F* refxsec   = (TH1F*) xsecfile->Get("h_c1n2_xsec");

  //--------------------------------------------------
  // preselection
  //--------------------------------------------------

  TCut sig("mini_whsig == 1");
  TCut njets("mini_njets == 2 && mini_njets_fwd == 0");
  TCut mt100("mini_mt > 100.");
  TCut mt2bl("mini_mt2bl>200.0");
  TCut met100("mini_met > 100.0");
  TCut met125("mini_met > 125.0");
  TCut met150("mini_met > 150.0");
  TCut met175("mini_met > 175.0");
  TCut testing("event%2==0");

  // weights
  TCut isrweight("mini_isrweight");
  TCut sltrigeff("mini_sltrigeff");
  TCut btagsf("mini_btagsf");
  TCut lepfastsimsf("mini_lepfastsimsf");

  TCut presel;

  //-------------------------------------------
  // THESE CUTS DEFINE PRESELECTION REGION
  //-------------------------------------------

  presel += sig;
  presel += njets;
  presel += mt100;
  presel += mt2bl;

  // if( doBDT ){
  //   presel += testing;
  // }
  if( !doISRWeight ){
    isrweight = TCut("1");
  }

  TCut BDTweight("1");
  //  if( doBDT ) BDTweight = TCut("2");

  TCut weight = sltrigeff * btagsf * BDTweight * lepfastsimsf;

  TCut whweight = "mini_whweight";

  cout << "Using pre-selection   " << presel.GetTitle()    << endl;
  cout << "Using weight          " << weight.GetTitle()    << endl;
  cout << "Using ISR weight      " << isrweight.GetTitle() << endl;

  *logfile << "Using pre-selection   " << presel.GetTitle()    << endl;
  *logfile << "Using weight          " << weight.GetTitle()    << endl;
  *logfile << "Using ISR weight      " << isrweight.GetTitle() << endl;

  //--------------------------------------------------
  // signal regions
  //--------------------------------------------------

  vector<TCut>    sigcuts;
  vector<string>  signames;
  vector<string>  labels;
  vector<float>   uls;

  if( doBDT ){

    // //-----------------------------
    // // T2tt BDT
    // //-----------------------------

    // if( TString(sample).Contains("T2tt") ){

    //   cout << "Doing T2tt BDT signal regions" << endl;  

    //   sigcuts.push_back(TCut(presel+"mini_bdt[1] > 0.30"));  signames.push_back("T2TT_BDT1L");  labels.push_back("T2TT_BDT1L");  uls.push_back(-1.0);
    //   sigcuts.push_back(TCut(presel+"mini_bdt[1] > 0.40"));  signames.push_back("T2TT_BDT1T");  labels.push_back("T2TT_BDT1T");  uls.push_back(-1.0);
    //   sigcuts.push_back(TCut(presel+"mini_bdt[2] > 0.55"));  signames.push_back("T2TT_BDT2");   labels.push_back("T2TT_BDT2");   uls.push_back(-1.0);
    //   sigcuts.push_back(TCut(presel+"mini_bdt[3] > 0.65"));  signames.push_back("T2TT_BDT3");   labels.push_back("T2TT_BDT3");   uls.push_back(-1.0);
    //   sigcuts.push_back(TCut(presel+"mini_bdt[4] > 0.50"));  signames.push_back("T2TT_BDT4");   labels.push_back("T2TT_BDT4");   uls.push_back(-1.0);
    //   sigcuts.push_back(TCut(presel+"mini_bdt[5] > 0.30"));  signames.push_back("T2TT_BDT5");   labels.push_back("T2TT_BDT5");   uls.push_back(-1.0);
    // }

  }

  else{

    //-----------------------------
    // TChiWH cut-and-count
    //-----------------------------

    if( TString(sample).Contains("TChiWH") ){

      cout << "Doing TChiWH cut-based signal regions" << endl;  

      // signal regions defined by met cuts
      // XX fixme: where do uls come from??? is it num events? - must be
      sigcuts.push_back(presel+met100);        signames.push_back("TChiWH_MET100");   labels.push_back("TChiWH_MET100");   uls.push_back(79.6);
      sigcuts.push_back(presel+met125);        signames.push_back("TChiWH_MET125");   labels.push_back("TChiWH_MET125");   uls.push_back(42.1);
      sigcuts.push_back(presel+met150);        signames.push_back("TChiWH_MET150");   labels.push_back("TChiWH_MET150");   uls.push_back(19.0);
      sigcuts.push_back(presel+met175);        signames.push_back("TChiWH_MET175");   labels.push_back("TChiWH_MET175");   uls.push_back( 9.9);

    }

  }

  const unsigned int nsig = sigcuts.size();

  //--------------------------------------------------
  // make efficiency and xsec TH2's
  //--------------------------------------------------
  
  TH2F* heff[nsig];
  TH2F* hnsig[nsig];
  TH2F* hnevents[nsig];
  TH2F* heff_noisr[nsig];
  TH2F* heffup[nsig];
  TH2F* heffdn[nsig];
  TH2F* heffbup[nsig];
  TH2F* heffbdn[nsig];
  TH2F* hxsec[nsig];
  TH2F* hxsec_exp[nsig];
  TH2F* hxsec_expp1[nsig];
  TH2F* hxsec_expm1[nsig];
  TH2F* hexcl[nsig];
  TH2F* hexcl_exp[nsig];
  TH2F* hexcl_expp1[nsig];
  TH2F* hexcl_expm1[nsig];
  TH2F* hexcl_obsp1[nsig];
  TH2F* hexcl_obsm1[nsig];
  TH2F* hjes[nsig];
  TH2F* hbtagerr[nsig];
  TH2F* hstaterr[nsig];
  TH2F* htoterr[nsig];
  TH2F* hisrerr[nsig];
  
  TCanvas *ctemp = new TCanvas();
  ctemp->cd();

  for( unsigned int i = 0 ; i < nsig ; ++i ){

    //------------------------------------------------
    // calculate point-by-point signal uncertainties
    //------------------------------------------------

    TString jesup(sigcuts.at(i));
    jesup.ReplaceAll("mini_njets"         , "mini_njetsup"    );
    jesup.ReplaceAll("mini_njets_fwd"     , "mini_njets_fwdup"    );
    jesup.ReplaceAll("mini_met"           , "mini_metup"      );
    jesup.ReplaceAll(" mini_mt"           , "mini_mtup"       );
    jesup.ReplaceAll("mini_mt2bl"         , "mini_mt2blup"     );
    //    jesup.ReplaceAll("mini_bdt"           , "mini_bdtup"      );

    TString jesdown(sigcuts.at(i));
    jesdown.ReplaceAll("mini_njets"       , "mini_njetsdown"  );
    jesup.ReplaceAll("mini_njets_fwd"     , "mini_njets_fwddown"    );
    jesdown.ReplaceAll("mini_met"         , "mini_metdown"    );
    jesdown.ReplaceAll(" mini_mt"         , "mini_mtdown"     );
    jesdown.ReplaceAll("mini_mt2bl"       , "mini_mt2bldown"   );
    //    jesdown.ReplaceAll("mini_bdt"         , "mini_bdtdown"    );

    // TString btagup(sigcuts.at(i));
    // btagup.ReplaceAll("mini_nb"           , "mini_nbupBC"     );
    // btagup.ReplaceAll("mini_chi2"         , "mini_chi2bup"    );
    // btagup.ReplaceAll("mini_mt2w"         , "mini_mt2wbup"    );
    // btagup.ReplaceAll("mini_pt_b"         , "mini_pt_b_bup"   );
    // btagup.ReplaceAll("mini_bdt"          , "mini_bdtbup"     );

    // TString btagdn(sigcuts.at(i));
    // btagdn.ReplaceAll("mini_nb"         , "mini_nbdownBC"   );
    // btagdn.ReplaceAll("mini_chi2"       , "mini_chi2bdown"  );
    // btagdn.ReplaceAll("mini_mt2w"       , "mini_mt2wbdown"  );
    // btagdn.ReplaceAll("mini_pt_b"       , "mini_pt_b_bdown" );
    // btagdn.ReplaceAll("mini_bdt"        , "mini_bdtbdown"   );

    TCut jesupcut(jesup);
    TCut jesdncut(jesdown);
    // TCut btagupcut(btagup);
    // TCut btagdncut(btagdn);

    cout << endl << endl << endl;
    cout << "Signal region       : " << labels.at(i)   << endl << endl;
    cout << "Selection           : " << sigcuts.at(i).GetTitle()  << endl << endl;
    cout << "Selection JES up    : " << jesupcut.GetTitle()       << endl << endl;
    cout << "Selection JES down  : " << jesdncut.GetTitle()       << endl << endl;
    // cout << "Selection btag up   : " << btagupcut      << endl << endl;
    // cout << "Selection btag down : " << btagdncut      << endl << endl;

    *logfile << "Signal region       : " << labels.at(i)   << endl << endl;
    *logfile << "Selection           : " << sigcuts.at(i).GetTitle()  << endl << endl;
    *logfile << "Selection JES up    : " << jesupcut.GetTitle()       << endl << endl;
    *logfile << "Selection JES down  : " << jesdncut.GetTitle()       << endl << endl;
    // *logfile << "Selection btag up   : " << btagupcut      << endl << endl;
    // *logfile << "Selection btag down : " << btagdncut      << endl << endl;

    int   nbinsx  =      21;
    float xmin    =   -12.5;
    float xmax    =   512.5;
    int   nbinsy  =      21;
    float ymin    =   -12.5;
    float ymax    =   512.5;

    heff[i]        = new TH2F(Form("heff_%i",i)          , Form("heff_%i",i)         , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hnsig[i]       = new TH2F(Form("hnsig_%i",i)         , Form("hnsig_%i",i)        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hnevents[i]    = new TH2F(Form("hnevents_%i",i)      , Form("hnevents_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heff_noisr[i]  = new TH2F(Form("heff_noisr_%i",i)    , Form("heff_noisr_%i",i)   , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffup[i]      = new TH2F(Form("heffup_%i",i)        , Form("heffup_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    heffdn[i]      = new TH2F(Form("heffdn_%i",i)        , Form("heffdn_%i",i)       , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    // heffbup[i]     = new TH2F(Form("heffbup_%i",i)       , Form("heffbup_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    // heffbdn[i]     = new TH2F(Form("heffbdn_%i",i)       , Form("heffbdn_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec[i]       = new TH2F(Form("hxsec_%i",i)         , Form("hxsec_%i",i)        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec_exp[i]   = new TH2F(Form("hxsec_exp_%i",i)     , Form("hxsec_exp_%i",i)    , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec_expp1[i] = new TH2F(Form("hxsec_expp1_%i",i)   , Form("hxsec_expp1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hxsec_expm1[i] = new TH2F(Form("hxsec_expm1_%i",i)   , Form("hxsec_expm1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
 
    hexcl[i]       = new TH2F(Form("hexcl_%i",i)         , Form("hexcl_%i",i)        , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_obsp1[i] = new TH2F(Form("hexcl_obsp1_%i",i)   , Form("hexcl_obsp1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_obsm1[i] = new TH2F(Form("hexcl_obsm1_%i",i)   , Form("hexcl_obsm1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_exp[i]   = new TH2F(Form("hexcl_exp_%i",i)     , Form("hexcl_exp_%i",i)    , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_expp1[i] = new TH2F(Form("hexcl_expp1_%i",i)   , Form("hexcl_expp1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hexcl_expm1[i] = new TH2F(Form("hexcl_expm1_%i",i)   , Form("hexcl_expm1_%i",i)  , nbinsx , xmin , xmax , nbinsy , ymin , ymax );

    hjes[i]        = new TH2F(Form("hjes_%i",i)          , Form("hjes_%i",i)         , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    htoterr[i]     = new TH2F(Form("htoterr_%i",i)       , Form("htoterr_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hisrerr[i]     = new TH2F(Form("hisrerr_%i",i)       , Form("hisrerr_%i",i)      , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    //    hbtagerr[i]    = new TH2F(Form("hbtagerr_%i",i)      , Form("hbtagerr_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );
    hstaterr[i]    = new TH2F(Form("hstaterr_%i",i)      , Form("hstaterr_%i",i)     , nbinsx , xmin , xmax , nbinsy , ymin , ymax );

    heff[i]->Sumw2();
    hnsig[i]->Sumw2();

    // nominal
    ch->Draw(Form("mini_mlsp:mini_mchargino>>heff_%i",i),sigcuts.at(i)*weight*isrweight);
    heff[i]->Divide(hdenom);
    ch->Draw(Form("mini_mlsp:mini_mchargino>>hnsig_%i",i),sigcuts.at(i)*weight*isrweight*whweight);

    // raw number of signal events
    ch->Draw(Form("mini_mlsp:mini_mchargino>>hnevents_%i",i),sigcuts.at(i));

    // // JES up
    // ch->Draw(Form("mini_mlsp:mini_mchargino>>heffup_%i",i),jesupcut*weight*isrweight);
    // heffup[i]->Divide(hdenom);

    // // JES down
    // ch->Draw(Form("mini_mlsp:mini_mchargino>>heffdn_%i",i),jesdncut*weight*isrweight);
    // heffdn[i]->Divide(hdenom);

    // // btag up
    // ch->Draw(Form("mini_mlsp:mini_mchargino>>heffbup_%i",i),btagupcut*weight*isrweight);
    // heffbup[i]->Divide(hdenom);

    // // btag down
    // ch->Draw(Form("mini_mlsp:mini_mchargino>>heffbdn_%i",i),btagdncut*weight*isrweight);
    // heffbdn[i]->Divide(hdenom);

    // remove ISR weight
    ch->Draw(Form("mini_mlsp:mini_mchargino>>heff_noisr_%i",i),sigcuts.at(i)*weight);
    heff_noisr[i]->Divide(hdenom);

    for( int ibin = 1 ; ibin <= nbinsx ; ibin++ ){
      for( int jbin = 1 ; jbin <= nbinsy ; jbin++ ){

	//	hjes[i]->SetBinContent(ibin,jbin,0.0);
	//	hbtagerr[i]->SetBinContent(ibin,jbin,0.0);
	htoterr[i]->SetBinContent(ibin,jbin,0.0);
	hisrerr[i]->SetBinContent(ibin,jbin,0.0);
	hstaterr[i]->SetBinContent(ibin,jbin,0.0);

	float mg = heff[i]->GetXaxis()->GetBinCenter(ibin);
	//float ml = heff[i]->GetYaxis()->GetBinCenter(jbin);

	// nominal
	float eff          = heff[i]->GetBinContent(ibin,jbin);
	float eff_staterr  = 1.0/sqrt(hnevents[i]->GetBinContent(ibin,jbin));

	// JES up/down
	// float effup      = heffup[i]->GetBinContent(ibin,jbin);
	// float effdn      = heffdn[i]->GetBinContent(ibin,jbin);

	// // btag up/down
	// float effbup     = heffbup[i]->GetBinContent(ibin,jbin);
	// float effbdn     = heffbdn[i]->GetBinContent(ibin,jbin);

	// no ISR weighting
	float eff_noisr  = heff_noisr[i]->GetBinContent(ibin,jbin);

	if( eff   < 1e-20 ) continue;

	// JES uncertainty
	// float dup    = fabs(effup/eff-1);
	// float ddn    = fabs(1-effdn/eff);
	// float djes   = 0.5 * (dup+ddn);

	// // btag uncertainty
	// float dbup    = fabs(effbup/eff-1);
	// float dbdn    = fabs(1-effbdn/eff);
	// float dbtag   = 0.5 * (dbup+dbdn);

	// ISR uncertainty
	float disr   = fabs(1-eff/eff_noisr);

	// lumi (4.4%), trigger (5%), lepton selection (5%), btagging (5%), JES, ISR, 
	//	float toterr  = sqrt( 0.044*0.044 + 0.05*0.05 + 0.05*0.05 + 0.05*0.05 + djes*djes + disr*disr + eff_staterr*eff_staterr );
	//	float toterr  = sqrt( 0.044*0.044 + 0.05*0.05 + 0.05*0.05 + 0.05*0.05 + djes*djes + disr*disr );
	// use 0.05 temporarily for JES
	float toterr  = sqrt( 0.044*0.044 + 0.05*0.05 + 0.05*0.05 + 0.05*0.05 + 0.05*0.05 + disr*disr );
	htoterr[i]->SetBinContent(ibin,jbin,toterr);
	hisrerr[i]->SetBinContent(ibin,jbin,disr);
	//	hbtagerr[i]->SetBinContent(ibin,jbin,dbtag);
	//	hjes[i]->SetBinContent(ibin,jbin,djes);
	hstaterr[i]->SetBinContent(ibin,jbin,eff_staterr);

	float this_ul;
	float this_ul_exp;
	float this_ul_expp1;
	float this_ul_expm1;

	// this_ul       = uls.at(i);
	// this_ul_exp   = uls.at(i);
	// this_ul_expp1 = uls.at(i);
	// this_ul_expm1 = uls.at(i);
	
	//------------------------------------------
	// TChiWH cut-and-count
	//------------------------------------------

	if( TString(labels.at(i)).Contains("TChiWH_MET100") ){

	  this_ul       = getUpperLimit_TChiWH_MET100( toterr );
	  this_ul_exp   = getExpectedUpperLimit_TChiWH_MET100( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_TChiWH_MET100( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_TChiWH_MET100( toterr );
	}

	else if( TString(labels.at(i)).Contains("TChiWH_MET125") ){

	  this_ul       = getUpperLimit_TChiWH_MET125( toterr );
	  this_ul_exp   = getExpectedUpperLimit_TChiWH_MET125( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_TChiWH_MET125( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_TChiWH_MET125( toterr );
	}

	else if( TString(labels.at(i)).Contains("TChiWH_MET150") ){

	  this_ul       = getUpperLimit_TChiWH_MET150( toterr );
	  this_ul_exp   = getExpectedUpperLimit_TChiWH_MET150( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_TChiWH_MET150( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_TChiWH_MET150( toterr );
	}

	else if( TString(labels.at(i)).Contains("TChiWH_MET175") ){

	  this_ul       = getUpperLimit_TChiWH_MET175( toterr );
	  this_ul_exp   = getExpectedUpperLimit_TChiWH_MET175( toterr );
	  this_ul_expp1 = getExpectedP1UpperLimit_TChiWH_MET175( toterr );
	  this_ul_expm1 = getExpectedM1UpperLimit_TChiWH_MET175( toterr );
	}

	else{
	  cout << "ERROR! UNRECOGNIZED SIGNAL REGION " << labels.at(i) << endl;
	  exit(0);
	}

	// cout << endl;
	// cout << "toterr " << toterr << endl;
	// cout << "this_ul        " << this_ul       << endl;
	// cout << "this_ul_exp    " << this_ul_exp   << endl;
	// cout << "this_ul_expp1  " << this_ul_expp1 << endl;
	// cout << "this_ul_expm1  " << this_ul_expm1 << endl;

	float xsecul        = this_ul       / ( lumi * eff );
	float xsecul_exp    = this_ul_exp   / ( lumi * eff );
	float xsecul_expp1  = this_ul_expp1 / ( lumi * eff );
	float xsecul_expm1  = this_ul_expm1 / ( lumi * eff );

	if( eff > 0 ){
	  hxsec[i]->SetBinContent(ibin,jbin, xsecul );
	  hxsec_exp[i]->SetBinContent(ibin,jbin, xsecul_exp );
	  hxsec_expp1[i]->SetBinContent(ibin,jbin, xsecul_expp1 );
	  hxsec_expm1[i]->SetBinContent(ibin,jbin, xsecul_expm1 );
	}

	int   bin     = refxsec->FindBin(mg);
	// need to scale by BR(W->lv) * BR(H->bb) = 0.33 * 0.56
	float xsec    = refxsec->GetBinContent(bin) * 0.33 * 0.56;
	float xsec_up = (refxsec->GetBinContent(bin) + refxsec->GetBinError(bin)) * 0.33 * 0.56;
	float xsec_dn = (refxsec->GetBinContent(bin) - refxsec->GetBinError(bin)) * 0.33 * 0.56;

	hexcl[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul )   hexcl[i]->SetBinContent(ibin,jbin,1);

	hexcl_exp[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul_exp )   hexcl_exp[i]->SetBinContent(ibin,jbin,1);

	hexcl_expp1[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul_expp1 )   hexcl_expp1[i]->SetBinContent(ibin,jbin,1);

	hexcl_expm1[i]->SetBinContent(ibin,jbin,0);
	if( xsec > xsecul_expm1 )   hexcl_expm1[i]->SetBinContent(ibin,jbin,1);

	hexcl_obsp1[i]->SetBinContent(ibin,jbin,0);
	if( xsec_up > xsecul )   hexcl_obsp1[i]->SetBinContent(ibin,jbin,1);

	hexcl_obsm1[i]->SetBinContent(ibin,jbin,0);
	if( xsec_dn > xsecul )   hexcl_obsm1[i]->SetBinContent(ibin,jbin,1);

	//cout << "ibin jbin mg xsec " << ibin << " " << jbin << " " << mg << " " << xsec << endl;
      }
    }
  }

  delete ctemp;

  cout << endl << endl;

  //--------------------------------------------------
  // make pretty pictures
  //--------------------------------------------------
  
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextSize(0.04);

  TCanvas* can[nsig];
  //TCanvas* can_exclusion[nsig];

  TGraph* gr[nsig];
  TGraph* gr_exp[nsig];
  //TGraph* gr_expp1[nsig];
  //TGraph* gr_expm1[nsig];


  for( unsigned int i = 0 ; i < nsig ; ++i ){

    //TGraph *gr     = getGraph( sample , "observed" , signames.at(i) );
    //TGraph *gr_exp = getGraph( sample , "expected" , signames.at(i) );

    // need to multiply by BR for W->lv (0.33) * H->bb (0.56)
    // gr[i]      = getRefXsecGraph(hxsec[i]     , (char*) "TChiWH", 0.33 * 0.56);
    // gr_exp[i]  = getRefXsecGraph(hxsec_exp[i] , (char*) "TChiWH", 0.33 * 0.56);

    // gr[i]->SetLineWidth(3);
    // gr_exp[i]->SetLineWidth(3);
    // gr_exp[i]->SetLineStyle(2);

    // gr[i]->SetName(Form("gr_%i",i));
    // gr_exp[i]->SetName(Form("gr_exp_%i",i));
    // gr[i]->SetTitle(Form("gr_%i",i));
    // gr_exp[i]->SetTitle(Form("gr_exp_%i",i));

    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    //can[i]->Divide(2,1);
    //can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1800,600);
    //can[i]->Divide(3,1);

    can[i] = new TCanvas(Form("can_%i",i),Form("can_%i",i),1200,600);
    can[i]->Divide(2,1);

    //-------------------------------
    // efficiency
    //-------------------------------

    can[i]->cd(1);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    heff[i]->Scale(100);
    heff[i]->GetXaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetLabelSize(0.035);
    heff[i]->GetYaxis()->SetTitle("M_{#tilde{#chi}^{0}_{1}} [GeV]");
    heff[i]->GetYaxis()->SetTitleOffset(1.15);
    heff[i]->GetXaxis()->SetTitle("M_{#tilde{#chi}^{#pm}_{1}} [GeV]");
    heff[i]->GetZaxis()->SetTitle("efficiency (%)");
    heff[i]->GetZaxis()->SetTitleOffset(1.2);
    heff[i]->GetXaxis()->SetRangeUser(100,800);
    heff[i]->GetYaxis()->SetRangeUser(0,600);
    heff[i]->Draw("colz");
    heff[i]->SetMinimum(0.0);
    heff[i]->SetMaximum(4.0);
    //heff[i]->Draw("sametext");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.78,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.6 fb^{-1}");

    //-------------------------------
    // cross section
    //-------------------------------
  
    can[i]->cd(2);
    gPad->SetTopMargin(0.1);
    gPad->SetRightMargin(0.2);
    gPad->SetLogz();

    hxsec[i]->GetXaxis()->SetLabelSize(0.035);
    hxsec[i]->GetYaxis()->SetLabelSize(0.035);
    hxsec[i]->GetYaxis()->SetTitle("M_{#tilde{#chi}^{0}_{1}} [GeV]");
    hxsec[i]->GetYaxis()->SetTitleOffset(1.15);
    hxsec[i]->GetXaxis()->SetTitle("M_{#tilde{#chi}^{#pm}_{1}} [GeV]");
    hxsec[i]->GetZaxis()->SetTitle("#sigma upper limit [pb]");
    hxsec[i]->GetZaxis()->SetTitleOffset(1.2);
    hxsec[i]->Draw("colz");
    //hxsec[i]->Draw("sametext");
    hxsec[i]->SetMinimum(0.005);
    hxsec[i]->SetMaximum(100);
    hxsec[i]->GetXaxis()->SetRangeUser(100,800);
    hxsec[i]->GetYaxis()->SetRangeUser(0,600);
    hexcl[i]->Draw("samebox");

    // gr[i]->Draw("same");
    // gr_exp[i]->Draw("same");

    // TLegend *leg = new TLegend(0.2,0.6,0.4,0.75);
    // leg->AddEntry(gr[i],    "observed" ,"l");
    // leg->AddEntry(gr_exp[i],"expected" ,"l");
    // leg->SetFillColor(0);
    // leg->SetBorderSize(0);
    // leg->Draw();


    // TGraph* gr_excl      = getRefXsecGraph(hxsec[i], "T5zz", 1.0);
    // TGraph* gr_excl_down = getRefXsecGraph(hxsec[i], "T5zz", 1./3.);
    // TGraph* gr_excl_up   = getRefXsecGraph(hxsec[i], "T5zz", 3.);

    // gr_excl->SetLineWidth(2);
    // gr_excl_up->SetLineWidth(2);
    // gr_excl_down->SetLineWidth(2);
    // gr_excl_up->SetLineStyle(2);
    // gr_excl_down->SetLineStyle(3);
    // gr_excl->Draw("same");
    // gr_excl_up->Draw("same");
    // gr_excl_down->Draw("same");


    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.78,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.6 fb^{-1}");

    //-------------------------------
    // excluded points
    //-------------------------------
    /*
    can_exclusion[i] = new TCanvas(Form("can_exclusion_%i",i),Form("can_exclusion_%i",i),1200,600);
    can_exclusion[i]->Divide(2,1);

    can_exclusion[i]->cd(1);    
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hexcl[i]->SetMinimum(0);
    hexcl[i]->SetMaximum(1);
    hexcl[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hexcl[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hexcl[i]->GetZaxis()->SetTitle("observed excluded points");
    hexcl[i]->Draw("colz");
    gr[i]->Draw("l");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary  #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.6 fb^{-1}");

    //-------------------------------
    // JES uncertainty
    //-------------------------------

    can_exclusion[i]->cd(2);
    gPad->SetRightMargin(0.2);
    gPad->SetTopMargin(0.1);
    hexcl_exp[i]->SetMinimum(0);
    hexcl_exp[i]->SetMaximum(1);
    hexcl_exp[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hexcl_exp[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hexcl_exp[i]->GetZaxis()->SetTitle("expected excluded points");
    hexcl_exp[i]->Draw("colz");
    gr_exp[i]->Draw("l");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.15,0.92,"CMS Preliminary   #sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.6 fb^{-1}");
    */
    
    /*
    hjes[i]->GetYaxis()->SetTitle("#chi^{0}_{1} mass (GeV)");
    hjes[i]->GetXaxis()->SetTitle("#tilde{t} mass (GeV)");
    hjes[i]->GetZaxis()->SetTitle("JES uncertainty");
    hjes[i]->Draw("colz");

    t->DrawLatex(0.2,0.83,label);
    //t->DrawLatex(0.2,0.77,"m(#tilde{q}) >> m(#tilde{g})");
    t->DrawLatex(0.2,0.71,signames.at(i).c_str());
    t->DrawLatex(0.18,0.92,"CMS Preliminary            #sqrt{s} = 7 TeV, #scale[0.6]{#int}Ldt = 4.98 fb^{-1}");
    */

    if( print ){
      can[i]->          Print(Form("plots/%s%s_%s%s.eps"       ,sample,suffix,labels.at(i).c_str(),BDTchar));
      can[i]->          Print(Form("plots/%s%s_%s%s.pdf"       ,sample,suffix,labels.at(i).c_str(),BDTchar));
      //can_exclusion[i]->Print(Form("../../plots/%s%s_%s_points%s.pdf",sample,suffix,labels.at(i).c_str(),pol));
    }

    int bin = heff[i]->FindBin(300,50);

    if( TString(sample).Contains("T2bw") ) bin = heff[i]->FindBin(500,50);

    if( TString(sample).Contains("T2bw") ){ 
      cout << "efficiency (500,50)  " << heff[i]->GetBinContent(bin)      << endl;
    }
    else{
      cout << "efficiency (300,50)  " << heff[i]->GetBinContent(bin)      << endl;
    }
    cout << "nevents              " << hnevents[i]->GetBinContent(bin)  << endl;
    cout << "xsec UL              " << hxsec[i]->GetBinContent(bin)     << endl;
    cout << "xsec UL exp          " << hxsec_exp[i]->GetBinContent(bin) << endl;
    cout << "JES                  " << hjes[i]->GetBinContent(bin)      << endl;
    cout << "tot err              " << htoterr[i]->GetBinContent(bin)   << endl;
    //    cout << "btag err             " << hbtagerr[i]->GetBinContent(bin)  << endl;
    cout << "ISR err              " << hisrerr[i]->GetBinContent(bin)   << endl;
    cout << "stat err             " << hstaterr[i]->GetBinContent(bin)  << endl;
    cout << endl << endl;

    if( TString(sample).Contains("T2bw") ){ 
      *logfile << "efficiency (500,50)  " << heff[i]->GetBinContent(bin)      << endl;
    }
    else{
      *logfile << "efficiency (300,50)  " << heff[i]->GetBinContent(bin)      << endl;
    }
    *logfile << "nevents              " << hnevents[i]->GetBinContent(bin)  << endl;
    *logfile << "xsec UL              " << hxsec[i]->GetBinContent(bin)     << endl;
    *logfile << "xsec UL exp          " << hxsec_exp[i]->GetBinContent(bin) << endl;
    *logfile << "JES                  " << hjes[i]->GetBinContent(bin)      << endl;
    *logfile << "tot err              " << htoterr[i]->GetBinContent(bin)   << endl;
    //    *logfile << "btag err             " << hbtagerr[i]->GetBinContent(bin)  << endl;
    *logfile << "ISR err              " << hisrerr[i]->GetBinContent(bin)   << endl;
    *logfile << "stat err             " << hstaterr[i]->GetBinContent(bin)  << endl;
    *logfile << endl << endl;

  }


  TFile *outfile = TFile::Open(Form("rootfiles/%s%s%s_histos.root",sample,suffix,BDTchar),"RECREATE");

  outfile->cd();
  for( unsigned int i = 0 ; i < nsig ; ++i ){
    hxsec[i]->Write();
    hxsec_exp[i]->Write();
    hxsec_expp1[i]->Write();
    hxsec_expm1[i]->Write();
    hexcl[i]->Write();
    hexcl_obsp1[i]->Write();
    hexcl_obsm1[i]->Write();
    hexcl_exp[i]->Write();
    hexcl_expp1[i]->Write();
    hexcl_expm1[i]->Write();
    heff[i]->Write();
    hnsig[i]->Write();
    hnevents[i]->Write();
    hjes[i]->Write();
    htoterr[i]->Write();
    hisrerr[i]->Write();
    //    hbtagerr[i]->Write();
    hstaterr[i]->Write();
    // gr[i]->Write();
    // gr_exp[i]->Write();
  }
  outfile->Close();

}


void doAll(){

  // SMS("T2bw_MG",75,true ,"T2BW_SS",true);
  // SMS("T2bw_MG",75,false,"T2BW_SS",true);

  // SMS("T2bw_MG",25,true ,"T2BW_SS",true);
  // SMS("T2bw_MG",25,false,"T2BW_SS",true);

  // SMS("T2bw_MG",50,true ,"T2BW_SS",true);
  // SMS("T2bw_MG",50,false,"T2BW_SS",true);



  //--------------------------
  // T2tt
  //--------------------------

  // SMS("T2tt", 1,false,""     ,true);
  // SMS("T2tt", 1,false,"left" ,true);
  // SMS("T2tt", 1,false,"right",true);

  // SMS("T2tt", 1,true,""     ,true);
  // SMS("T2tt", 1,true,"left" ,true);
  // SMS("T2tt", 1,true,"right",true);

  /*
  //--------------------------
  // T2bw madgraph
  //--------------------------

  SMS("T2bw_MG",25,true,"",true);
  SMS("T2bw_MG",50,true,"",true);
  SMS("T2bw_MG",75,true,"",true);

  SMS("T2bw_MG",25,false,"",true);
  SMS("T2bw_MG",50,false,"",true);
  SMS("T2bw_MG",75,false,"",true);

  //--------------------------
  // T2bw pythia
  //--------------------------

  SMS("T2bw",25,false,"",true);
  SMS("T2bw",50,false,"",true);
  SMS("T2bw",75,false,"",true);

  SMS("T2bw",25,true,"",true);
  SMS("T2bw",50,true,"",true);
  SMS("T2bw",75,true,"",true);

  //--------------------------
  // T2bw madgraph reweighted
  //--------------------------

  char* weights[9]={
    "T2BW_LR",
    "T2BW_LS",
    "T2BW_LL",
    "T2BW_SR",
    "T2BW_SS",
    "T2BW_SL",
    "T2BW_RR",
    "T2BW_RS",
    "T2BW_RL"
  };

  for( int i = 0 ; i < 9 ; i++ ){
    SMS("T2bw_MG",25,true,weights[i],true);
    SMS("T2bw_MG",50,true,weights[i],true);
    SMS("T2bw_MG",75,true,weights[i],true);

    SMS("T2bw_MG",25,false,weights[i],true);
    SMS("T2bw_MG",50,false,weights[i],true);
    SMS("T2bw_MG",75,false,weights[i],true);
  }
  */

  // SMS("T2bw_MG",25,true,"T2BW_SS",true);
  // SMS("T2bw_MG",50,true,"T2BW_SS",true);
  // SMS("T2bw_MG",75,true,"T2BW_SS",true);

  // SMS("T2bw_MG",25,false,"T2BW_SS",true);
  // SMS("T2bw_MG",50,false,"T2BW_SS",true);
  // SMS("T2bw_MG",75,false,"T2BW_SS",true);

}

void doAllPlots(){

  // //--------------------------
  // // T2bw MG
  // //--------------------------

  // char* weights[5]={
  //   "T2BW_LR",
  //   "T2BW_LL",
  //   "T2BW_SS",
  //   "T2BW_RR",
  //   "T2BW_RL"
  // };

  // for( int i = 0 ; i < 5 ; i++ ){
  //   SMS("T2bw_MG",25,true,weights[i],true);
  //   SMS("T2bw_MG",50,true,weights[i],true);
  //   SMS("T2bw_MG",75,true,weights[i],true);

  //   SMS("T2bw_MG",25,false,weights[i],true);
  //   SMS("T2bw_MG",50,false,weights[i],true);
  //   SMS("T2bw_MG",75,false,weights[i],true);
  // }

  // //--------------------------
  // // T2tt
  // //--------------------------

  // SMS("T2tt", 1,false,""     ,true);
  // SMS("T2tt", 1,false,"left" ,true);
  // SMS("T2tt", 1,false,"right",true);

  // SMS("T2tt", 1,true,""     ,true);
  // SMS("T2tt", 1,true,"left" ,true);
  // SMS("T2tt", 1,true,"right",true);

}
