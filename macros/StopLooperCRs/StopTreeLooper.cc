#include "StopTreeLooper.h"

//#include "../../CORE/jetSmearingTools.h"
//#include "../../CORE/Thrust.h"
//#include "../../CORE/EventShape.h"
#include "TStopwatch.h"

#include "Math/VectorUtil.h"
#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"
#include "../Core/MT2Utility.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"
#include "../Core/PartonCombinatorics.h"

#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TChainElement.h"
#include "TFile.h"
#include "TMath.h"
#include "TChain.h"
#include "Riostream.h"
#include "TFitter.h"
#include "TRandom.h"

#include <algorithm>
#include <utility>
#include <map>
#include <set>
#include <list>

using namespace Stop;

std::set<DorkyEventIdentifier> already_seen; 
std::set<DorkyEventIdentifier> events_lasercalib; 
std::set<DorkyEventIdentifier> events_hcallasercalib; 

StopTreeLooper::StopTreeLooper()
{
  m_outfilename_ = "histos.root";
  t1metphicorr = -9999.;
  t1metphicorrphi = -9999.;
  t1metphicorrmt = -9999.;
  dphimjmin = -9999.;
  dphimj1 = -9999.;
  dphimj2 = -9999.;
  pt_b = -9999.;
  htssl = -9999.;
  htosl = -9999.;
  htratiol = -9999.;
  htssm = -9999.;
  htosm = -9999.;
  htratiom = -9999.;
  min_mtpeak = -9999.;
  max_mtpeak = -9999.; 
  n_jets  = -9999;
  n_bjets = -9999;
  n_ljets = -9999;
  chi2min_ = -9999.;
  chi2minprob_ = -9999.;
  mt2bmin_ = -9999.;
  mt2blmin_ = -9999.;
  mt2wmin_ = -9999.;

}

StopTreeLooper::~StopTreeLooper()
{
}

void StopTreeLooper::setOutFileName(string filename)
{
  m_outfilename_ = filename;

}

void StopTreeLooper::loop(TChain *chain, TString name)
{

  TStopwatch stwatch;
  //------------------------------
  // check for valid chain
  //------------------------------

  printf("[StopTreeLooper::loop] %s\n", name.Data());

  load_badlaserevents("../Core/badlaser_events.txt", events_lasercalib);
  load_badlaserevents("../Core/badhcallaser_events.txt", events_hcallasercalib);

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  if (listOfFiles->GetEntries() == 0) {
    cout << "[StopTreeLooper::loop] no files in chain" << endl;
    return;
  }

  //------------------------------
  // set up histograms
  //------------------------------

  gROOT->cd();

  cout << "[StopTreeLooper::loop] setting up histos" << endl;

  std::map<std::string, TH1F*> h_1d;
  //also for control regions
  std::map<std::string, TH1F*> h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr5;
  //for signal region 
  std::map<std::string, TH1F*> h_1d_sig;
  //for ttbar dilepton njets distribution
  std::map<std::string, TH1F*> h_1d_nj;
  //z sample for yields etc
  std::map<std::string, TH1F*> h_1d_z;

  //------------------------------
  // vtx reweighting
  //------------------------------

  // TFile* vtx_file = TFile::Open("../vtxreweight/vtxreweight_Summer12_DR53X-PU_S10_9p7ifb_Zselection.root");
  // if( vtx_file == 0 ){
  //   cout << "vtxreweight error, couldn't open vtx file. Quitting!"<< endl;
  //   exit(0);
  // }

  // TH1F* h_vtx_wgt = (TH1F*)vtx_file->Get("hratio");
  // h_vtx_wgt->SetName("h_vtx_wgt");

  //------------------------------
  // file loop
  //------------------------------

  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  ULong64_t nEventsTotal = 0;

  bool isData = name.Contains("data") ? true : false;

  //Define peak region 
  min_mtpeak = 50.; 
  max_mtpeak = 80.; 
  printf("[StopTreeLooper::loop] MT PEAK definition %.0f - %.0f GeV \n", min_mtpeak, max_mtpeak);

  cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

  stwatch.Start();

  while (TChainElement *currentFile = (TChainElement*)fileIter.Next()) {
  
    //----------------------------
    // load the stop baby tree
    //----------------------------

    TFile *file = new TFile( currentFile->GetTitle() );
    TTree *tree = (TTree*)file->Get("t");
    stopt.Init(tree);

    //----------------------------
    // event loop
    //----------------------------

    ULong64_t nEvents = tree->GetEntriesFast();
    for(ULong64_t event = 0; event < nEvents; ++event) {
      stopt.GetEntry(event);

      //----------------------------
      // increment counters
      //----------------------------

      ++nEventsTotal;
      if (nEventsTotal%10000==0) {
	ULong64_t i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
	//if (i_permille != i_permille_old) {//this prints too often!
	// xterm magic from L. Vacavant and A. Cerri
	if (isatty(1)) {
	  printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
		 "\033[0m\033[32m <---\033[0m\015", i_permille/10.);
	  fflush(stdout);
	  stwatch.Stop();
	  if (i_permille%100<0.0001)
	    cout<<"At "<<i_permille/10.<<"% time is "<<stwatch.RealTime()<<" cpu: "<<stwatch.CpuTime()<<endl;
	  stwatch.Start();//stwatch.Continue();
	  
	}
	//	i_permille_old = i_permille;
      }

      //---------------------
      // skip duplicates
      //---------------------

      if( isData ) {
        DorkyEventIdentifier id = {stopt.run(), stopt.event(), stopt.lumi() };
        if (is_duplicate(id, already_seen) ){
          continue;
        }
	if (is_badLaserEvent(id,events_lasercalib) ){
	  //std::cout<< "Removed bad laser calibration event:" << run << "   " << event<<"\n";
	  continue;
	}
	if (is_badLaserEvent(id,events_hcallasercalib) ){
	  //std::cout<< "Removed bad hcal laser calibration event:" << run << "   " << event<<"\n";
	  continue;
	}
      }

      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      float evtweight = isData ? 1. : 
	( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );
      // to reweight from file - also need to comment stuff before
      //      float vtxweight = vtxweight_n( nvtx, h_vtx_wgt, isData );

      plot1D("h_vtx",       stopt.nvtx(),       evtweight, h_1d, 40, 0, 40);
      plot1D("h_vtxweight", stopt.nvtxweight(), evtweight, h_1d, 41, -4., 4.);

      //----------------------------------------------------------------------------
      // apply preselection:
      // rho 0-40 GeV, MET filters, >=1 good lepton, veto 2 leptons dR < 0.1
      //----------------------------------------------------------------------------

      if ( !passEvtSelection(name) ) continue;

      //----------------------------------------------------------------------------
      // Function to perform MET phi corrections on-the-fly
      // Default branches are: tree->t1metphicorr_ and tree->t1metphicorrmt_
      //----------------------------------------------------------------------------

      pair<float, float> p_t1metphicorr = 
      	getPhiCorrMET( stopt.t1met10(), stopt.t1met10phi(), stopt.nvtx(), !isData);
      t1metphicorr    = p_t1metphicorr.first;
      t1metphicorrphi = p_t1metphicorr.second;
      t1metphicorrmt  = getMT( stopt.lep1().Pt() , stopt.lep1().Phi() , t1metphicorr , t1metphicorrphi );  

      //----------------------------------------------------------------------------
      // get jet information
      //----------------------------------------------------------------------------

      jets.clear();
      btag.clear();
      sigma_jets.clear();
      mc.clear();
      n_jets  = 0;
      n_bjets = 0;
      n_ljets = 0;
      htssl = 0.;
      htosl = 0.;
      htssm = 0.;
      htosm = 0.;

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){
	
	if( stopt.pfjets().at(i).pt()<30 )  continue;
	if( fabs(stopt.pfjets().at(i).eta())>2.4 )  continue;
	
	n_jets++;
	n_ljets++;
      	jets.push_back( stopt.pfjets().at(i)    );
	if (stopt.pfjets_csv().at(i) > 0.679) n_bjets++;
      	btag.push_back( stopt.pfjets_csv().at(i) );

      	if ( !isData ) mc.push_back  ( stopt.pfjets_mc3().at(i) );
      	else mc.push_back  ( 0 );

      	float sigma_i = stopt.pfjets_sigma().at(i);
        if ( isData ) sigma_i *= getDataMCRatio(stopt.pfjets().at(i).eta());
      	sigma_jets.push_back(sigma_i);

	float dPhiL = getdphi(stopt.lep1().Phi(), stopt.pfjets().at(i).phi() );
	float dPhiM = getdphi(t1metphicorr, stopt.pfjets().at(i).phi() );    
	if(dPhiL<(3.14/2))  htssl += stopt.pfjets().at(i).pt();
	if(dPhiL>=(3.14/2)) htosl += stopt.pfjets().at(i).pt();
	if(dPhiM<(3.14/2))  htssm += stopt.pfjets().at(i).pt();
	if(dPhiM>=(3.14/2)) htosm += stopt.pfjets().at(i).pt();

        //count jets that are not overlapping with second lepton
	if (isData) continue;
	if (stopt.nleps()!=2) continue;
	if (stopt.mclep2().pt() < 30.) continue;
	if (ROOT::Math::VectorUtil::DeltaR(stopt.mclep2(), stopt.pfjets().at(i)) > 0.4 ) continue;
	n_ljets--;

      } 

      dphimjmin= (n_jets>1) ? getMinDphi(t1metphicorrphi, jets.at(0),jets.at(1)) : -9999.;
      //b-pt 
      vector<int> indexBJets=getBJetIndex(0.679,-1,-1);
      if(indexBJets.size()>0) pt_b = stopt.pfjets().at(indexBJets.at(0)).pt();

      //maria variables
      htratiol = htssl / (htosl + htssl);
      htratiom = htssm / (htosm + htssm);

      // get list of candidates 
     PartonCombinatorics pc (stopt.pfjets(), stopt.pfjets_csv(), stopt.pfjets_sigma(), stopt.pfjets_mc3(), 
			     stopt.lep1(), t1metphicorr, t1metphicorrphi, isData);
     MT2CHI2 mt2c2 = pc.getMt2Chi2();

     // chi2 and MT2 variables
     chi2min_= mt2c2.one_chi2;               // minimum chi2 
     chi2minprob_= TMath::Prob(chi2min_,1);   // probability of minimum chi2

     mt2bmin_= mt2c2.three_mt2b;             // minimum MT2b
     mt2blmin_= mt2c2.three_mt2bl;            // minimum MT2bl
     mt2wmin_= mt2c2.three_mt2w;             // minimum MT2w

      //----------------------------------------------------------------------------
      // histogram tags
      //----------------------------------------------------------------------------

      //jet multiplicity - inclusive above 4
      string tag_njets = Form("_nj%i", (n_jets<4) ? n_jets : 4);

      //b-tagging
      string tag_btag = (n_bjets<1) ? "_bveto" : "";

      //iso-trk-veto
      string tag_isotrk = passIsoTrkVeto() ? "" : "_wisotrk";

      // tag_mt2w
      string tag_mt2w = (mt2wmin_<175) ? "_failmt2w" : "_mt2w"; 

      // tag_chi2
      string tag_chi2 = (chi2minprob_>0.1) ? "_failchi2" : "_chi2";

      //z-peak/veto
      string tag_zcut;
      if ( fabs( stopt.dilmass() - 91.) > 15. ) tag_zcut = "_zveto";
      else if  ( fabs( stopt.dilmass() - 91.) < 10. ) tag_zcut = "_zpeak";
      else tag_zcut = "_ignore";

      //Corrected jet counting with simple overlap removal
      string tag_kbin = (n_ljets<4) ? "_K3" : "_K4";

      //event with true truth-level track
      bool hastruetrk = false;
      if (stopt.nleps()==2 && abs(stopt.mclep2().Eta())<2.5)  {
	//check if second lepton is e/mu pT>10GeV
	if (abs(stopt.mcid2())<14 && stopt.mclep2().Pt()>10.) hastruetrk = true;
	//if second lepton is tau 
	//check if daughter lepton or single track has pT>10GeV
	if (abs(stopt.mcid2())>14 && stopt.mctaudpt2()>10.) {  
	  if (stopt.mcdecay2()==2) hastruetrk = true;
	  if (stopt.mcdecay2()==1 && stopt.mcndec2()==1) hastruetrk = true;
	}
      }
      string tag_truetrk = hastruetrk ? "_wtruetrk" : "_notruetrk";

      //flavor types
      string flav_tag_sl;
      if ( abs(stopt.id1())==13 ) flav_tag_sl = "_muo";
      else if ( abs(stopt.id1())==11 ) flav_tag_sl = "_ele";
      else flav_tag_sl = "_mysterysl";
      string flav_tag_dl;
      if      ( abs(stopt.id1()) == abs(stopt.id2()) && abs(stopt.id1()) == 13 ) 
	flav_tag_dl = "_dimu";
      else if ( abs(stopt.id1()) == abs(stopt.id2()) && abs(stopt.id1()) == 11 ) 
	flav_tag_dl = "_diel";
      else if ( abs(stopt.id1()) != abs(stopt.id2()) && abs(stopt.id1()) == 13 ) 
	flav_tag_dl = "_muel";
      else if ( abs(stopt.id1()) != abs(stopt.id2()) && abs(stopt.id1()) == 11 ) 
	flav_tag_dl = "_elmu";
      else flav_tag_dl = "_mysterydl";
      string basic_flav_tag_dl = flav_tag_dl;
      if ( abs(stopt.id1()) != abs(stopt.id2()) ) basic_flav_tag_dl = "_mueg";

      //
      // SIGNAL REGION - single lepton + b-tag
      //

      // selection - 1 lepton 
      // Add iso track veto
      // Add b-tag
      if ( passSingleLeptonSelection(isData) && n_jets>=4 )
	{
	  float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
	  //default 
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag", 	         tag_kbin, flav_tag_sl, 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag,   	         tag_kbin, flav_tag_sl, 150. );
	  makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+tag_truetrk, tag_kbin, flav_tag_sl, 150. );
	  if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+tag_mt2w,   	    tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+tag_chi2,   	    tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 150. );

	  //met > 50 GeV requirement 
	  if ( t1metphicorr > 50. ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met50",  	    tag_kbin, flav_tag_sl, 150. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met50", 	    tag_kbin, flav_tag_sl, 150. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met50"+tag_truetrk, tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met50"+tag_mt2w, 	       tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met50"+tag_chi2, 	       tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met50"+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 150. );

	  }
	  //met > 100 GeV requirement 
	  if ( t1metphicorr > 100. ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met100",  		tag_kbin, flav_tag_sl, 150. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met100", 		tag_kbin, flav_tag_sl, 150. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met100"+tag_truetrk, tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met100"+tag_mt2w, 	tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met100"+tag_chi2, 	tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met100"+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 150. );

	  }
	  //met > 150 GeV requirement 
	  if ( t1metphicorr > 150. ) { 
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met150",  	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met150", 	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met150"+tag_truetrk, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met150"+tag_mt2w, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met150"+tag_chi2, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met150"+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 120. );

	  }
	  //met > 200 GeV requirement 
	  if ( t1metphicorr > 200. ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met200",  	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met200", 	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met200"+tag_truetrk, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met200"+tag_mt2w, 	 tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met200"+tag_chi2, 	 tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met200"+tag_chi2+tag_mt2w, tag_kbin, flav_tag_sl, 120. );
	  }

	  //met > 250 GeV requirement 
	  if ( t1metphicorr > 250. ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met250",  	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met250", 	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met250"+tag_truetrk, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met250"+tag_mt2w, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met250"+tag_chi2, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met250"+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 300 GeV requirement 
	  if ( t1metphicorr > 300. ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met300",  	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met300", 	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met300"+tag_truetrk, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met300"+tag_mt2w, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met300"+tag_chi2, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met300"+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 350 GeV requirement 
	  if ( t1metphicorr > 350. ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met350",  	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met350", 	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met350"+tag_truetrk, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met350"+tag_mt2w, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met350"+tag_chi2, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met350"+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 120. );

	  }
	  //met > 400 GeV requirement 
	  if ( t1metphicorr > 400. ) {
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+"_prebtag_met400",  	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met400", 	     tag_kbin, flav_tag_sl, 120. );
	    makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met400"+tag_truetrk, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met400"+tag_mt2w, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met400"+tag_chi2, 	tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeSIGPlots( evtweight*trigweight, h_1d_sig, tag_isotrk+tag_btag+"_met400"+tag_chi2+tag_mt2w,tag_kbin, flav_tag_sl, 120. );
	  }

	}

      //
      // CR1 - single lepton + b-veto
      //

      // selection - 1 lepton + iso track veto
      // Add b-tag veto
      if ( passOneLeptonSelection(isData) )
	{

	  float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
	  //default 
	  makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto", tag_njets, tag_kbin, flav_tag_sl, 150. );
	  if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w, tag_njets, tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2, tag_njets, tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w, tag_njets, tag_kbin, flav_tag_sl, 150. );
	  //met > 50 GeV requirement 
	  if ( t1metphicorr > 50. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met50", tag_njets, tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met50", tag_njets, tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met50", tag_njets, tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met50", tag_njets, tag_kbin, flav_tag_sl, 150. );
	  }
	  //met > 100 GeV requirement 
	  if ( t1metphicorr > 100. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met100", tag_njets, tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met100", tag_njets, tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met100", tag_njets, tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met100", tag_njets, tag_kbin, flav_tag_sl, 150. );
	  }
	  //met > 150 GeV requirement 
	  if ( t1metphicorr > 150. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met150", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met150", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met150", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met150", tag_njets, tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 200 GeV requirement 
	  if ( t1metphicorr > 200. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met200", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met200", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met200", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met200", tag_njets, tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 250 GeV requirement 
	  if ( t1metphicorr > 250. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met250", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met250", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met250", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met250", tag_njets, tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 300 GeV requirement 
	  if ( t1metphicorr > 300. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met300",                      tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met300",   	   tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met300", 	   tag_njets, tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met300", tag_njets, tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 350 GeV requirement 
	    if ( t1metphicorr > 350. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met350", tag_njets, tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met350", tag_njets, tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met350", tag_njets, tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met350", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 400 GeV requirement 
	    if ( t1metphicorr > 400. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto_met400", tag_njets, tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_mt2w+"_met400", tag_njets, tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+"_met400", tag_njets, tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto"+tag_chi2+tag_mt2w+"_met400", tag_njets, tag_kbin, flav_tag_sl, 120. );
	    }
	  
	  if ( n_bjets==0 ) {

	    //default 
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1,                "", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w, tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2, tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w, tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    //met > 50 GeV requirement 
	    if ( t1metphicorr > 50. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met50", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met50", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met50", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met50", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    }
	    //met > 100 GeV requirement 
	    if ( t1metphicorr > 100. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met100", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	     if (mt2wmin_>175)                      makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met100", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	     if (chi2minprob_<0.1)                  makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met100", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	     if (chi2minprob_<0.1 && mt2wmin_>175)  makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met100", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    }
	    //met > 150 GeV requirement 
	    if ( t1metphicorr > 150. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met150", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met150", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met150", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met150", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 200 GeV requirement 
	    if ( t1metphicorr > 200. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 250 GeV requirement 
	    if ( t1metphicorr > 250. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met250", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met250", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met250", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met250", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 300 GeV requirement 
	    if ( t1metphicorr > 300. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met300", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met300", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met300", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met300", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 350 GeV requirement 
	    if ( t1metphicorr > 350. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met350", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met350", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met350", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met350", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 400 GeV requirement 
	    if ( t1metphicorr > 400. ) {
	      makeCR1Plots( evtweight*trigweight, h_1d_cr1,                   "_met400", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_mt2w+"_met400", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr1,          tag_chi2+"_met400", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr1, tag_chi2+tag_mt2w+"_met400", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	  }
	}

      //
      // CR2 - Z-peak for yields and mT resolution studies
      // 

      // selection - SF dilepton, veto on isolated track in addition to 2 leptons, in z-peak
      //if ( passTwoLeptonSelection(isData) )
      if ( passDileptonSelection(isData) )
	{
	  float trigweight = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());
	  //invariant mass - basic check of inclusive distribution
	  plot1D("h_z_dilmass"          +flav_tag_dl, stopt.dilmass(), evtweight*trigweight, h_1d_z,  30 , 76 , 106);
	  plot1D("h_z_dilmass"+tag_njets+flav_tag_dl, stopt.dilmass(), evtweight*trigweight, h_1d_z,  30 , 76 , 106);

	  if ( fabs( stopt.dilmass() - 91.) < 10. ) 
	    {

	      // if (n_jets>8) 
	      // 	cout<<"NJETS: "<<n_jets<<" * dataset: "<<stopt.dataset()
	      // 	    <<" run: "<<stopt.run()<<" lumi: "<<stopt.lumi()<<" event: "<<stopt.event()<<endl;
	      
	      //z peak plots
	      plot1D("h_z_njets"    +flav_tag_dl, min(n_jets,4),  evtweight*trigweight, h_1d_z, 5,0,5);
	      plot1D("h_z_njets_all"+flav_tag_dl, min(n_jets,9),  evtweight*trigweight, h_1d_z, 10, 0, 10);
	      plot1D("h_z_nbjets"   +flav_tag_dl, min(n_bjets,3), evtweight*trigweight, h_1d_z, 4, 0, 4);
	      makeZPlots( evtweight*trigweight, h_1d_z, "", tag_njets, flav_tag_dl );

	      // Add stricter 3rd lepton veto
	      // require at least 2 jets
	      if ( (stopt.trkpt10loose() <0.0001 || stopt.trkreliso10loose() > 0.1) && n_jets >= 2 ) {

		//find positive lepton
		bool isfirstp = (stopt.id1() > 0) ? true : false;
		
		//recalculate met
		float metx = t1metphicorr * cos( t1metphicorrphi );
		float mety = t1metphicorr * sin( t1metphicorrphi );
		
		//recalculate the MET with the positive lepton
		metx += isfirstp ? stopt.lep1().px() : stopt.lep2().px();
		mety += isfirstp ? stopt.lep1().py() : stopt.lep2().py();
		
		float t1metphicorr_lep    = sqrt(metx*metx + mety*mety);
		
		//default 
		makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto", tag_njets,  tag_kbin, basic_flav_tag_dl, 150. );
		//pseudomet > 50 GeV requirement 
		if ( t1metphicorr_lep > 50. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met50", tag_njets,  tag_kbin, basic_flav_tag_dl, 150. );
		//pseudomet > 100 GeV requirement 
		if ( t1metphicorr_lep > 100. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met100", tag_njets,  tag_kbin, basic_flav_tag_dl, 150. );
		//pseudomet > 150 GeV requirement 
		if ( t1metphicorr_lep > 150. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met150", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		//pseudomet > 200 GeV requirement 
		if ( t1metphicorr_lep > 200. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met200", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		//pseudomet > 250 GeV requirement 
		if ( t1metphicorr_lep > 250. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met250", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		//pseudomet > 300 GeV requirement 
		if ( t1metphicorr_lep > 300. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met300", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		//pseudomet > 350 GeV requirement 
		if ( t1metphicorr_lep > 350. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met350", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		//pseudomet > 400 GeV requirement 
		if ( t1metphicorr_lep > 400. ) 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_prebveto_met400", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );

		// Add b-tag veto 
		if ( n_bjets==0) {
		  //default 
		  makeCR2Plots( evtweight*trigweight, h_1d_cr2, "", tag_njets,  tag_kbin, basic_flav_tag_dl, 150. );
		  //pseudomet > 50 GeV requirement 
		  if ( t1metphicorr_lep > 50. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_met50", tag_njets,  tag_kbin, basic_flav_tag_dl, 150. );
		  //pseudomet > 100 GeV requirement 
		  if ( t1metphicorr_lep > 100. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_met100", tag_njets,  tag_kbin, basic_flav_tag_dl, 150. );
		  //pseudomet > 150 GeV requirement 
		  if ( t1metphicorr_lep > 150. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_met150", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		  //pseudomet > 200 GeV requirement 
		  if ( t1metphicorr_lep > 200. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_met200", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		  //pseudomet > 250 GeV requirement 
		  if ( t1metphicorr_lep > 250. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2, "_met250", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		  //pseudomet > 300 GeV requirement 
		  if ( t1metphicorr_lep > 300. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2,"_met300", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		  //pseudomet > 350 GeV requirement 
		  if ( t1metphicorr_lep > 350. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2,"_met350", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		  //pseudomet > 400 GeV requirement 
		  if ( t1metphicorr_lep > 400. ) 
		    makeCR2Plots( evtweight*trigweight, h_1d_cr2,"_met400", tag_njets,  tag_kbin, basic_flav_tag_dl, 120. );
		}
	      }
	    }
	}

      //
      // CR4 - ttbar dilepton sample with 2 good leptons
      //
      
      // selection - all dilepton, z-veto for SF dilepton
      // Add b-tag requirement
      if ( passDileptonSelection(isData) 
	   && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. ) 
	   && n_bjets>0 ) 
	{
	  float trigweight = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

	  //jet multiplicity distributions 
	  //store in separate file since this is used for njet reweighting
	  makeNJPlots( evtweight*trigweight, h_1d_nj, "", basic_flav_tag_dl);
	  //met > 50 GeV requirement 
	  if ( t1metphicorr > 50. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met50", basic_flav_tag_dl);
	  //met > 100 GeV requirement 
	  if ( t1metphicorr > 100. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met100", basic_flav_tag_dl);
	  //met > 150 GeV requirement 
	  if ( t1metphicorr > 150. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met150", basic_flav_tag_dl);	    
	  //met > 200 GeV requirement 
	  if ( t1metphicorr > 200. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met200", basic_flav_tag_dl);	    
	  //met > 250 GeV requirement 
	  if ( t1metphicorr > 250. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met250", basic_flav_tag_dl);	    
	  //met > 300 GeV requirement 
	  if ( t1metphicorr > 300. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met300", basic_flav_tag_dl);	    
	  //met > 350 GeV requirement 
	  if ( t1metphicorr > 350. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met350", basic_flav_tag_dl);	    
	  //met > 400 GeV requirement 
	  if ( t1metphicorr > 400. ) 
	    makeNJPlots( evtweight*trigweight, h_1d_nj, "_met400", basic_flav_tag_dl);	    

	  if ( n_jets < 2 ) continue; 
	  
	  //default 
	  makeCR4Plots( evtweight*trigweight, h_1d_cr4, "",                tag_njets,  tag_kbin, flav_tag_dl, 150. );
	  if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w,          tag_njets,  tag_kbin, flav_tag_dl, 150. );
	  if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2,          tag_njets,  tag_kbin, flav_tag_dl, 150. );
	  if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w, tag_njets,  tag_kbin, flav_tag_dl, 150. );
	  //met > 50 GeV requirement 
	  if ( t1metphicorr > 50. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met50",                   tag_njets,  tag_kbin, flav_tag_dl, 150. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met50",          tag_njets,  tag_kbin, flav_tag_dl, 150. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met50",          tag_njets,  tag_kbin, flav_tag_dl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met50", tag_njets,  tag_kbin, flav_tag_dl, 150. );
	  }
	  //met > 100 GeV requirement 
	  if ( t1metphicorr > 100. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met100",                   tag_njets,  tag_kbin, flav_tag_dl, 150. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met100",          tag_njets,  tag_kbin, flav_tag_dl, 150. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met100",          tag_njets,  tag_kbin, flav_tag_dl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met100", tag_njets,  tag_kbin, flav_tag_dl, 150. );
	  }
	  //met > 150 GeV requirement 
	  if ( t1metphicorr > 150. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met150",                   tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met150",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met150",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met150", tag_njets,  tag_kbin, flav_tag_dl, 120. );
	  }
	  //met > 200 GeV requirement 
	  if ( t1metphicorr > 200. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met200",                   tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met200",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met200",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met200", tag_njets,  tag_kbin, flav_tag_dl, 120. );
	  }
	  //met > 250 GeV requirement 
	  if ( t1metphicorr > 250. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met250",                   tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met250",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met250",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met250", tag_njets,  tag_kbin, flav_tag_dl, 120. );
	  }
	  //met > 300 GeV requirement 
	  if ( t1metphicorr > 300. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met300",                   tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met300",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met300",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met300", tag_njets,  tag_kbin, flav_tag_dl, 120. );
	  }
	  //met > 350 GeV requirement 
	  if ( t1metphicorr > 350. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met350",                   tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met350",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met350",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met350", tag_njets,  tag_kbin, flav_tag_dl, 120. );
	  }
	  //met > 400 GeV requirement 
	  if ( t1metphicorr > 400. ) {
	    makeCR4Plots( evtweight*trigweight, h_1d_cr4, "_met400",                   tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (mt2wmin_>175)                     makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_mt2w+"_met400",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1)                 makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+"_met400",          tag_njets,  tag_kbin, flav_tag_dl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR4Plots( evtweight*trigweight, h_1d_cr4, tag_chi2+tag_mt2w+"_met400", tag_njets,  tag_kbin, flav_tag_dl, 120. );
	  }
	}

      ////////////////////////////////////////////////////////////////////////////////////////////////////
      // Ask for at least 2 jets from now on
      if ( n_jets < 2 ) continue;


      //
      // Sample before isolated track requirement - for fake rate of requirement
      //

      // selection - at least 1 lepton
      // Add b-tag requirement
      if ( passSingleLeptonSelection(isData) 
	   && n_bjets>0 ) 
	{
	  float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
	  //inclusive sample
	  //default 
	  makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto",                   tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w,          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2,          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w, tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  //met > 50 GeV requirement  
	  if ( t1metphicorr > 50. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met50",                      tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met50",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met50",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_met50", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  }
	  //met > 100 GeV requirement 
	  if ( t1metphicorr > 100. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met100",                      tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met100",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met100",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_met100", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  }
	  //met > 150 GeV requirement 
	  if ( t1metphicorr > 150. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met150",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met150",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met150",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_met150", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 200 GeV requirement 
	  if ( t1metphicorr > 200. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met200",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met200",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met200",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_met200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 250 GeV requirement 
	  if ( t1metphicorr > 250. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met250",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met250",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met250",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_met250", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 300 GeV requirement 
	  if ( t1metphicorr > 300. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met300",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met300",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met300",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_me3200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 350 GeV requirement 
	  if ( t1metphicorr > 350. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met350",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met350",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met350",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_met350", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 400 GeV requirement 
	  if ( t1metphicorr > 400. ) {
	    makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto_met400",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_mt2w+"_met400",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+"_met400",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR1Plots( evtweight*trigweight, h_1d_cr5, "_preveto"+tag_chi2+tag_mt2w+"_met400", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	}
      
      //
      // CR5 - lepton + isolated track
      //

      // selection - lepton + isolated track
      // Add b-tag requirement
      if ( passLepPlusIsoTrkSelection(isData) 
	   && n_bjets>0 ) 
	{
	  float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
	  //inclusive sample
	  //default 
	  makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all",                   tag_njets, tag_kbin, flav_tag_sl, 150. );
	  if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w,          tag_njets, tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2,          tag_njets, tag_kbin, flav_tag_sl, 150. );
	  if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w, tag_njets, tag_kbin, flav_tag_sl, 150. );
	  //met > 50 GeV requirement 
	  if ( t1metphicorr > 50. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met50",                      tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met50",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met50",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met50", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  }
	  //met > 100 GeV requirement 
	  if ( t1metphicorr > 100. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met100",                      tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met100",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met100",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met100", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	  }
	  //met > 150 GeV requirement 
	  if ( t1metphicorr > 150. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met150",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met150",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met150",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met150", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 200 GeV requirement 
	  if ( t1metphicorr > 200. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met200",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met200",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met200",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 250 GeV requirement 
	  if ( t1metphicorr > 250. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met250",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met250",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met250",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met250", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 300 GeV requirement 
	  if ( t1metphicorr > 300. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met300",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met300",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met300",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met300", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 350 GeV requirement 
	  if ( t1metphicorr > 350. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met350",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met350",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met350",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met350", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }
	  //met > 400 GeV requirement 
	  if ( t1metphicorr > 400. ) {
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all_met400",                      tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_mt2w+"_met400",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+"_met400",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_all"+tag_chi2+tag_mt2w+"_met400", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	  }

	  // sample with only 1 lepton - this is the true CR5
	  if ( stopt.ngoodlep() == 1 ) {

	    //default 
	    makeCR5Plots( evtweight*trigweight, h_1d_cr5, "",                tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w,          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2,          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w, tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    //met > 50 GeV requirement 
	    if ( t1metphicorr > 50. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met50",                   tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met50",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met50",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met50", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    }
	    //met > 100 GeV requirement 
	    if ( t1metphicorr > 100. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met100",                   tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met100",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met100",          tag_njets,  tag_kbin, flav_tag_sl, 150. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met100", tag_njets,  tag_kbin, flav_tag_sl, 150. );
	    }
	    //met > 150 GeV requirement 
	    if ( t1metphicorr > 150. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met150",                   tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met150",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met150",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met150", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 200 GeV requirement 
	    if ( t1metphicorr > 200. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met200",                   tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met200",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met200",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met200", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 250 GeV requirement 
	    if ( t1metphicorr > 250. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met250",                   tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met250",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met250",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met250", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 300 GeV requirement 
	    if ( t1metphicorr > 300. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met300",                   tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met300",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met300",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met300", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 350 GeV requirement 
	    if ( t1metphicorr > 350. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met350",                   tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met350",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met350",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met350", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	    //met > 400 GeV requirement 
	    if ( t1metphicorr > 400. ) {
	      makeCR5Plots( evtweight*trigweight, h_1d_cr5, "_met400",                   tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (mt2wmin_>175)                     makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_mt2w+"_met400",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1)                 makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+"_met400",          tag_njets,  tag_kbin, flav_tag_sl, 120. );
	      if (chi2minprob_<0.1 && mt2wmin_>175) makeCR5Plots( evtweight*trigweight, h_1d_cr5, tag_chi2+tag_mt2w+"_met400", tag_njets,  tag_kbin, flav_tag_sl, 120. );
	    }
	  }
	}

    } // end event loop

    // delete tree;
   
    stwatch.Stop();
    cout<<"time is "<<stwatch.CpuTime()<<endl;
    stwatch.Start();

  } // end file loop
  
    //
    // finish
    //
  
  TFile outfile(m_outfilename_.c_str(),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());
  
  std::map<std::string, TH1F*>::iterator it1d;
  for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
    it1d->second->Write(); 
    delete it1d->second;
  }
  
  outfile.Write();
  outfile.Close();

  TFile outfile_sig(Form("SIG%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving SIG histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_sig;
  for(it1d_sig=h_1d_sig.begin(); it1d_sig!=h_1d_sig.end(); it1d_sig++) {
    it1d_sig->second->Write(); 
    delete it1d_sig->second;
  }

  outfile_sig.Write();
  outfile_sig.Close();

  //control regions
  //h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr5

  TFile outfile_cr1(Form("CR1%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR1 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr1;
  for(it1d_cr1=h_1d_cr1.begin(); it1d_cr1!=h_1d_cr1.end(); it1d_cr1++) {
    it1d_cr1->second->Write(); 
    delete it1d_cr1->second;
  }

  outfile_cr1.Write();
  outfile_cr1.Close();

  TFile outfile_cr2(Form("CR2%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR2 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr2;
  for(it1d_cr2=h_1d_cr2.begin(); it1d_cr2!=h_1d_cr2.end(); it1d_cr2++) {
    it1d_cr2->second->Write(); 
    delete it1d_cr2->second;
  }

  outfile_cr2.Write();
  outfile_cr2.Close();
    
  TFile outfile_cr4(Form("CR4%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR4 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr4;
  for(it1d_cr4=h_1d_cr4.begin(); it1d_cr4!=h_1d_cr4.end(); it1d_cr4++) {
    it1d_cr4->second->Write(); 
    delete it1d_cr4->second;
  }

  outfile_cr4.Write();
  outfile_cr4.Close();

  TFile outfile_cr5(Form("CR5%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR5 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr5;
  for(it1d_cr5=h_1d_cr5.begin(); it1d_cr5!=h_1d_cr5.end(); it1d_cr5++) {
    it1d_cr5->second->Write(); 
    delete it1d_cr5->second;
  }

  outfile_cr5.Write();
  outfile_cr5.Close();
    
  TFile outfile_nj(Form("NJ%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving NJ histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_nj;
  for(it1d_nj=h_1d_nj.begin(); it1d_nj!=h_1d_nj.end(); it1d_nj++) {
    it1d_nj->second->Write(); 
    delete it1d_nj->second;
  }

  outfile_nj.Write();
  outfile_nj.Close();

  TFile outfile_z(Form("Z%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving Z histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_z;
  for(it1d_z=h_1d_z.begin(); it1d_z!=h_1d_z.end(); it1d_z++) {
    it1d_z->second->Write(); 
    delete it1d_z->second;
  }

  outfile_z.Write();
  outfile_z.Close();

  already_seen.clear();

  gROOT->cd();

}



void StopTreeLooper::makeCR2Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag_dl, float mtcut ) 
{

  //find positive lepton - this is the one that is combined with the pseudomet to form the mT
  bool isfirstp = (stopt.id1() > 0) ? true : false;
  
  //recalculate met
  float metx = t1metphicorr * cos( t1metphicorrphi );
  float mety = t1metphicorr * sin( t1metphicorrphi );
          
  //recalculate the MET with the positive lepton
  metx += isfirstp ? stopt.lep1().px() : stopt.lep2().px();
  mety += isfirstp ? stopt.lep1().py() : stopt.lep2().py();
  
  float t1met10_lep    = sqrt(metx*metx + mety*mety);
  float t1met10phi_lep = atan2( mety , metx );
  
  //recalculate the MT with the negative lepton
  float t1met10mt_lep = isfirstp ?
    getMT( stopt.lep2().Pt() , stopt.lep2().Phi() , t1met10_lep , t1met10phi_lep ) :
    getMT( stopt.lep1().Pt() , stopt.lep1().Phi() , t1met10_lep , t1met10phi_lep );

  //binning for mT plots
  int nbins = 30;
  float h_xmin = 0.;
  float h_xmax = 300.;
  float x_ovflw = h_xmax-0.001;
  
  float pseudomt_count = -1.;
  if ( t1met10mt_lep > min_mtpeak
       && t1met10mt_lep < max_mtpeak )    pseudomt_count = 0.5;
  else if ( t1met10mt_lep > mtcut ) pseudomt_count = 1.5;
  
  //default met
  plot1D("h_cr2_met"+tag_selection          +flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_met"+tag_selection+tag_njets+flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_met"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //pseudo-met
  plot1D("h_cr2_pseudomet"+tag_selection	  +flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomet"+tag_selection+tag_njets+flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomet"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1met10_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //positive lepton pt - enters mT calculation
  float leppt = isfirstp ? stopt.lep1().Pt() : stopt.lep2().Pt();
  plot1D("h_cr2_leppt"+tag_selection          +flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_leppt"+tag_selection+tag_njets+flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(leppt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between pos-lep and pseudopseudomet
  float dphi_pseudometlep = isfirstp ?
    getdphi( stopt.lep2().Phi() , t1met10phi_lep ) :
    getdphi( stopt.lep1().Phi() , t1met10phi_lep );
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection          +flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection+tag_njets+flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr2_dphi_pseudometlep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, dphi_pseudometlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //pseudo-mt
  plot1D("h_cr2_pseudomt"      +tag_selection          +flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt"      +tag_selection+tag_njets+flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt"      +tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1met10mt_lep, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr2_pseudomt_count"+tag_selection          +flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr2_pseudomt_count"+tag_selection+tag_njets+flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr2_pseudomt_count"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, pseudomt_count, evtweight, h_1d, 2, 0, 2);
  //boson pt 
  float pt_boson = (stopt.lep1()+stopt.lep2()).pt();
  plot1D("h_cr2_pt_dilep"+tag_selection                   +flav_tag_dl, pt_boson, evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr2_pt_dilep"+tag_selection+tag_njets         +flav_tag_dl, pt_boson, evtweight, h_1d, 100, 0., 500);

}

void StopTreeLooper::makeCR4Plots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag_dl, float mtcut ) 
{
  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_cr4_met"+tag_selection                   +flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr4_met"+tag_selection+tag_njets         +flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr4_met"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  //leading lepton pt - enters mT calculation
  plot1D("h_cr4_leppt"+tag_selection                   +flav_tag_dl, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_leppt"+tag_selection+tag_njets         +flav_tag_dl, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //leading lepton eta
  plot1D("h_cr4_lepeta"+tag_selection                   +flav_tag_dl, stopt.lep1().Eta(), evtweight, h_1d, 21, -2.1, 2.1);
  plot1D("h_cr4_lepeta"+tag_selection+tag_njets         +flav_tag_dl, stopt.lep1().Eta(), evtweight, h_1d, 21, -2.1, 2.1);
  plot1D("h_cr4_lepeta"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, stopt.lep1().Eta(), evtweight, h_1d, 21, -2.1, 2.1);
  //subleading lepton pt
  plot1D("h_cr4_subleadleppt"+tag_selection                   +flav_tag_dl, min(stopt.lep2().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_subleadleppt"+tag_selection+tag_njets         +flav_tag_dl, min(stopt.lep2().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_subleadleppt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(stopt.lep2().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lead-lep and met
  float dphi_metlep = getdphi( stopt.lep1().Phi() , t1metphicorrphi );
  plot1D("h_cr4_dphi_metlep"+tag_selection                   +flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_metlep"+tag_selection+tag_njets         +flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_metlep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //min dphi leading two jets - this should cut most of the ttbar single leptons
  plot1D("h_cr4_mindPhiJ12"+tag_selection                   +flav_tag_dl, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_mindPhiJ12"+tag_selection+tag_njets         +flav_tag_dl, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_mindPhiJ12"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  //b-pT
  plot1D("h_cr4_bpt"+tag_selection                   +flav_tag_dl, pt_b, evtweight, h_1d, 50, 30., 400.);
  plot1D("h_cr4_bpt"+tag_selection+tag_njets         +flav_tag_dl, pt_b, evtweight, h_1d, 50, 30., 400.);
  plot1D("h_cr4_bpt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, pt_b, evtweight, h_1d, 50, 30., 400.);
  //maria variables
  plot1D("h_cr4_htratiom"+tag_selection                   +flav_tag_dl, htratiom, evtweight, h_1d, 50, 0., 1.);
  plot1D("h_cr4_htratiom"+tag_selection+tag_njets         +flav_tag_dl, htratiom, evtweight, h_1d, 50, 0., 1.);
  plot1D("h_cr4_htratiom"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, htratiom, evtweight, h_1d, 50, 0., 1.);
  /// MT2 and chi2
  plot1D("h_cr4_mt2wmin"+tag_selection         		 +flav_tag_dl, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2wmin"+tag_selection+tag_njets         +flav_tag_dl, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2wmin"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2bmin"+tag_selection                   +flav_tag_dl, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2bmin"+tag_selection+tag_njets         +flav_tag_dl, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2bmin"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2blmin"+tag_selection                   +flav_tag_dl, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2blmin"+tag_selection+tag_njets         +flav_tag_dl, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_mt2blmin"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_chi2minprob"+tag_selection                   +flav_tag_dl, chi2minprob_ , evtweight, h_1d, 100, 0., 1);
  plot1D("h_cr4_chi2minprob"+tag_selection+tag_njets         +flav_tag_dl, chi2minprob_ , evtweight, h_1d, 100, 0., 1);
  plot1D("h_cr4_chi2minprob"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, chi2minprob_ , evtweight, h_1d, 100, 0., 1);

  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_cr4_mt"+tag_selection                   +flav_tag_dl, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt"+tag_selection+tag_njets         +flav_tag_dl, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr4_mt_count"+tag_selection                   +flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr4_mt_count"+tag_selection+tag_njets         +flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr4_mt_count"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, mt_count, evtweight, h_1d, 2, 0, 2);

  //Plot more angles between various objects
  //angle between 2 leptons
  float dphi_dilep = getdphi( stopt.lep1().Phi() ,  stopt.lep2().Phi() );
  plot1D("h_cr4_dphi_dilep"+tag_selection                   +flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_dilep"+tag_selection+tag_njets         +flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr4_dphi_dilep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, dphi_dilep, evtweight, h_1d, 15, 0., TMath::Pi());
  //dR between 2 leptons
  float dR_dilep = dRbetweenVectors( stopt.lep1() ,  stopt.lep2() );
  plot1D("h_cr4_dR_dilep"+tag_selection                   +flav_tag_dl, min(dR_dilep, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr4_dR_dilep"+tag_selection+tag_njets         +flav_tag_dl, min(dR_dilep, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr4_dR_dilep"+tag_selection+tag_njets+tag_kbin+flav_tag_dl, min(dR_dilep, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  //boson pt (here the dilepton pt is proportional to ISR)
  float pt_boson = (stopt.lep1()+stopt.lep2()).pt();
  plot1D("h_cr4_pt_dilep"+tag_selection                   +flav_tag_dl, pt_boson, evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr4_pt_dilep"+tag_selection+tag_njets         +flav_tag_dl, pt_boson, evtweight, h_1d, 100, 0., 500);

}

void StopTreeLooper::makeCR5Plots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag, float mtcut ) 
{

  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_cr5_met"+tag_selection                   +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr5_met"+tag_selection+tag_njets         +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr5_met"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  //leading lepton pt - enters mT calculation
  plot1D("h_cr5_leppt"+tag_selection                   +flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_leppt"+tag_selection+tag_njets         +flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //isolated track pt
  plot1D("h_cr5_isotrkpt"+tag_selection                   +flav_tag, min(stopt.pfcand10().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_isotrkpt"+tag_selection+tag_njets         +flav_tag, min(stopt.pfcand10().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_isotrkpt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(stopt.pfcand10().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lead-lep and met
  float dphi_metlep = getdphi( stopt.lep1().Phi() , t1metphicorrphi );
  plot1D("h_cr5_dphi_metlep"+tag_selection                   +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_metlep"+tag_selection+tag_njets         +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_metlep"+tag_selection+tag_njets+tag_kbin+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //min dphi leading two jets - this should cut most of the ttbar single leptons
  plot1D("h_cr5_mindPhiJ12"+tag_selection                   +flav_tag, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_mindPhiJ12"+tag_selection+tag_njets         +flav_tag, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_mindPhiJ12"+tag_selection+tag_njets+tag_kbin+flav_tag, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  //b-pT
  plot1D("h_cr5_bpt"+tag_selection                   +flav_tag, pt_b, evtweight, h_1d, 50, 30., 400.);
  plot1D("h_cr5_bpt"+tag_selection+tag_njets         +flav_tag, pt_b, evtweight, h_1d, 50, 30., 400.);
  plot1D("h_cr5_bpt"+tag_selection+tag_njets+tag_kbin+flav_tag, pt_b, evtweight, h_1d, 50, 30., 400.);
  //maria vari5bles
  plot1D("h_cr5_htratiom"+tag_selection                   +flav_tag, htratiom, evtweight, h_1d, 50, 0., 1.);
  plot1D("h_cr5_htratiom"+tag_selection+tag_njets         +flav_tag, htratiom, evtweight, h_1d, 50, 0., 1.);
  plot1D("h_cr5_htratiom"+tag_selection+tag_njets+tag_kbin+flav_tag, htratiom, evtweight, h_1d, 50, 0., 1.);
  /// MT2 and chi2
  plot1D("h_cr5_mt2wmin"+tag_selection                   +flav_tag, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2wmin"+tag_selection+tag_njets         +flav_tag, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2wmin"+tag_selection+tag_njets+tag_kbin+flav_tag, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2bmin"+tag_selection                   +flav_tag, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2bmin"+tag_selection+tag_njets         +flav_tag, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2bmin"+tag_selection+tag_njets+tag_kbin+flav_tag, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2blmin"+tag_selection                   +flav_tag, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2blmin"+tag_selection+tag_njets         +flav_tag, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_mt2blmin"+tag_selection+tag_njets+tag_kbin+flav_tag, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr5_chi2minprob"+tag_selection                   +flav_tag, chi2minprob_ , evtweight, h_1d, 100, 0., 1);
  plot1D("h_cr5_chi2minprob"+tag_selection+tag_njets         +flav_tag, chi2minprob_ , evtweight, h_1d, 100, 0., 1);
  plot1D("h_cr5_chi2minprob"+tag_selection+tag_njets+tag_kbin+flav_tag, chi2minprob_ , evtweight, h_1d, 100, 0., 1);

  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_cr5_mt"+tag_selection                   +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt"+tag_selection+tag_njets         +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr5_mt_count"+tag_selection                   +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr5_mt_count"+tag_selection+tag_njets         +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr5_mt_count"+tag_selection+tag_njets+tag_kbin+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

  //Plot more angles between various objects
  //angle between lepton and isolated track
  float dphi_leptrk = getdphi( stopt.lep1().Phi() ,  stopt.pfcand10().Phi() );
  plot1D("h_cr5_dphi_leptrk"+tag_selection                   +flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_leptrk"+tag_selection+tag_njets         +flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr5_dphi_leptrk"+tag_selection+tag_njets+tag_kbin+flav_tag, dphi_leptrk, evtweight, h_1d, 15, 0., TMath::Pi());
  //dR between lepton and isolated track
  float dR_leptrk = dRbetweenVectors( stopt.lep1() , stopt.pfcand10() );
  plot1D("h_cr5_dR_leptrk"+tag_selection                   +flav_tag, min(dR_leptrk, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr5_dR_leptrk"+tag_selection+tag_njets         +flav_tag, min(dR_leptrk, (float)4.999), evtweight, h_1d, 15, 0., 5.);
  plot1D("h_cr5_dR_leptrk"+tag_selection+tag_njets+tag_kbin+flav_tag, min(dR_leptrk, (float)4.999), evtweight, h_1d, 15, 0., 5.);

}

void StopTreeLooper::makeSIGPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_kbin, string flav_tag, float mtcut ) 
{

  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;
  
  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  //default met
  plot1D("h_sig_met"+tag_selection         +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50, h_xmax);
  plot1D("h_sig_met"+tag_selection+tag_kbin+flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50, h_xmax);
  //lepton pt - enters mT calculation
  plot1D("h_sig_leppt"+tag_selection         +flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_leppt"+tag_selection+tag_kbin+flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //angle between lepton and met
  float dphi_metlep = getdphi( stopt.lep1().Phi() , t1metphicorrphi );
  plot1D("h_sig_dphi_metlep"+tag_selection         +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_sig_dphi_metlep"+tag_selection+tag_kbin+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //min dphi leading two jets - this should cut most of the ttbar single leptons
  plot1D("h_sig_mindPhiJ12"+tag_selection      +flav_tag, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  //b-pT
  plot1D("h_sig_bpt"+tag_selection             +flav_tag, pt_b, evtweight, h_1d, 50, 30., 400.);
  //maria variables
  plot1D("h_sig_htratiom"+tag_selection        +flav_tag, htratiom, evtweight, h_1d, 50, 0., 1.);
  /// MT2 and chi2
  plot1D("h_sig_mt2wmin"+tag_selection         +flav_tag, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_sig_mt2bmin"+tag_selection         +flav_tag, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_sig_mt2blmin"+tag_selection        +flav_tag, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_sig_chi2minprob"+tag_selection     +flav_tag, chi2minprob_ , evtweight, h_1d, 100, 0., 1);

  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_sig_mt"+tag_selection         +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mt"+tag_selection+tag_kbin+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_sig_mt_count"+tag_selection         +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_sig_mt_count"+tag_selection+tag_kbin+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

  if(dphimjmin>0.8) {

    plot1D("h_sig_mt_mindPhiJ12"+tag_selection         +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
    plot1D("h_sig_mt_mindPhiJ12"+tag_selection+tag_kbin+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
    plot1D("h_sig_mt_count_mindPhiJ12"+tag_selection         +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
    plot1D("h_sig_mt_count_mindPhiJ12"+tag_selection+tag_kbin+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
    
  }


}

void StopTreeLooper::makeCR1Plots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string tag_njets, string tag_kbin, string flav_tag, float mtcut ) 
{

  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;
  float x_ovflw = h_xmax-0.001;

  float mt_count = -1.;
  if ( t1metphicorrmt > min_mtpeak 
       && t1metphicorrmt < max_mtpeak )    mt_count = 0.5;
  else if ( t1metphicorrmt > mtcut ) mt_count = 1.5;
  
  plot1D("h_cr1_njets"    +tag_selection+flav_tag, min(n_jets,4),  evtweight, h_1d, 5,0,5);
  plot1D("h_cr1_njets_all"+tag_selection+flav_tag, min(n_jets,9),  evtweight, h_1d, 10, 0, 10);
  //default met
  plot1D("h_cr1_met"+tag_selection                   +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr1_met"+tag_selection+tag_njets         +flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  plot1D("h_cr1_met"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorr, x_ovflw), evtweight, h_1d, nbins-5, 50., h_xmax);
  //lepton pt - enters mT calculation
  plot1D("h_cr1_leppt"+tag_selection                   +flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_leppt"+tag_selection+tag_njets         +flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_leppt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(stopt.lep1().Pt(), x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //lepton phi
  plot1D("h_cr1_lepphi"+tag_selection                   +flav_tag, stopt.lep1().Phi(), evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_lepphi"+tag_selection+tag_njets         +flav_tag, stopt.lep1().Phi(), evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_lepphi"+tag_selection+tag_njets+tag_kbin+flav_tag, stopt.lep1().Phi(), evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  //met phi
  plot1D("h_cr1_metphi"+tag_selection                   +flav_tag, t1metphicorrphi, evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_metphi"+tag_selection+tag_njets         +flav_tag, t1metphicorrphi, evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  plot1D("h_cr1_metphi"+tag_selection+tag_njets+tag_kbin+flav_tag, t1metphicorrphi, evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  //angle between lepton and met
  float dphi_metlep = getdphi( stopt.lep1().Phi() , t1metphicorrphi );
  plot1D("h_cr1_dphi_metlep"+tag_selection                   +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metlep"+tag_selection+tag_njets         +flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr1_dphi_metlep"+tag_selection+tag_njets+tag_kbin+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., TMath::Pi());
  //min dphi leading two jets - this should cut most of the ttbar single leptons
  plot1D("h_cr1_mindPhiJ12"+tag_selection                   +flav_tag, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr1_mindPhiJ12"+tag_selection+tag_njets         +flav_tag, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  plot1D("h_cr1_mindPhiJ12"+tag_selection+tag_njets+tag_kbin+flav_tag, dphimjmin, evtweight, h_1d, 15, 0., TMath::Pi());
  //b-pT
  plot1D("h_cr1_bpt"+tag_selection                   +flav_tag, pt_b, evtweight, h_1d, 50, 30., 400.);
  plot1D("h_cr1_bpt"+tag_selection+tag_njets         +flav_tag, pt_b, evtweight, h_1d, 50, 30., 400.);
  plot1D("h_cr1_bpt"+tag_selection+tag_njets+tag_kbin+flav_tag, pt_b, evtweight, h_1d, 50, 30., 400.);
  //maria variables
  plot1D("h_cr1_htratiom"+tag_selection                   +flav_tag, htratiom, evtweight, h_1d, 50, 0., 1.);
  plot1D("h_cr1_htratiom"+tag_selection+tag_njets         +flav_tag, htratiom, evtweight, h_1d, 50, 0., 1.);
  plot1D("h_cr1_htratiom"+tag_selection+tag_njets+tag_kbin+flav_tag, htratiom, evtweight, h_1d, 50, 0., 1.);
  /// MT2 and chi2
  plot1D("h_cr1_mt2wmin"+tag_selection                   +flav_tag, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2wmin"+tag_selection+tag_njets         +flav_tag, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2wmin"+tag_selection+tag_njets+tag_kbin+flav_tag, mt2wmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2bmin"+tag_selection                   +flav_tag, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2bmin"+tag_selection+tag_njets         +flav_tag, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2bmin"+tag_selection+tag_njets+tag_kbin+flav_tag, mt2bmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2blmin"+tag_selection                   +flav_tag, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2blmin"+tag_selection+tag_njets         +flav_tag, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_mt2blmin"+tag_selection+tag_njets+tag_kbin+flav_tag, mt2blmin_ , evtweight, h_1d, 100, 0., 500);
  plot1D("h_cr1_chi2minprob"+tag_selection                   +flav_tag, chi2minprob_ , evtweight, h_1d, 100, 0., 1);
  plot1D("h_cr1_chi2minprob"+tag_selection+tag_njets         +flav_tag, chi2minprob_ , evtweight, h_1d, 100, 0., 1);
  plot1D("h_cr1_chi2minprob"+tag_selection+tag_njets+tag_kbin+flav_tag, chi2minprob_ , evtweight, h_1d, 100, 0., 1);

  //MT
  //binning for mT plots
  nbins = 30;
  h_xmin = 0.;
  h_xmax = 300.;
  x_ovflw = h_xmax-0.001;
  plot1D("h_cr1_mt"+tag_selection                   +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mt"+tag_selection+tag_njets         +flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mt"+tag_selection+tag_njets+tag_kbin+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_cr1_mt_count"+tag_selection                   +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr1_mt_count"+tag_selection+tag_njets         +flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);
  plot1D("h_cr1_mt_count"+tag_selection+tag_njets+tag_kbin+flav_tag, mt_count, evtweight, h_1d, 2, 0, 2);

}

void StopTreeLooper::makeNJPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				  string tag_selection, string flav_tag ) 
{

  plot1D("h_njets"    +tag_selection,          min(n_jets,4), evtweight, h_1d, 4,1,5);
  plot1D("h_njets"    +tag_selection+flav_tag, min(n_jets,4), evtweight, h_1d, 4,1,5);
  plot1D("h_njets_all"+tag_selection,          min(n_jets,8), evtweight, h_1d, 7,1,8);
  plot1D("h_njets_all"+tag_selection+flav_tag, min(n_jets,8), evtweight, h_1d, 7,1,8);

}

void StopTreeLooper::makeZPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				 string tag_selection, string tag_njets, string flav_tag ) 
{

  int nbins = 30;
  float h_xmin = 0.;
  float h_xmax = 300.;
  float x_ovflw = h_xmax-0.001;

  plot1D("h_z_met"+tag_selection+flav_tag, min(t1metphicorr,x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_z_met"+tag_selection+tag_njets+flav_tag, min(t1metphicorr,x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);

  string lep1type =  abs(stopt.id1())==13 ? "h_muo" : "h_ele";
  string lep2type =  abs(stopt.id2())==13 ? "h_muo" : "h_ele";
  plot1D(lep1type+"pt"+tag_selection+flav_tag, min(stopt.lep1().Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep2type+"pt"+tag_selection+flav_tag, min(stopt.lep2().Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep1type+"pt"+tag_selection+tag_njets+flav_tag, min(stopt.lep1().Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep2type+"pt"+tag_selection+tag_njets+flav_tag, min(stopt.lep2().Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  
  plot1D("h_z_leppt"  +tag_selection+flav_tag, min(stopt.lep1().Pt(),(float)299.99), evtweight, h_1d, 50, 20., 300.);
  plot1D("h_z_lepeta" +tag_selection+flav_tag, stopt.lep1().Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_lep2pt" +tag_selection+flav_tag, min(stopt.lep2().Pt(),(float)199.99), evtweight, h_1d, 50, 20., 200.);
  plot1D("h_z_lep2eta"+tag_selection+flav_tag, stopt.lep2().Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_leppt"  +tag_selection+tag_njets+flav_tag, min(stopt.lep1().Pt(),(float)299.99), evtweight, h_1d, 50, 20., 300.);
  plot1D("h_z_lepeta" +tag_selection+tag_njets+flav_tag, stopt.lep1().Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_lep2pt" +tag_selection+tag_njets+flav_tag, min(stopt.lep2().Pt(),(float)199.99), evtweight, h_1d, 50, 20., 200.);
  plot1D("h_z_lep2eta"+tag_selection+tag_njets+flav_tag, stopt.lep2().Eta(), evtweight, h_1d, 24, -2.4, 2.4);

  float dphi_metlep = getdphi(stopt.lep1().Phi(), t1metphicorrphi);
  plot1D("h_z_dphi_metl"+tag_selection+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);
  plot1D("h_z_dphi_metl"+tag_selection+tag_njets+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);
  plot1D("h_z_mt"+tag_selection+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);
  plot1D("h_z_mt"+tag_selection+tag_njets+flav_tag, min(t1metphicorrmt, x_ovflw), evtweight, h_1d, nbins, h_xmin, h_xmax);

  if ( n_jets<1 ) return;
  plot1D("h_z_j1pt" +tag_selection+flav_tag, min(jets.at(0).Pt(), (float)399.99), evtweight, h_1d, 20, 30., 400.);
  plot1D("h_z_j1eta"+tag_selection+flav_tag, jets.at(0).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( n_jets<2 ) return;
  plot1D("h_z_j2pt" +tag_selection+flav_tag, min(jets.at(1).Pt(), (float)299.99), evtweight, h_1d, 20, 30., 300.);
  plot1D("h_z_j2eta"+tag_selection+flav_tag, jets.at(1).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( n_jets<3 ) return;
  plot1D("h_z_j3pt" +tag_selection+flav_tag, min(jets.at(2).Pt(), (float)199.99), evtweight, h_1d, 20, 30., 200.);
  plot1D("h_z_j3eta"+tag_selection+flav_tag, jets.at(2).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  if ( n_jets<4 ) return;
  plot1D("h_z_j4pt" +tag_selection+flav_tag, min(jets.at(3).Pt(), (float)119.99), evtweight, h_1d, 20, 30., 120.);
  plot1D("h_z_j4eta"+tag_selection+flav_tag, jets.at(3).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  
}


