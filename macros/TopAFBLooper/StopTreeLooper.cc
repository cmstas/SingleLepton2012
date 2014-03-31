#include "StopTreeLooper.h"

//#include "../../CORE/jetSmearingTools.h"
//#include "../../CORE/Thrust.h"
//#include "../../CORE/EventShape.h"
#include "TStopwatch.h"

#include "Math/VectorUtil.h"
#include "../Core/STOPT.h"
#include "../Core/stopUtils.h"
#include "../Plotting/PlotUtilities.h"
#include "../../Tools/BTagReshaping/BTagReshaping.h"

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
  min_njets = -9999;
  t1metphicorr = -9999.;
  t1metphicorrphi = -9999.;
  n_jets  = -9999;
  n_bjets = -9999;
  n_ljets = -9999;
  pfcalo_metratio = -9999.;
  pfcalo_metdphi  = -9999.;
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

  //------------------------------------------------------------------------------------------------------
  // set csv discriminator reshaping
  //------------------------------------------------------------------------------------------------------
  
  BTagShapeInterface * nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0);

  //------------------------------
  // set up histograms
  //------------------------------

  gROOT->cd();

  cout << "[StopTreeLooper::loop] setting up histos" << endl;

  std::map<std::string, TH1F*> h_1d;
  //also for control regions
  std::map<std::string, TH1F*> h_1d_cr1, h_1d_cr2, h_1d_cr3, h_1d_cr4;
  //for signal region 
  std::map<std::string, TH1F*> h_1d_sig;
  //for ttbar dilepton njets distribution
  std::map<std::string, TH1F*> h_1d_nj;
  //z sample for yields etc
  std::map<std::string, TH1F*> h_1d_z;

  //-----------------------------------
  // PU reweighting based on true PU
  //-----------------------------------

  TFile* pu_file = TFile::Open("../vtxreweight/puWeights_Summer12_53x_True_19p5ifb.root");
  if( pu_file == 0 ){
    cout << "vtxreweight error, couldn't open vtx file. Quitting!"<< endl;
    exit(0);
  }

  TH1F* h_pu_wgt = (TH1F*)pu_file->Get("puWeights");
  h_pu_wgt->SetName("h_pu_wgt");

  //------------------------------
  // file loop
  //------------------------------

  unsigned int nEventsChain=0;
  unsigned int nEvents = chain->GetEntries();
  nEventsChain = nEvents;
  ULong64_t nEventsTotal = 0;
  int nevt_check = 0;

  bool isData = name.Contains("data") ? true : false;

  //Define jet multiplicity requirement
  min_njets = 2;
  printf("[StopTreeLooper::loop] N JET min. requirement for signal %i \n", min_njets);

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

      nevt_check++;
      //---------------------------------------------------------------------------- 
      // determine event weight
      // make 2 example histograms of nvtx and corresponding weight
      //---------------------------------------------------------------------------- 

      // to reweight from the nvtx distribution
      // float evtweight = isData ? 1. : 
      // 	( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );
      float puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt, isData );
      float evtweight = isData ? 1. : 
	( stopt.weight() * 19.5 * puweight );
      if (!name.Contains("lmg")) evtweight *= stopt.mgcor();

      plot1D("h_vtx",       stopt.nvtx(), evtweight, h_1d, 40, 0, 40);
      plot1D("h_vtxweight",     puweight, evtweight, h_1d, 41, -4., 4.);

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

      pfcalo_metratio = t1metphicorr/stopt.calomet();
      pfcalo_metdphi  = getdphi(t1metphicorrphi, stopt.calometphi());

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

      for( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i ){
	
	if( stopt.pfjets().at(i).pt()<30 )  continue;
	if( fabs(stopt.pfjets().at(i).eta())>2.4 )  continue;  //note, this was 2.5 for 7 TeV AFB
	//	if( stopt.pfjets_beta2_0p5().at(i)<0.2 )  continue;

	//	bool passMediumPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mvaPUid().at(i),1);
	bool passTightPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mva5xPUid().at(i),0);

	if(!passTightPUid) continue;

	n_jets++;
	n_ljets++;

  jets.push_back( stopt.pfjets().at(i) );

	//to not use reshaped discriminator
	//float csv_nominal= stopt.pfjets_csv().at(i);
	float csv_nominal=isData ? stopt.pfjets_csv().at(i)
	  : nominalShape->reshape( stopt.pfjets().at(i).eta(),
				   stopt.pfjets().at(i).pt(),
				   stopt.pfjets_csv().at(i),
				   stopt.pfjets_mcflavorAlgo().at(i) ); 
	if (csv_nominal > 0.679) {
	  n_bjets++;
	}
	btag.push_back( csv_nominal );

      	if ( !isData ) mc.push_back  ( stopt.pfjets_mc3().at(i) );
      	else mc.push_back  ( 0 );

      	sigma_jets.push_back(stopt.pfjets_sigma().at(i));
 
        //count jets that are not overlapping with second lepton
	if (isData) continue;
	if (stopt.nleps()<2) continue;
	if (stopt.mclep2().pt() < 30.) continue;
	if (ROOT::Math::VectorUtil::DeltaR(stopt.mclep2(), stopt.pfjets().at(i)) > 0.4 ) continue;
	n_ljets--;

      } 

      
      //----------------------------------------------------------------------------
      // histogram tags
      //----------------------------------------------------------------------------

      //b-tagging
      string tag_btag = (n_bjets<1) ? "_bveto" : "";
 
      //iso-trk-veto & tau veto
      bool passisotrk = passIsoTrkVeto_v4() && passTauVeto();
      string tag_isotrk = passisotrk ? "" : "_wisotrk";
      //string tag_isotrk = (passLepPlusIsoTrkSelection(isData)) ? "" : "_wisotrk";   

      //z-peak/veto
      string tag_zcut;
      if ( fabs( stopt.dilmass() - 91.) > 15. ) tag_zcut = "_zveto";
      else if  ( fabs( stopt.dilmass() - 91.) < 10. ) tag_zcut = "_zpeak";
      else tag_zcut = "_ignore";

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
      if ( abs(stopt.id1()) != abs(stopt.id2()) && flav_tag_dl != "_mysterydl" ) basic_flav_tag_dl = "_mueg";



      //------------------------------------------ 
      // datasets bit
      //------------------------------------------ 
      
      bool dataset_1l=false;

      if((isData) && name.Contains("muo") 
	 && (abs(stopt.id1()) == 13 ))  dataset_1l=true;
      if((isData) && name.Contains("ele") 
	 && (abs(stopt.id1()) == 11 ))  dataset_1l=true;

      if(!isData) dataset_1l=true;

      bool dataset_2l=false;

      if((isData) && name.Contains("dimu") 
	 && (abs(stopt.id1()) == 13 ) 
	 && (abs(stopt.id2())==13)) dataset_2l=true;
      if((isData) && name.Contains("diel") 
	 && (abs(stopt.id1()) == 11 ) 
	 && (abs(stopt.id2())==11)) dataset_2l=true;
      if((isData) && name.Contains("mueg") 
	 && abs(stopt.id1()) != abs(stopt.id2())) 
	dataset_2l=true;

      if(!isData) dataset_2l=true;

      float trigweight = isData ? 1. : getsltrigweight(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta());
      float trigweight_dl = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());

      //
      // SIGNAL REGION - dilepton + b-tag
      //
      // selection - 1 lepton 
      // Add iso track veto
      // Add b-tag

      //need proper signal region cuts here, including 3rd lepton veto etc.

      if ( dataset_2l && passDileptonSelection(isData) 
     && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. ) 
     && n_bjets>0 ) 
  {
      makeNJPlots( evtweight*trigweight_dl, h_1d_nj, "", basic_flav_tag_dl);
      if ( n_jets < min_njets  ) makeSIGPlots( evtweight*trigweight_dl, h_1d_sig, tag_isotrk+tag_btag  , basic_flav_tag_dl );
  }

/*
      //
      // CR1 - single lepton + b-veto
      //

      // selection - 1 lepton + iso track veto
      // Add b-tag veto
      if ( dataset_1l && passSingleLeptonSelection(isData) 
	   && passisotrk
	   && n_jets>=min_njets )
	{
	    //pre b-tag veto
	    makeCR1Plots( evtweight*trigweight, h_1d_cr1, "_prebveto", flav_tag_sl );
	    //b-veto
	    if ( n_bjets==0 ) makeCR1Plots( evtweight*trigweight, h_1d_cr1,, flav_tag_sl );

	}//end CR1 selection
      
      //
      // CR2 - Z-peak for yields and mT resolution studies
      // 
      
      // selection - SF dilepton, veto on isolated track in addition to 2 leptons, in z-peak
      if ( dataset_2l && passDileptonSelection(isData) )
	{

	  //invariant mass - basic check of inclusive distribution
	  plot1D("h_z_dilmass"+flav_tag_dl, stopt.dilmass(), evtweight*trigweight_dl, h_1d_z,  30 , 76 , 106);

	  if ( fabs( stopt.dilmass() - 91.) < 10. ) 
	    {

	      // if (n_jets>8) 
	      // 	cout<<"NJETS: "<<n_jets<<" * dataset: "<<stopt.dataset()
	      // 	    <<" run: "<<stopt.run()<<" lumi: "<<stopt.lumi()<<" event: "<<stopt.event()<<endl;
	      
	      //z peak plots
	      plot1D("h_z_njets"    +flav_tag_dl, min(n_jets,4),  evtweight*trigweight_dl, h_1d_z, 5,0,5);
	      plot1D("h_z_njets_all"+flav_tag_dl, min(n_jets,9),  evtweight*trigweight_dl, h_1d_z, 10, 0, 10);
	      plot1D("h_z_nbjets"   +flav_tag_dl, min(n_bjets,3), evtweight*trigweight_dl, h_1d_z, 4, 0, 4);
	      makeZPlots( evtweight*trigweight_dl, h_1d_z, "", flav_tag_dl );

	      // Add stricter 3rd lepton veto
	      // require at least 2 jets
	      // Add b-tag veto 
	      if ( (stopt.trkpt10loose() <0.0001 || stopt.trkreliso10loose() > 0.1) 
		   && n_jets>=min_njets 
		   && n_bjets==0 ) {

          makeCR2Plots( evtweight*trigweight_dl, h_1d_cr2,, basic_flav_tag_dl );

	      }
	    }
	}//end CR2 selection



      //
      // CR3 - single lepton control REGION - single lepton + b-tag
      //
      // selection - 1 lepton 
      // Add iso track veto
      // Add b-tag

      if ( dataset_1l && passSingleLeptonSelection(isData) 
     && n_jets>=min_njets )
  {
      makeCR3Plots( evtweight*trigweight, h_1d_cr3, tag_isotrk+tag_btag  , flav_tag_sl );
  }
*/

    } // end event loop
    
   // delete tree;
    
    stwatch.Stop();
    cout<<"time is "<<stwatch.CpuTime()<<endl;
    stwatch.Start();
    
  } // end file loop
  
    //
    // finish
    //
  
  cout<<"N EVENT CHECK "<<nevt_check<<endl;
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
/*
  //control regions
  //h_1d_cr1, h_1d_cr2, h_1d_cr4, h_1d_cr3

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

  TFile outfile_cr3(Form("CR3%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving CR3 histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_cr3;
  for(it1d_cr3=h_1d_cr3.begin(); it1d_cr3!=h_1d_cr3.end(); it1d_cr3++) {
    it1d_cr3->second->Write(); 
    delete it1d_cr3->second;
  }

  outfile_cr3.Write();
  outfile_cr3.Close();
*/
  TFile outfile_nj(Form("NJ%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving NJ histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_nj;
  for(it1d_nj=h_1d_nj.begin(); it1d_nj!=h_1d_nj.end(); it1d_nj++) {
    it1d_nj->second->Write(); 
    delete it1d_nj->second;
  }

  outfile_nj.Write();
  outfile_nj.Close();
/*
  TFile outfile_z(Form("Z%s",m_outfilename_.c_str()),"RECREATE") ; 
  printf("[StopTreeLooper::loop] Saving Z histograms to %s\n", m_outfilename_.c_str());

  std::map<std::string, TH1F*>::iterator it1d_z;
  for(it1d_z=h_1d_z.begin(); it1d_z!=h_1d_z.end(); it1d_z++) {
    it1d_z->second->Write(); 
    delete it1d_z->second;
  }

  outfile_z.Write();
  outfile_z.Close();
*/
  already_seen.clear();

  gROOT->cd();

}



void StopTreeLooper::makeSIGPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
				   string tag_selection, string flav_tag) 
{

  int nbins = 50;

  //default met
  plot1D("h_sig_met"+tag_selection+flav_tag, t1metphicorr, evtweight, h_1d, nbins-5, 50, 500);

  //check HO and TOB/TEC cleanup cut variables
  plot1D("h_sig_pfcaloMET"+tag_selection+flav_tag, min(pfcalo_metratio, (float)3.9999) , evtweight, h_1d, 100, 0, 4.);
  plot1D("h_sig_pfcalodPhi"+tag_selection+flav_tag, pfcalo_metdphi , evtweight, h_1d, 100, 0, TMath::Pi());

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
				 string tag_selection, string flav_tag ) 
{

  int nbins = 30;
  float h_xmin = 0.;
  float h_xmax = 300.;

  plot1D("h_z_met"+tag_selection+flav_tag, t1metphicorr, evtweight, h_1d, nbins, h_xmin, h_xmax);

  string lep1type =  abs(stopt.id1())==13 ? "h_muo" : "h_ele";
  string lep2type =  abs(stopt.id2())==13 ? "h_muo" : "h_ele";
  plot1D(lep1type+"pt"+tag_selection+flav_tag, min(stopt.lep1().Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  plot1D(lep2type+"pt"+tag_selection+flav_tag, min(stopt.lep2().Pt(),(float)199.99), evtweight, h_1d, 40, 20, 200);
  
  plot1D("h_z_leppt"  +tag_selection+flav_tag, min(stopt.lep1().Pt(),(float)299.99), evtweight, h_1d, 50, 20., 300.);
  plot1D("h_z_lepeta" +tag_selection+flav_tag, stopt.lep1().Eta(), evtweight, h_1d, 24, -2.4, 2.4);
  plot1D("h_z_lep2pt" +tag_selection+flav_tag, min(stopt.lep2().Pt(),(float)199.99), evtweight, h_1d, 50, 20., 200.);
  plot1D("h_z_lep2eta"+tag_selection+flav_tag, stopt.lep2().Eta(), evtweight, h_1d, 24, -2.4, 2.4);

  float dphi_metlep = getdphi(stopt.lep1().Phi(), t1metphicorrphi);
  plot1D("h_z_dphi_metl"+tag_selection+flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);

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


  
