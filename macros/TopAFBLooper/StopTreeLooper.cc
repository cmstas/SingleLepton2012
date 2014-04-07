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
#include "TPython.h"

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
    if (listOfFiles->GetEntries() == 0)
    {
        cout << "[StopTreeLooper::loop] no files in chain" << endl;
        return;
    }

    //------------------------------------------------------------------------------------------------------
    // load Betchart solver and AMWT
    //------------------------------------------------------------------------------------------------------

    TPython::LoadMacro("loadBetchart.py");
    //d_llsol = new ttdilepsolve;

    //------------------------------------------------------------------------------------------------------
    // set csv discriminator reshaping
    //------------------------------------------------------------------------------------------------------

    BTagShapeInterface *nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0);

    //------------------------------
    // set up histograms
    //------------------------------

    gROOT->cd();

    cout << "[StopTreeLooper::loop] setting up histos" << endl;

    std::map<std::string, TH1F *> h_1d;
    //also for control regions
    std::map<std::string, TH1F *> h_1d_cr1, h_1d_cr2, h_1d_cr3, h_1d_cr4;
    //for signal region
    std::map<std::string, TH1F *> h_1d_sig;
    //for ttbar dilepton njets distribution
    std::map<std::string, TH1F *> h_1d_nj;
    //z sample for yields etc
    std::map<std::string, TH1F *> h_1d_z;

    //-----------------------------------
    // PU reweighting based on true PU
    //-----------------------------------

    TFile *pu_file = TFile::Open("../vtxreweight/puWeights_Summer12_53x_True_19p5ifb.root");
    if ( pu_file == 0 )
    {
        cout << "vtxreweight error, couldn't open vtx file. Quitting!" << endl;
        exit(0);
    }

    TH1F *h_pu_wgt = (TH1F *)pu_file->Get("puWeights");
    h_pu_wgt->SetName("h_pu_wgt");

    //------------------------------
    // file loop
    //------------------------------

    unsigned int nEventsChain = 0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    ULong64_t nEventsTotal = 0;
    int nevt_check = 0;
    int nevt_nlep[10] = {0};

    bool isData = name.Contains("data") ? true : false;

    //Define jet multiplicity requirement
    min_njets = 2;
    printf("[StopTreeLooper::loop] N JET min. requirement for signal %i \n", min_njets);

    cout << "[StopTreeLooper::loop] running over chain with total entries " << nEvents << endl;

    stwatch.Start();

    while (TChainElement *currentFile = (TChainElement *)fileIter.Next())
    {

        //----------------------------
        // load the stop baby tree
        //----------------------------

        TFile *file = new TFile( currentFile->GetTitle() );
        TTree *tree = (TTree *)file->Get("t");
        stopt.Init(tree);

        //----------------------------
        // event loop
        //----------------------------

        ULong64_t nEvents = tree->GetEntriesFast();
        for (ULong64_t event = 0; event < nEvents; ++event)
        {
            stopt.GetEntry(event);

            //----------------------------
            // increment counters
            //----------------------------

            ++nEventsTotal;
            if (nEventsTotal % 10000 == 0)
            {
                ULong64_t i_permille = (int)floor(1000 * nEventsTotal / float(nEventsChain));
                //if (i_permille != i_permille_old) {//this prints too often!
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1))
                {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                           "\033[0m\033[32m <---\033[0m\015", i_permille / 10.);
                    fflush(stdout);
                    stwatch.Stop();
                    if (i_permille % 100 < 0.0001)
                        cout << "At " << i_permille / 10. << "% time is " << stwatch.RealTime() << " cpu: " << stwatch.CpuTime() << endl;
                    stwatch.Start();//stwatch.Continue();

                }
                //  i_permille_old = i_permille;
            }

            //---------------------
            // skip duplicates
            //---------------------

            if ( isData )
            {
                DorkyEventIdentifier id = {stopt.run(), stopt.event(), stopt.lumi() };
                if (is_duplicate(id, already_seen) )
                {
                    continue;
                }
                if (is_badLaserEvent(id, events_lasercalib) )
                {
                    //std::cout<< "Removed bad laser calibration event:" << run << "   " << event<<"\n";
                    continue;
                }
                if (is_badLaserEvent(id, events_hcallasercalib) )
                {
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
            //    ( stopt.weight() * 19.5 * stopt.nvtxweight() * stopt.mgcor() );
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
            nevt_nlep[stopt.ngoodlep()]++;
            if (stopt.ngoodlep() > 2) continue; //veto the third lepton

            //----------------------------------------------------------------------------
            // Function to perform MET phi corrections on-the-fly
            // Default branches are: tree->t1metphicorr_ and tree->t1metphicorrmt_
            //----------------------------------------------------------------------------

            pair<float, float> p_t1metphicorr =
                getPhiCorrMET( stopt.t1met10(), stopt.t1met10phi(), stopt.nvtx(), !isData);
            t1metphicorr    = p_t1metphicorr.first;
            t1metphicorrphi = p_t1metphicorr.second;

            pfcalo_metratio = t1metphicorr / stopt.calomet();
            pfcalo_metdphi  = getdphi(t1metphicorrphi, stopt.calometphi());

            //----------------------------------------------------------------------------
            // get jet information
            //----------------------------------------------------------------------------

            jets.clear();
            bjets.clear();
            nonbjets.clear();
            bcandidates.clear();
            btag.clear();
            sigma_jets.clear();
            mc.clear();

            n_jets  = 0;
            n_bjets = 0;
            n_ljets = 0;

            //stopt.pfjets() are already sorted by pT
            for ( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i )
            {

                if ( stopt.pfjets().at(i).pt() < 30 )  continue;
                if ( fabs(stopt.pfjets().at(i).eta()) > 2.4 )  continue; //note, this was 2.5 for 7 TeV AFB
                //  if( stopt.pfjets_beta2_0p5().at(i)<0.2 )  continue;

                //  bool passMediumPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mvaPUid().at(i),1);
                bool passTightPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(), stopt.pfjets_mva5xPUid().at(i), 0);

                if (!passTightPUid) continue;

                n_jets++;
                n_ljets++;

                jets.push_back( stopt.pfjets().at(i) );

                //to not use reshaped discriminator
                //float csv_nominal= stopt.pfjets_csv().at(i);
                float csv_nominal = isData ? stopt.pfjets_csv().at(i)
                                    : nominalShape->reshape( stopt.pfjets().at(i).eta(),
                                            stopt.pfjets().at(i).pt(),
                                            stopt.pfjets_csv().at(i),
                                            stopt.pfjets_mcflavorAlgo().at(i) );
                if (csv_nominal > 0.679)
                {
                    bjets.push_back( stopt.pfjets().at(i) );
                    //if more than 2 b-jets, use the ones with highest pT as b candidates
                    if (bcandidates.size()<2) bcandidates.push_back( stopt.pfjets().at(i) );
                    n_bjets++;
                }
                else
                {
                    nonbjets.push_back( stopt.pfjets().at(i) );
                }
                btag.push_back( csv_nominal );

                if ( !isData ) mc.push_back  ( stopt.pfjets_mc3().at(i) );
                else mc.push_back  ( 0 );

                sigma_jets.push_back(stopt.pfjets_sigma().at(i));

                //count jets that are not overlapping with second lepton
                if (isData) continue;
                if (stopt.nleps() < 2) continue;
                if (stopt.mclep2().pt() < 30.) continue;
                if (ROOT::Math::VectorUtil::DeltaR(stopt.mclep2(), stopt.pfjets().at(i)) > 0.4 ) continue;
                n_ljets--;

            }

            //if <2 btagged jets, take the highest pT light jets as b candidates
            if (bcandidates.size()<2 && nonbjets.size()>0) {
                bcandidates.push_back( nonbjets.at(0) );
            }
            if (bcandidates.size()<2 && nonbjets.size()>1) {
                bcandidates.push_back( nonbjets.at(1) );
            }

            //require at least 2 jets and leptons
            if(!(n_jets > 1 && stopt.id2() != -999)) continue;

            //----------------------------------------------------------------------------
            // ttbar solver
            //----------------------------------------------------------------------------  

            nu1_vecs.clear();
            nu2_vecs.clear();
            top1_vecs.clear();
            top2_vecs.clear();
            AMWT_weights.clear();

            //lepPlus.SetPtEtaPhiE(0,0,0,0);
            //lepMinus.SetPtEtaPhiE(0,0,0,0);
            //jet1.SetPtEtaPhiE(0,0,0,0);
            //jet2.SetPtEtaPhiE(0,0,0,0);

            if(stopt.id1() > 0) {
                lepPlus.SetPtEtaPhiE(stopt.lep1().Pt(),stopt.lep1().Eta(),stopt.lep1().Phi(),stopt.lep1().E());
                lepMinus.SetPtEtaPhiE(stopt.lep2().Pt(),stopt.lep2().Eta(),stopt.lep2().Phi(),stopt.lep2().E());
            }
            else {
                lepMinus.SetPtEtaPhiE(stopt.lep1().Pt(),stopt.lep1().Eta(),stopt.lep1().Phi(),stopt.lep1().E());
                lepPlus.SetPtEtaPhiE(stopt.lep2().Pt(),stopt.lep2().Eta(),stopt.lep2().Phi(),stopt.lep2().E());
            }
            jet1.SetPtEtaPhiE(bcandidates.at(0).Pt(),bcandidates.at(0).Eta(),bcandidates.at(0).Phi(),bcandidates.at(0).E());
            jet2.SetPtEtaPhiE(bcandidates.at(1).Pt(),bcandidates.at(1).Eta(),bcandidates.at(1).Phi(),bcandidates.at(1).E());

            top1_p4.SetPtEtaPhiE(0,0,0,0);
            top2_p4.SetPtEtaPhiE(0,0,0,0);
            nusum.SetPtEtaPhiE(0,0,0,0);
            cms.SetPtEtaPhiE(0,0,0,0);

            m_top = -999;
            m_top_B = -999;
            closestDeltaMET_maxwcombo = -999;
            closestDeltaMET_othercombo = -999;
            closestDeltaMET_bestcombo = -999;
            imaxweight = -1;
            closestApproach = false;

            solvettbar();

            m_top = m_top_B;

            //----------------------------------------------------------------------------
            // asymmetry calculations
            //----------------------------------------------------------------------------

            lep_charge_asymmetry = -999.0;
            lep_azimuthal_asymmetry = -999.0;
            lep_azimuthal_asymmetry_2 = -999.0;
            top_rapiditydiff_cms = -999.0;
            top_pseudorapiditydiff_cms = -999.0;
            top_rapiditydiff_Marco = -999.0;
            top_costheta_cms = -999.0;
            lepPlus_costheta_cms = -999.0;
            lepMinus_costheta_cms = -999.0;
            top_spin_correlation = -999.0;
            lep_cos_opening_angle = -999.0;
            tt_mass = -999.0;
            ttRapidity2 = -999.0;
            tt_pT = -999.0;
            top1_pt = -999.0;
            top2_pt = -999.0;
            top1_p_CM = -999.0;
            top2_p_CM = -999.0;
            top_rapiditydiffsigned_cms = -999.0;

            //fully leptonic asymmetries
            lep_charge_asymmetry = abs(lepPlus.Eta()) - abs(lepMinus.Eta());
            lep_azimuthal_asymmetry = lepPlus.DeltaPhi(lepMinus);  //lep_azimuthal_asymmetry is same as lep_azimuthal_asymmetry_2 but from -pi to pi instead of folding it over from 0 to pi
            lep_azimuthal_asymmetry_2 = acos(cos(lepPlus.DeltaPhi(lepMinus)));

            //more variables for plots
            float lep_pseudorap_diff = (lepPlus.Eta()) - (lepMinus.Eta());
            float lep_cosalpha =  lepPlus.Vect().Dot( lepMinus.Vect() ) / (lepPlus.Vect().Mag() * lepMinus.Vect().Mag());
            float lepPlus_phi = lepPlus.Phi();
            float lepMinus_phi = lepMinus.Phi();
            float lepPlus_Eta = lepPlus.Eta();
            float lepMinus_Eta = lepMinus.Eta();
            float lepPlus_Pt = lepPlus.Pt();
            float lepMinus_Pt = lepMinus.Pt();

            float jet_azimuthal_asymmetry = jet1.DeltaPhi(jet2);
            float jet_pseudorap_diff = jet1.Eta() - jet2.Eta();
            float jet_cosalpha =  jet1.Vect().Dot( jet2.Vect() ) / (jet1.Vect().Mag() * jet2.Vect().Mag());
            float jet1_phi = jet1.Phi();
            float jet2_phi = jet2.Phi();

            //variables that require ttbar solution
            if ( m_top > 0 ) {

                tt_mass = (top1_p4 + top2_p4).M();
                ttRapidity2 = (top1_p4 + top2_p4).Rapidity();
                top_rapiditydiff_cms = (top1_p4.Rapidity() - top2_p4.Rapidity()) * (top1_p4.Rapidity() + top2_p4.Rapidity());
                top_pseudorapiditydiff_cms = abs(top1_p4.Eta()) - abs(top2_p4.Eta());
                top_rapiditydiff_Marco = abs(top1_p4.Rapidity()) - abs(top2_p4.Rapidity());
                top_rapiditydiffsigned_cms = (top1_p4.Rapidity() - top2_p4.Rapidity());

                top1_pt =  top1_p4.Pt();
                top2_pt =  top2_p4.Pt();
                cms = top1_p4 + top2_p4;
                tt_pT = cms.Pt();

                top1_p4.Boost(-cms.BoostVector());
                top2_p4.Boost(-cms.BoostVector());

                top1_p_CM =  top1_p4.P();
                top2_p_CM =  top2_p4.P();

                top_costheta_cms = top1_p4.Vect().Dot(cms.Vect()) / (top1_p4.Vect().Mag() * cms.Vect().Mag());

                lepPlus.Boost(-cms.BoostVector());
                lepMinus.Boost(-cms.BoostVector());

                lepPlus.Boost(-top1_p4.BoostVector());
                lepMinus.Boost(-top2_p4.BoostVector());

                lepPlus_costheta_cms = lepPlus.Vect().Dot(top1_p4.Vect()) / (lepPlus.Vect().Mag() * top1_p4.Vect().Mag());
                lepMinus_costheta_cms = lepMinus.Vect().Dot(top2_p4.Vect()) / (lepMinus.Vect().Mag() * top2_p4.Vect().Mag());
                top_spin_correlation = lepPlus_costheta_cms * lepMinus_costheta_cms;
                lep_cos_opening_angle = lepPlus.Vect().Dot(lepMinus.Vect()) / (lepPlus.Vect().Mag() * lepMinus.Vect().Mag());

            }
            

            //----------------------------------------------------------------------------
            // histogram tags
            //----------------------------------------------------------------------------

            //b-tagging
            string tag_btag = (n_bjets < 1) ? "_bveto" : "";

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
            if ( abs(stopt.id1()) == 13 ) flav_tag_sl = "_muo";
            else if ( abs(stopt.id1()) == 11 ) flav_tag_sl = "_ele";
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
            if (stopt.id1()*stopt.id2() > 0) basic_flav_tag_dl += "_SS";



            //------------------------------------------
            // datasets bit
            //------------------------------------------

            bool dataset_1l = false;

            if ((isData) && name.Contains("muo")
                    && (abs(stopt.id1()) == 13 ))  dataset_1l = true;
            if ((isData) && name.Contains("ele")
                    && (abs(stopt.id1()) == 11 ))  dataset_1l = true;

            if (!isData) dataset_1l = true;

            bool dataset_2l = false;

            if ((isData) && name.Contains("dimu")
                    && (abs(stopt.id1()) == 13 )
                    && (abs(stopt.id2()) == 13)) dataset_2l = true;
            if ((isData) && name.Contains("diel")
                    && (abs(stopt.id1()) == 11 )
                    && (abs(stopt.id2()) == 11)) dataset_2l = true;
            if ((isData) && name.Contains("mueg")
                    && abs(stopt.id1()) != abs(stopt.id2()))
                dataset_2l = true;

            if (!isData) dataset_2l = true;

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
                    && n_bjets > 0
                    && (stopt.lep1() + stopt.lep2()).M() >= 20.0 
                    && (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. ) )
            {
                //check cuts
                if ( (stopt.lep1().Pt() < 20.0 || stopt.lep2().Pt() < 20.0 ) ) cout<<stopt.lep1().Pt()<<" "<<stopt.lep2().Pt()<<endl;
                if ( (stopt.lep1() + stopt.lep2()).M() < 20.0 )  cout<<(stopt.lep1() + stopt.lep2()).M()<<endl;
                if ( t1metphicorr < 40. && basic_flav_tag_dl != "_mueg" )  cout<<"MET: "<<t1metphicorr<<" "+basic_flav_tag_dl<<endl;

                makeNJPlots( evtweight * trigweight_dl, h_1d_nj, "", basic_flav_tag_dl);
                if ( n_jets >= min_njets  ) makeSIGPlots( evtweight * trigweight_dl, h_1d_sig,  tag_btag  , basic_flav_tag_dl );
                if ( n_jets >= min_njets  ) makeSIGPlots( evtweight * trigweight_dl, h_1d_sig,  tag_btag  , "_all" );
                //if ( n_jets >= min_njets  ) makettPlots( evtweight * trigweight_dl, h_1d_sig,  tag_btag  , basic_flav_tag_dl );
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
                      //    cout<<"NJETS: "<<n_jets<<" * dataset: "<<stopt.dataset()
                      //        <<" run: "<<stopt.run()<<" lumi: "<<stopt.lumi()<<" event: "<<stopt.event()<<endl;

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
        cout << "time is " << stwatch.CpuTime() << endl;
        stwatch.Start();

    } // end file loop

    //
    // finish
    //

    cout << "N EVENT CHECK " << nevt_check << endl;
    TFile outfile(m_outfilename_.c_str(), "RECREATE") ;
    printf("[StopTreeLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());

    std::map<std::string, TH1F *>::iterator it1d;
    for (it1d = h_1d.begin(); it1d != h_1d.end(); it1d++)
    {
        it1d->second->Write();
        delete it1d->second;
    }

    outfile.Write();
    outfile.Close();

    TFile outfile_sig(Form("SIG%s", m_outfilename_.c_str()), "RECREATE") ;
    printf("[StopTreeLooper::loop] Saving SIG histograms to %s\n", m_outfilename_.c_str());

    std::map<std::string, TH1F *>::iterator it1d_sig;
    for (it1d_sig = h_1d_sig.begin(); it1d_sig != h_1d_sig.end(); it1d_sig++)
    {
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
    TFile outfile_nj(Form("NJ%s", m_outfilename_.c_str()), "RECREATE") ;
    printf("[StopTreeLooper::loop] Saving NJ histograms to %s\n", m_outfilename_.c_str());

    printf("[StopTreeLooper::loop] nevt_nlep %i %i %i %i %i %i %i\n",nevt_nlep[0],nevt_nlep[1],nevt_nlep[2],nevt_nlep[3],nevt_nlep[4],nevt_nlep[5],nevt_nlep[6]) ;

    std::map<std::string, TH1F *>::iterator it1d_nj;
    for (it1d_nj = h_1d_nj.begin(); it1d_nj != h_1d_nj.end(); it1d_nj++)
    {
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



void StopTreeLooper::makeSIGPlots( float evtweight, std::map<std::string, TH1F *> &h_1d,
                                   string tag_selection, string flav_tag)
{

    int nbins = 80;

    plot1DUnderOverFlow("h_sig_lep_charge_asymmetry" + tag_selection + flav_tag, lep_charge_asymmetry , evtweight, h_1d, nbins, -4, 4);    
    plot1DUnderOverFlow("h_sig_lep_azimuthal_asymmetry" + tag_selection + flav_tag, lep_azimuthal_asymmetry , evtweight, h_1d, nbins, -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow("h_sig_lep_azimuthal_asymmetry_2" + tag_selection + flav_tag, lep_azimuthal_asymmetry_2 , evtweight, h_1d, nbins, 0, TMath::Pi());
    plot1DUnderOverFlow("h_sig_m_top" + tag_selection + flav_tag, m_top , evtweight, h_1d, nbins-1, 171.4, 173.6);

    if(m_top > 0) {
        plot1DUnderOverFlow("h_sig_top_rapiditydiff_cms" + tag_selection + flav_tag, top_rapiditydiff_cms , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow("h_sig_top_pseudorapiditydiff_cms" + tag_selection + flav_tag, top_pseudorapiditydiff_cms , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow("h_sig_top_rapiditydiff_Marco" + tag_selection + flav_tag, top_rapiditydiff_Marco , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow("h_sig_top_costheta_cms" + tag_selection + flav_tag, top_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_sig_lepPlus_costheta_cms" + tag_selection + flav_tag, lepPlus_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_sig_lepMinus_costheta_cms" + tag_selection + flav_tag, lepMinus_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_sig_top_spin_correlation" + tag_selection + flav_tag, top_spin_correlation , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_sig_lep_cos_opening_angle" + tag_selection + flav_tag, lep_cos_opening_angle , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_sig_tt_mass" + tag_selection + flav_tag, tt_mass , evtweight, h_1d, nbins, 0, 1600);
        plot1DUnderOverFlow("h_sig_ttRapidity2" + tag_selection + flav_tag, ttRapidity2 , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow("h_sig_tt_pT" + tag_selection + flav_tag, tt_pT , evtweight, h_1d, nbins, 0, 400);
        plot1DUnderOverFlow("h_sig_top1_pt" + tag_selection + flav_tag, top1_pt , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow("h_sig_top2_pt" + tag_selection + flav_tag, top2_pt , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow("h_sig_top1_p_CM" + tag_selection + flav_tag, top1_p_CM , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow("h_sig_top2_p_CM" + tag_selection + flav_tag, top2_p_CM , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow("h_sig_top_rapiditydiffsigned_cms" + tag_selection + flav_tag, top_rapiditydiffsigned_cms , evtweight, h_1d, nbins, -4, 4);
    }


    //default met
    plot1DUnderOverFlow("h_sig_met" + tag_selection + flav_tag, t1metphicorr, evtweight, h_1d, nbins , 0, 500);

    //check HO and TOB/TEC cleanup cut variables
    plot1DUnderOverFlow("h_sig_pfcaloMET" + tag_selection + flav_tag, min(pfcalo_metratio, (float)3.9999) , evtweight, h_1d, 100, 0, 4.);
    plot1DUnderOverFlow("h_sig_pfcalodPhi" + tag_selection + flav_tag, pfcalo_metdphi , evtweight, h_1d, 100, 0, TMath::Pi());

}



void StopTreeLooper::makeAccPlots( float evtweight, std::map<std::string, TH1F *> &h_1d,
                                   string tag_selection, string flav_tag)
{

    int nbins = 80;


}







/*

void StopTreeLooper::makettPlots( float evtweight, std::map<std::string, TH1F *> &h_1d,
                                   string tag_selection, string flav_tag)
{

    double top1dotgen = top1_vecs[imaxweight].Vect().Dot( topplus_genp_p4.Vect() ) / top1_vecs[imaxweight].Vect().Mag() / topplus_genp_p4.Vect().Mag();
    double top1dotgent2 = top1_vecs[imaxweight].Vect().Dot( topminus_genp_p4.Vect() ) / top1_vecs[imaxweight].Vect().Mag() / topminus_genp_p4.Vect().Mag();
    double top2dotgent1 = top2_vecs[imaxweight].Vect().Dot( topplus_genp_p4.Vect() ) / top2_vecs[imaxweight].Vect().Mag() / topplus_genp_p4.Vect().Mag();
    double top2dotgen = top2_vecs[imaxweight].Vect().Dot( topminus_genp_p4.Vect() ) / top2_vecs[imaxweight].Vect().Mag() / topminus_genp_p4.Vect().Mag();
    double top1Pratio = ( top1_vecs[imaxweight].Vect().Mag() - topplus_genp_p4.Vect().Mag() ) / ( top1_vecs[imaxweight].Vect().Mag() + topplus_genp_p4.Vect().Mag() );
    double top2Pratio = ( top2_vecs[imaxweight].Vect().Mag() - topminus_genp_p4.Vect().Mag() ) / ( top2_vecs[imaxweight].Vect().Mag() + topminus_genp_p4.Vect().Mag() );


    double nu1dotgen = nu1_vecs[imaxweight].Vect().Dot( nuPlus_gen.Vect() ) / nu1_vecs[imaxweight].Vect().Mag() / nuPlus_gen.Vect().Mag();
    double nu1dotgennu2 = nu1_vecs[imaxweight].Vect().Dot( nuMinus_gen.Vect() ) / nu1_vecs[imaxweight].Vect().Mag() / nuMinus_gen.Vect().Mag();
    double nu2dotgennu1 = nu2_vecs[imaxweight].Vect().Dot( nuPlus_gen.Vect() ) / nu2_vecs[imaxweight].Vect().Mag() / nuPlus_gen.Vect().Mag();
    double nu2dotgen = nu2_vecs[imaxweight].Vect().Dot( nuMinus_gen.Vect() ) / nu2_vecs[imaxweight].Vect().Mag() / nuMinus_gen.Vect().Mag();
    double nu1Pratio = ( nu1_vecs[imaxweight].Vect().Mag() - nuPlus_gen.Vect().Mag() ) / ( nu1_vecs[imaxweight].Vect().Mag() + nuPlus_gen.Vect().Mag() );
    double nu2Pratio = ( nu2_vecs[imaxweight].Vect().Mag() - nuMinus_gen.Vect().Mag() ) / ( nu2_vecs[imaxweight].Vect().Mag() + nuMinus_gen.Vect().Mag() );


    TLorentzVector nusum_gen = nuPlus_gen + nuMinus_gen;

    double met_x = p_met.first*cos(p_met.second);
    double met_y = p_met.first*sin(p_met.second);

    double DeltaMETsol_gen = sqrt( pow( nusum.Px() - nusum_gen.Px() , 2 ) + pow( nusum.Py() - nusum_gen.Py() , 2 ) );
    double DeltaMETmeas_gen = sqrt( pow( met_x - nusum_gen.Px() , 2 ) + pow( met_y - nusum_gen.Py() , 2 ) );

    //cout<<"about to plot2D"<<endl;

    if(myType == 0 ) {
        plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET_ee", acos(top1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET_ee", acos(top2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);                                
    }
    if(myType == 1 ) {
        plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET_mm", acos(top1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET_mm", acos(top2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);                                
    }
    if(myType == 2 ) {
        plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET_em", acos(top1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET_em", acos(top2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);                                
    }

    plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET", acos(top1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
    plot2DUnderOverFlow(prefix+"_topdotgen_vs_MET", acos(top2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);

    plot2DUnderOverFlow(prefix+"_topdotgen_vs_closestDeltaMET", acos(top1dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);
    plot2DUnderOverFlow(prefix+"_topdotgen_vs_closestDeltaMET", acos(top2dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);

    plot2DUnderOverFlow(prefix+"_topdotgen_vs_njets", acos(top1dotgen), nJets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_topdotgen_vs_njets", acos(top2dotgen), nJets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_topdotgen_vs_nbtag", acos(top1dotgen), nBtagJets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);
    plot2DUnderOverFlow(prefix+"_topdotgen_vs_nbtag", acos(top2dotgen), nBtagJets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);

    plot2DUnderOverFlow(prefix+"_topdotgen_vs_ntaugen", acos(top1dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);
    plot2DUnderOverFlow(prefix+"_topdotgen_vs_ntaugen", acos(top2dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);

    plot2DUnderOverFlow(prefix+"_topdotgen_vs_mttgen", acos(top1dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);
    plot2DUnderOverFlow(prefix+"_topdotgen_vs_mttgen", acos(top2dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);

    if(myType == 0 ) {
        plot2DUnderOverFlow(prefix+"_topPratio_vs_MET_ee", fabs(top1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_topPratio_vs_MET_ee", fabs(top2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);                                
    }
    if(myType == 1 ) {
        plot2DUnderOverFlow(prefix+"_topPratio_vs_MET_mm", fabs(top1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_topPratio_vs_MET_mm", fabs(top2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);                                
    }
    if(myType == 2 ) {
        plot2DUnderOverFlow(prefix+"_topPratio_vs_MET_em", fabs(top1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_topPratio_vs_MET_em", fabs(top2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);                                
    }

    plot2DUnderOverFlow(prefix+"_topPratio_vs_MET", fabs(top1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
    plot2DUnderOverFlow(prefix+"_topPratio_vs_MET", fabs(top2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);

    plot2DUnderOverFlow(prefix+"_topPratio_vs_closestDeltaMET", fabs(top1Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);
    plot2DUnderOverFlow(prefix+"_topPratio_vs_closestDeltaMET", fabs(top2Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);

    plot2DUnderOverFlow(prefix+"_topPratio_vs_njets", fabs(top1Pratio), nJets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_topPratio_vs_njets", fabs(top2Pratio), nJets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_topPratio_vs_nbtag", fabs(top1Pratio), nBtagJets, 1., h_2d, 40, 0., 1., 4, 0., 4.);
    plot2DUnderOverFlow(prefix+"_topPratio_vs_nbtag", fabs(top2Pratio), nBtagJets, 1., h_2d, 40, 0., 1., 4, 0., 4.);

    plot2DUnderOverFlow(prefix+"_topPratio_vs_ntaugen", fabs(top1Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);
    plot2DUnderOverFlow(prefix+"_topPratio_vs_ntaugen", fabs(top2Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);

    plot2DUnderOverFlow(prefix+"_topPratio_vs_mttgen", fabs(top1Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);
    plot2DUnderOverFlow(prefix+"_topPratio_vs_mttgen", fabs(top2Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);



    if(myType == 0 ) {
        plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET_ee", acos(nu1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET_ee", acos(nu2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);                                
    }
    if(myType == 1 ) {
        plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET_mm", acos(nu1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET_mm", acos(nu2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);                                
    }
    if(myType == 2 ) {
        plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET_em", acos(nu1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET_em", acos(nu2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);                                
    }

    plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET", acos(nu1dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
    plot2DUnderOverFlow(prefix+"_nudotgen_vs_MET", acos(nu2dotgen), p_met.first, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);

    plot2DUnderOverFlow(prefix+"_nudotgen_vs_closestDeltaMET", acos(nu1dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);
    plot2DUnderOverFlow(prefix+"_nudotgen_vs_closestDeltaMET", acos(nu2dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);

    plot2DUnderOverFlow(prefix+"_nudotgen_vs_njets", acos(nu1dotgen), nJets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_nudotgen_vs_njets", acos(nu2dotgen), nJets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_nudotgen_vs_nbtag", acos(nu1dotgen), nBtagJets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);
    plot2DUnderOverFlow(prefix+"_nudotgen_vs_nbtag", acos(nu2dotgen), nBtagJets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);

    plot2DUnderOverFlow(prefix+"_nudotgen_vs_ntaugen", acos(nu1dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);
    plot2DUnderOverFlow(prefix+"_nudotgen_vs_ntaugen", acos(nu2dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);

    plot2DUnderOverFlow(prefix+"_nudotgen_vs_mttgen", acos(nu1dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);
    plot2DUnderOverFlow(prefix+"_nudotgen_vs_mttgen", acos(nu2dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);

    if(myType == 0 ) {
        plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET_ee", fabs(nu1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET_ee", fabs(nu2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);                                
    }
    if(myType == 1 ) {
        plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET_mm", fabs(nu1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET_mm", fabs(nu2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);                                
    }
    if(myType == 2 ) {
        plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET_em", fabs(nu1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
        plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET_em", fabs(nu2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);                                
    }

    plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET", fabs(nu1Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);
    plot2DUnderOverFlow(prefix+"_nuPratio_vs_MET", fabs(nu2Pratio), p_met.first, 1., h_2d, 40, 0., 1., 40, 0., 200);

    plot2DUnderOverFlow(prefix+"_nuPratio_vs_closestDeltaMET", fabs(nu1Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);
    plot2DUnderOverFlow(prefix+"_nuPratio_vs_closestDeltaMET", fabs(nu2Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);

    plot2DUnderOverFlow(prefix+"_nuPratio_vs_njets", fabs(nu1Pratio), nJets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_nuPratio_vs_njets", fabs(nu2Pratio), nJets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
    plot2DUnderOverFlow(prefix+"_nuPratio_vs_nbtag", fabs(nu1Pratio), nBtagJets, 1., h_2d, 40, 0., 1., 4, 0., 4.);
    plot2DUnderOverFlow(prefix+"_nuPratio_vs_nbtag", fabs(nu2Pratio), nBtagJets, 1., h_2d, 40, 0., 1., 4, 0., 4.);

    plot2DUnderOverFlow(prefix+"_nuPratio_vs_ntaugen", fabs(nu1Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);
    plot2DUnderOverFlow(prefix+"_nuPratio_vs_ntaugen", fabs(nu2Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);

    plot2DUnderOverFlow(prefix+"_nuPratio_vs_mttgen", fabs(nu1Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);
    plot2DUnderOverFlow(prefix+"_nuPratio_vs_mttgen", fabs(nu2Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);





    if (top1sdp > -998) plot1DUnderOverFlow(prefix+"_top1sdp", acos(top1sdp), 1., h_1d, 40, 0., 3.1415926536);
    if (top2sdp > -998) plot1DUnderOverFlow(prefix+"_top2sdp", acos(top2sdp), 1., h_1d, 40, 0., 3.1415926536);

    plot1DUnderOverFlow(prefix+"_top1dotgen", acos(top1dotgen), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_top2dotgen", acos(top2dotgen), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_top2dotgent1", acos(top2dotgent1), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_top1dotgent2", acos(top1dotgent2), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_top1Pratio", top1Pratio, 1., h_1d, 40, -1., 1.);
    plot1DUnderOverFlow(prefix+"_top2Pratio", top2Pratio, 1., h_1d, 40, -1., 1.);

    plot1DUnderOverFlow(prefix+"_nu1dotgen", acos(nu1dotgen), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_nu2dotgen", acos(nu2dotgen), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_nu2dotgennu1", acos(nu2dotgennu1), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_nu1dotgennu2", acos(nu1dotgennu2), 1., h_1d, 40, 0., 3.1415926536);
    plot1DUnderOverFlow(prefix+"_nu1Pratio", nu1Pratio, 1., h_1d, 40, -1., 1.);
    plot1DUnderOverFlow(prefix+"_nu2Pratio", nu2Pratio, 1., h_1d, 40, -1., 1.);

    plot1DUnderOverFlow(prefix+"_DeltaMETsol_gen",  DeltaMETsol_gen , 1., h_1d, 40, 0., 200.);
    plot1DUnderOverFlow(prefix+"_DeltaMETmeas_gen",  DeltaMETmeas_gen , 1., h_1d, 40, 0., 200.);

    if(closestApproach) {
        plot1DUnderOverFlow(prefix+"_top1dotgen_closest_bestsol", acos(top1dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_top2dotgen_closest_bestsol", acos(top2dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_top1Pratio_closest_bestsol", top1Pratio, 1., h_1d, 40, -1., 1.);
        plot1DUnderOverFlow(prefix+"_top2Pratio_closest_bestsol", top2Pratio, 1., h_1d, 40, -1., 1.);

        plot1DUnderOverFlow(prefix+"_nu1dotgen_closest_bestsol", acos(nu1dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_nu2dotgen_closest_bestsol", acos(nu2dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_nu1Pratio_closest_bestsol", nu1Pratio, 1., h_1d, 40, -1., 1.);
        plot1DUnderOverFlow(prefix+"_nu2Pratio_closest_bestsol", nu2Pratio, 1., h_1d, 40, -1., 1.);

        plot1DUnderOverFlow(prefix+"_closestDeltaMET_closest_maxwsol",      closestDeltaMET_maxwcombo , 1., h_1d, 40, 0., 200.);
        plot1DUnderOverFlow(prefix+"_closestDeltaMET_closest_bestsol",      closestDeltaMET_bestcombo , 1., h_1d, 40, 0., 200.);
        if(closestDeltaMET_othercombo>0) plot1DUnderOverFlow(prefix+"_closestDeltaMET_closest_othersol",      closestDeltaMET_othercombo , 1., h_1d, 40, 0., 200.);

        plot1DUnderOverFlow(prefix+"_DeltaMETsol_gen_closest",  DeltaMETsol_gen , 1., h_1d, 40, 0., 200.);
        plot1DUnderOverFlow(prefix+"_DeltaMETmeas_gen_closest",  DeltaMETmeas_gen , 1., h_1d, 40, 0., 200.);

        for (int i = 0; i < top1_vecs.size(); ++i)
        {
            if (i!=imaxweight) {
                double top1dotgen_i = top1_vecs[i].Vect().Dot( topplus_genp_p4.Vect() ) / top1_vecs[i].Vect().Mag() / topplus_genp_p4.Vect().Mag();
                double top2dotgen_i = top2_vecs[i].Vect().Dot( topminus_genp_p4.Vect() ) / top2_vecs[i].Vect().Mag() / topminus_genp_p4.Vect().Mag();
                double top1Pratio_i = ( top1_vecs[i].Vect().Mag() - topplus_genp_p4.Vect().Mag() ) / ( top1_vecs[i].Vect().Mag() + topplus_genp_p4.Vect().Mag() );
                double top2Pratio_i = ( top2_vecs[i].Vect().Mag() - topminus_genp_p4.Vect().Mag() ) / ( top2_vecs[i].Vect().Mag() + topminus_genp_p4.Vect().Mag() );
                plot1DUnderOverFlow(prefix+"_top1dotgen_closest_othersol", acos(top1dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                plot1DUnderOverFlow(prefix+"_top2dotgen_closest_othersol", acos(top2dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                plot1DUnderOverFlow(prefix+"_top1Pratio_closest_othersol", top1Pratio_i, 1., h_1d, 40, -1., 1.);
                plot1DUnderOverFlow(prefix+"_top2Pratio_closest_othersol", top2Pratio_i, 1., h_1d, 40, -1., 1.);
            }
        }


    }
    else {
        plot1DUnderOverFlow(prefix+"_top1dotgen_max", acos(top1dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_top2dotgen_max", acos(top2dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_top1Pratio_max", top1Pratio, 1., h_1d, 40, -1., 1.);
        plot1DUnderOverFlow(prefix+"_top2Pratio_max", top2Pratio, 1., h_1d, 40, -1., 1.);

        plot1DUnderOverFlow(prefix+"_nu1dotgen_max", acos(nu1dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_nu2dotgen_max", acos(nu2dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow(prefix+"_nu1Pratio_max", nu1Pratio, 1., h_1d, 40, -1., 1.);
        plot1DUnderOverFlow(prefix+"_nu2Pratio_max", nu2Pratio, 1., h_1d, 40, -1., 1.);

        plot1DUnderOverFlow(prefix+"_closestDeltaMET_max",      closestDeltaMET_maxwcombo , 1., h_1d, 40, 0., 200.);

        plot1DUnderOverFlow(prefix+"_DeltaMETsol_gen_max",  DeltaMETsol_gen , 1., h_1d, 40, 0., 200.);
        plot1DUnderOverFlow(prefix+"_DeltaMETmeas_gen_max",  DeltaMETmeas_gen , 1., h_1d, 40, 0., 200.);

        for (int i = 0; i < top1_vecs.size(); ++i)
        {
            if (i!=imaxweight) {
                double top1dotgen_i = top1_vecs[i].Vect().Dot( topplus_genp_p4.Vect() ) / top1_vecs[i].Vect().Mag() / topplus_genp_p4.Vect().Mag();
                double top2dotgen_i = top2_vecs[i].Vect().Dot( topminus_genp_p4.Vect() ) / top2_vecs[i].Vect().Mag() / topminus_genp_p4.Vect().Mag();
                double top1Pratio_i = ( top1_vecs[i].Vect().Mag() - topplus_genp_p4.Vect().Mag() ) / ( top1_vecs[i].Vect().Mag() + topplus_genp_p4.Vect().Mag() );
                double top2Pratio_i = ( top2_vecs[i].Vect().Mag() - topminus_genp_p4.Vect().Mag() ) / ( top2_vecs[i].Vect().Mag() + topminus_genp_p4.Vect().Mag() );
                plot1DUnderOverFlow(prefix+"_top1dotgen_othersols", acos(top1dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                plot1DUnderOverFlow(prefix+"_top2dotgen_othersols", acos(top2dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                plot1DUnderOverFlow(prefix+"_top1Pratio_othersols", top1Pratio_i, 1., h_1d, 40, -1., 1.);
                plot1DUnderOverFlow(prefix+"_top2Pratio_othersols", top2Pratio_i, 1., h_1d, 40, -1., 1.);
            }
        }



    }


}

*/
















void StopTreeLooper::makeNJPlots( float evtweight, std::map<std::string, TH1F *> &h_1d,
                                  string tag_selection, string flav_tag )
{

    plot1D("h_njets"    + tag_selection,          min(n_jets, 4), evtweight, h_1d, 4, 1, 5);
    plot1D("h_njets"    + tag_selection + flav_tag, min(n_jets, 4), evtweight, h_1d, 4, 1, 5);
    plot1D("h_njets_all" + tag_selection,          min(n_jets, 8), evtweight, h_1d, 7, 1, 8);
    plot1D("h_njets_all" + tag_selection + flav_tag, min(n_jets, 8), evtweight, h_1d, 7, 1, 8);

}

void StopTreeLooper::makeZPlots( float evtweight, std::map<std::string, TH1F *> &h_1d,
                                 string tag_selection, string flav_tag )
{

    int nbins = 30;
    float h_xmin = 0.;
    float h_xmax = 300.;

    plot1DUnderOverFlow("h_z_met" + tag_selection + flav_tag, t1metphicorr, evtweight, h_1d, nbins, h_xmin, h_xmax);

    string lep1type =  abs(stopt.id1()) == 13 ? "h_muo" : "h_ele";
    string lep2type =  abs(stopt.id2()) == 13 ? "h_muo" : "h_ele";
    plot1D(lep1type + "pt" + tag_selection + flav_tag, min(stopt.lep1().Pt(), (float)199.99), evtweight, h_1d, 40, 20, 200);
    plot1D(lep2type + "pt" + tag_selection + flav_tag, min(stopt.lep2().Pt(), (float)199.99), evtweight, h_1d, 40, 20, 200);

    plot1D("h_z_leppt"  + tag_selection + flav_tag, min(stopt.lep1().Pt(), (float)299.99), evtweight, h_1d, 50, 20., 300.);
    plot1D("h_z_lepeta" + tag_selection + flav_tag, stopt.lep1().Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    plot1D("h_z_lep2pt" + tag_selection + flav_tag, min(stopt.lep2().Pt(), (float)199.99), evtweight, h_1d, 50, 20., 200.);
    plot1D("h_z_lep2eta" + tag_selection + flav_tag, stopt.lep2().Eta(), evtweight, h_1d, 24, -2.4, 2.4);

    float dphi_metlep = getdphi(stopt.lep1().Phi(), t1metphicorrphi);
    plot1D("h_z_dphi_metl" + tag_selection + flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);

    if ( n_jets < 1 ) return;
    plot1D("h_z_j1pt" + tag_selection + flav_tag, min(jets.at(0).Pt(), (float)399.99), evtweight, h_1d, 20, 30., 400.);
    plot1D("h_z_j1eta" + tag_selection + flav_tag, jets.at(0).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    if ( n_jets < 2 ) return;
    plot1D("h_z_j2pt" + tag_selection + flav_tag, min(jets.at(1).Pt(), (float)299.99), evtweight, h_1d, 20, 30., 300.);
    plot1D("h_z_j2eta" + tag_selection + flav_tag, jets.at(1).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    if ( n_jets < 3 ) return;
    plot1D("h_z_j3pt" + tag_selection + flav_tag, min(jets.at(2).Pt(), (float)199.99), evtweight, h_1d, 20, 30., 200.);
    plot1D("h_z_j3eta" + tag_selection + flav_tag, jets.at(2).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    if ( n_jets < 4 ) return;
    plot1D("h_z_j4pt" + tag_selection + flav_tag, min(jets.at(3).Pt(), (float)119.99), evtweight, h_1d, 20, 30., 120.);
    plot1D("h_z_j4eta" + tag_selection + flav_tag, jets.at(3).Eta(), evtweight, h_1d, 24, -2.4, 2.4);

}





void StopTreeLooper::solvettbar()
{
    //configure solver
    bool useMaxCombo = false; //use lepton-jet combo with maximum sum of weights. false seems to give slightly better resolution.
    bool useClosestDeltaMET = true; //when both combos have only closest-approach solutions, take the one closest to the measured MET instead of the one with the highest weight. true seems to give slightly better resolution.
    bool doDeltaMETcut = true; //reject events where the minimum difference between the solved MET and measured MET exceeds the cut below
    double deltaMETcut = 80.;

    //internal variables
    double maxweight = -1;
    int imaxweightcombos[2] = {-1,-1};
    double maxweightcombos[2] = {-1,-1};
    double avgweightcombos[2] = {0,0};

    //now repeat using the Betchart solver
    double nusols[4][2][3];
    double met_x = t1metphicorr*cos(t1metphicorrphi);
    double met_y = t1metphicorr*sin(t1metphicorrphi);

    //create lines of python code to transmit the input values
    TString l0 = Form("l0 = lv(%0.8f,%0.8f,%0.8f,%0.8f)",lepPlus.Pt(),lepPlus.Eta(),lepPlus.Phi(),lepPlus.E());
    TString l1 = Form("l1 = lv(%0.8f,%0.8f,%0.8f,%0.8f)",lepMinus.Pt(),lepMinus.Eta(),lepMinus.Phi(),lepMinus.E());
    TString j0 = Form("j0 = lv(%0.8f,%0.8f,%0.8f,%0.8f)",bcandidates.at(0).Pt(),bcandidates.at(0).Eta(),bcandidates.at(0).Phi(),bcandidates.at(0).E());
    TString j1 = Form("j1 = lv(%0.8f,%0.8f,%0.8f,%0.8f)",bcandidates.at(1).Pt(),bcandidates.at(1).Eta(),bcandidates.at(1).Phi(),bcandidates.at(1).E());
    TString metxy = Form("metx, mety = %0.8f, %0.8f",met_x,met_y);

    //cout<<l0<<endl;
    //cout<<l1<<endl;
    //cout<<j0<<endl;
    //cout<<j1<<endl;
    //cout<<metxy<<endl;

    TPython::Exec(l0);
    TPython::Exec(l1);
    TPython::Exec(j0);
    TPython::Exec(j1);
    TPython::Exec(metxy);

    int ncombo0 = 0;
    int ncombo1 = 0;
    int num_err_sols[2] = {0,0};

    for (int icombo = 0; icombo < 2; ++icombo)
    {
        if(icombo == 0) TPython::Exec("dnsC = doubleNeutrinoSolutionsCheckLinAlg((j0, j1), (l0, l1), (metx, mety))");
        if(icombo == 1) TPython::Exec("dnsC = doubleNeutrinoSolutionsCheckLinAlg((j1, j0), (l0, l1), (metx, mety))");
        TPython::Exec("dns = dnsC.dns");
        //TPython::Exec("dns = doubleNeutrinoSolutions((j0, j1), (l0, l1), (metx, mety))");
        TPython::Exec("soltest = 1");
        TPython::Exec("if dns==0: soltest = 0");
        //TPython::Exec("print soltest");
        int soltest  = TPython::Eval("soltest");


        if(soltest) {


            TPython::Exec("solutions = dns.nunu_s");
            //TPython::Exec("print solutions");
            TPython::Exec("nSolB = len(solutions)");

            const int nSolB  = TPython::Eval("nSolB");
            //cout<<"nSolB: "<<nSolB<<endl;

            for (int is = 0; is < nSolB; ++is)
            {
                for (int inu = 0; inu < 2; ++inu)
                {
                    for (int ix = 0; ix < 3; ++ix)
                    {
                        TString sols = Form("solutions[%0d][%0d][%0d]",is,inu,ix);
                        //cout<<sols<<endl;
                        nusols[is][inu][ix]  = TPython::Eval(sols);
                        //if(nusols[is][inu][ix] == -1) cout<<nusols[is][inu][ix]<<endl;
                    }
                }
                TLorentzVector nu1_vec , nu2_vec, lvTop1, lvTop2;
                nu1_vec.SetXYZM( nusols[is][0][0] , nusols[is][0][1] , nusols[is][0][2] , 0 );
                nu2_vec.SetXYZM( nusols[is][1][0] , nusols[is][1][1] , nusols[is][1][2] , 0 );

                map<double, double >  mapJetPhi2Discr;
                double sol_weight = -1;
                //double sol_weight_check = -1;

                if(icombo==0) { 
                    lvTop1 = lepPlus + nu1_vec + jet1;
                    lvTop2 = lepMinus + nu2_vec + jet2;
                    //sol_weight = d_llsol->get_weight(jet1 , jet2, lepPlus, lepMinus, nu1_vec, nu2_vec, AMWTmass, mapJetPhi2Discr);
                    sol_weight = 1.;
                    ncombo0++;
                }
                if(icombo==1) {
                    lvTop1 = lepPlus + nu1_vec + jet2;
                    lvTop2 = lepMinus + nu2_vec + jet1;
                    //sol_weight = d_llsol->get_weight(jet2 , jet1, lepPlus, lepMinus, nu1_vec, nu2_vec, AMWTmass, mapJetPhi2Discr);
                    sol_weight = 1.;
                    ncombo1++;
                }


                TLorentzVector lvW1 = lepPlus + nu1_vec;
                TLorentzVector lvW2 = lepMinus + nu2_vec;
                //cout<<"combo "<<icombo<<" solution "<<is<<" weight "<<sol_weight<<" masses: "<<lvTop1.M()<<" "<<lvTop2.M()<<" "<<lvW1.M()<<" "<<lvW2.M()<<endl;

                //don't use solutions with numerical error in solution (output masses don't match input). The input masses used are hard-coded in nuSolutions.py.
                if (  (fabs(172.5 - lvTop1.M()) > 1.0 || fabs(172.5 - lvTop2.M()) > 1.0) ) {
                    num_err_sols[icombo]++;
                    continue;
                }

                if (  (fabs(80.385 - lvW1.M()) > 1.0 || fabs(80.385 - lvW2.M()) > 1.0) ) {
                    //cout<<"mW error: "<<lvW1.M()<<" "<<lvW2.M()<<" mt: "<<lvTop1.M()<<" "<<lvTop2.M()<<endl;
                    num_err_sols[icombo]++;
                    continue;
                }

                nu1_vecs.push_back(nu1_vec);
                nu2_vecs.push_back(nu2_vec);

                top1_vecs.push_back(lvTop1);
                top2_vecs.push_back(lvTop2);
                AMWT_weights.push_back(sol_weight);
                //cout<<"w: "<<sol_weight<<endl;

            }

        }
    }//icombo



    int ncombo0_filled = ncombo0 - num_err_sols[0];
    int ncombo1_filled = ncombo1 - num_err_sols[1];

    //cout<<AMWT_weights.size()<<endl;
    //if(AMWT_weights.size() < 1) cout<<AMWT_weights.size()<<endl;

    for (int is = 0; is < int(AMWT_weights.size()); ++is)
    {
        //cout<<"w: "<<AMWT_weights[is]<<endl;
        if (AMWT_weights[is]>maxweight) {
            imaxweight = is;
            maxweight = AMWT_weights[is];
        }
        if( is < ncombo0_filled ) {
            avgweightcombos[0] += AMWT_weights[is];
            if (AMWT_weights[is]>maxweightcombos[0]) {
                imaxweightcombos[0] = is;
                maxweightcombos[0] = AMWT_weights[is];
            }
        }
        else {
            avgweightcombos[1] += AMWT_weights[is];
            if (AMWT_weights[is]>maxweightcombos[1]) {
                imaxweightcombos[1] = is;
                maxweightcombos[1] = AMWT_weights[is];
            }
        }
    }

    //using average instead of sum with useMaxCombo gives slightly worse resolution, so commented out
    //if(ncombo0_filled>0) avgweightcombos[0] /= ncombo0_filled;
    //if(ncombo1_filled>0) avgweightcombos[1] /= ncombo1_filled;

    //cout<<AMWT_weights.size() <<" "<<ncombo0_filled <<" "<<ncombo1_filled <<" "<<imaxweight<<" "<<maxweight<<" "<<imaxweightcombos[0]<<" "<<maxweightcombos[0]<<" "<<avgweightcombos[0]<<" "<<imaxweightcombos[1]<<" "<<maxweightcombos[1]<<" "<<avgweightcombos[1]<<endl;

    if(useMaxCombo && ncombo0_filled > 0 && ncombo1_filled > 0 ) imaxweight = (avgweightcombos[0] > avgweightcombos[1]) ? imaxweightcombos[0] : imaxweightcombos[1];

    //if( (ncombo0 == 1 || ncombo1 == 1) && (ncombo0 > 1 || ncombo1 > 1) && (ncombo0 - num_err_sols[0] == 0 || ncombo1 - num_err_sols[1] == 0)  )  cout<<"all exact solutions have numerr: "<<ncombo0<<" "<<num_err_sols[0]<<" "<<ncombo1<<" "<<num_err_sols[1]<<" imax: "<<imaxweight<<" "<<imaxweightcombos[0]<<" "<<imaxweightcombos[1]<<endl;

    //don't take "closest approach" solution if exact solutions are available
    if(ncombo0 == 1) {
        if(ncombo1 > 1 && ncombo1_filled > 0) imaxweight = imaxweightcombos[1];
        else closestApproach = true;
    }
    if(ncombo1 == 1) {
        if(ncombo0 > 1 && ncombo0_filled > 0) imaxweight = imaxweightcombos[0];
        else closestApproach = true;
    }

    //cout<<imaxweight<<endl;

    if( AMWT_weights.size() > 0) {
        nusum = nu1_vecs[imaxweight]+nu2_vecs[imaxweight];
        closestDeltaMET_maxwcombo = sqrt( pow( nusum.Px() - met_x , 2 ) + pow( nusum.Py() - met_y , 2 ) );

        if(closestApproach && ncombo0 == 1 && ncombo1 == 1) {
            TLorentzVector nusum_othercombo = nu1_vecs[1-imaxweight]+nu2_vecs[1-imaxweight];
            closestDeltaMET_othercombo = sqrt( pow( nusum_othercombo.Px() - met_x , 2 ) + pow( nusum_othercombo.Py() - met_y , 2 ) );
            if(useClosestDeltaMET && closestDeltaMET_othercombo<closestDeltaMET_maxwcombo ) imaxweight = 1-imaxweight;
            if(doDeltaMETcut && !useClosestDeltaMET && closestDeltaMET_maxwcombo > deltaMETcut && closestDeltaMET_othercombo<closestDeltaMET_maxwcombo ) imaxweight = 1-imaxweight;
        }

        m_top_B = (top1_vecs[imaxweight].M() + top2_vecs[imaxweight].M())/2.;
        nusum = nu1_vecs[imaxweight]+nu2_vecs[imaxweight];
        closestDeltaMET_bestcombo = sqrt( pow( nusum.Px() - met_x , 2 ) + pow( nusum.Py() - met_y , 2 ) );

        //cut on deltaMET (measured vs closest solution)
        if(doDeltaMETcut && closestDeltaMET_bestcombo > deltaMETcut) m_top_B = -999.;

        if(imaxweight < 0) cout<<"something went wrong choosing the best solution"<<endl;

    }

    if(m_top_B > 0) {
        top1_p4 = top1_vecs.at(imaxweight);
        top2_p4 = top2_vecs.at(imaxweight);
    }

}
