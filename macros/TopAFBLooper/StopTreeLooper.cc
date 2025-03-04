#include "StopTreeLooper.h"

//#include "../../CORE/jetSmearingTools.h"
//#include "../../CORE/Thrust.h"
//#include "../../CORE/EventShape.h"
#include "TStopwatch.h"

#include "Math/VectorUtil.h"
#include "STOPT.h"
#include "stopUtils.h"
#include "../Plotting/PlotUtilities.h"
#include "../../Tools/BTagReshaping/BTagReshaping.h"
#include "LHAPDF/LHAPDF.h"

#include "TROOT.h"
#include "TH1D.h"
#include "TH2D.h"
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

#include "TRandom3.h"

#include "denominator/acceptanceplots.h"

using namespace Stop;

std::set<DorkyEventIdentifier> already_seen;
std::set<DorkyEventIdentifier> events_lasercalib;
std::set<DorkyEventIdentifier> events_hcallasercalib;

struct indCSV{
    float CSVdisc;
    int ind;
};

typedef vector< indCSV > VofiCSV;

inline bool sortICSVByCSV(indCSV iCSV1, indCSV iCSV2) {
    return iCSV1.CSVdisc > iCSV2.CSVdisc;
}

    
bool makeCRplots = true; //For CR studies. If true don't apply the selection until making plots.
bool doTobTecVeto = true; //Veto events in TOB/TEC transition region with (charged multiplicity) - (neutral multiplity) > 40 (due to large number of spurious fake tracks due to algoritm problem in 5_X)

//control systematic variations
bool applyTopPtWeighting = false;
bool scaleJESMETUp              = false;
bool scaleJESMETDown            = false;
bool scaleJERUp                 = false;
bool scaleJERDown               = false;
bool scaleLeptonEnergyUp = false;
bool scaleLeptonEnergyDown = false;
//bool scaleUnclusteredEnergyUp = true;
//bool scaleUnclusteredEnergyDown = false;

bool systupBCShape = false;
bool systdownBCShape = false;
bool systupLShape = false;
bool systdownLShape = false;
bool scaleTrigSFup = false;
bool scaleTrigSFdown = false;

bool systupPUShape = false;
bool systdownPUShape = false;
bool systnoPUReweighting = false;
bool systidisoeffweighting = false;

bool weighttaudecay = false;
bool calculatePDFsystweights = false;

//these values are also hard-coded in nuSolutions.py
const double mb_solver = 4.8;
const double mW_solver = 80.385;
const double mt_solver = 172.5;
const double mlb_max = sqrt(mt_solver * mt_solver - mW_solver * mW_solver);

//temporary variables for plotting
double lep1b_mindR = 9999.;
double lep1b_mindPhi = 9999.;
double lep2b_mindR = 9999.;
double lep2b_mindPhi = 9999.;
double lep_dR = -9999;
double lep_dEta = -9999;
double Mll = -9999;
double Etall = -9999;
double Phill = -9999;
double Ptll = -9999;
double t1metphicorrsmeared = -9999;


StopTreeLooper::StopTreeLooper()
{
    m_outfilename_ = "histos.root";
    min_njets = -9999;
    t1metphicorr = -9999.;
    t1metphicorrphi = -9999.;
    met_x = -999999.;
    met_y = -999999.;
    n_jets  = -9999;
    n_bjets = -9999;
    n_ljets = -9999;
    pfcalo_metratio = -9999.;
    pfcalo_metdphi  = -9999.;
    pfcalo_deltamet  = -9999.;
}

StopTreeLooper::~StopTreeLooper()
{
    delete babyFile_;
    //delete babyTree_; //this causes a crash
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

    //split DY sample into ee,mm and tautau so we can use different SFs later
    TString origname = name;
    bool doDYtautau = false;
    bool doDYeemm = false;
    if(name.Contains("DY1to4J") && name.Contains("tautau")) doDYtautau = true; 
    if(name.Contains("DY1to4J") && name.Contains("eemm")) doDYeemm = true;
    if(name.Contains("DY1to4J")) name = "DY1to4Jtot";

    TString babyFilename = origname + "_baby.root";
    MakeBabyNtuple(babyFilename.Data());

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
    // load Betchart solver, and PDFs for AMWT weight
    //------------------------------------------------------------------------------------------------------

    TPython::LoadMacro("loadBetchart.py");
    //LHAPDF::initPDFSet("pdfs/cteq6mE.LHgrid");


    LHAPDF::initPDFSetM(1,"pdfs/cteq6mE.LHgrid");
    LHAPDF::initPDFM(1, 0);
    LHAPDF::initPDFSetM(2,"pdfs/cteq6mE.LHgrid");


    //------------------------------------------------------------------------------------------------------
    // set csv discriminator reshaping
    //------------------------------------------------------------------------------------------------------

    BTagShapeInterface *nominalShape = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", 0.0, 0.0);
    //systematic variations for payloads
    BTagShapeInterface *upBCShape    = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  1.0 ,  0.0);
    BTagShapeInterface *downBCShape  = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root", -1.0 ,  0.0);
    BTagShapeInterface *upLShape     = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  0.0 ,  1.0);
    BTagShapeInterface *downLShape   = new BTagShapeInterface("../../Tools/BTagReshaping/csvdiscr.root",  0.0 , -1.0);

    //------------------------------
    // set up histograms
    //------------------------------

    gROOT->cd();

    cout << "[StopTreeLooper::loop] setting up histos" << endl;

    std::map<std::string, TH1D *> h_1d;
    //also for control regions
    std::map<std::string, TH1D *> h_1d_cr1, h_1d_cr1v, h_1d_cr2, h_1d_cr2v, h_1d_cr0, h_1d_cr3, h_1d_cr3v, h_1d_cr4, h_1d_cr4v, h_1d_cr40, h_1d_cr5, h_1d_cr5v, h_1d_cr6, h_1d_cr6v;
    //for signal region
    std::map<std::string, TH1D *> h_1d_sig;
    std::map<std::string, TH2D *> h_2d_sig;
    //for ttbar dilepton njets distribution
    std::map<std::string, TH1D *> h_1d_nj;
    //z sample for yields etc
    std::map<std::string, TH1D *> h_1d_z;

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

    TH1F *h_pu_wgt_systup = (TH1F *)pu_file->Get("puWeights");
    h_pu_wgt_systup->SetName("h_pu_wgt_systup");

    TH1F *h_pu_wgt_systdown = (TH1F *)pu_file->Get("puWeights");
    h_pu_wgt_systdown->SetName("h_pu_wgt_systdown");


    for (Int_t i = 1; i <= h_pu_wgt->GetNbinsX(); i++)
    {

        double max_nvtx = 14.;
        double widthparameter = 5.;
        double errorweightparameter = (h_pu_wgt->GetBinCenter(i) - max_nvtx)/widthparameter;
        if(errorweightparameter > 1.) errorweightparameter=1.;
        if(errorweightparameter <-1.) errorweightparameter=-1.;
        errorweightparameter *= 2.; //multiply uncertainty by two to make it visible because it is so small

        if (h_pu_wgt_systup->GetBinContent(i) != 0)
        {
            h_pu_wgt_systup->SetBinContent(i, h_pu_wgt_systup->GetBinContent(i) + h_pu_wgt_systup->GetBinError(i)*errorweightparameter );
            h_pu_wgt_systup->SetBinError  (i, h_pu_wgt_systup->GetBinError(i));
        }

        if (h_pu_wgt_systdown->GetBinContent(i) != 0)
        {
            h_pu_wgt_systdown->SetBinContent(i, h_pu_wgt_systdown->GetBinContent(i) - h_pu_wgt_systdown->GetBinError(i)*errorweightparameter );
            h_pu_wgt_systdown->SetBinError  (i, h_pu_wgt_systdown->GetBinError(i));
        }

    }

    TCanvas *canPUweights = new TCanvas("canPUweights","canPUweights",800, 600);
    h_pu_wgt->Draw("hist");
    h_pu_wgt_systdown->SetLineColor(kRed);
    h_pu_wgt_systup->SetLineColor(kGreen-2);
    h_pu_wgt_systdown->Draw("hist sames"); 
    h_pu_wgt_systup->Draw("hist sames"); 
    canPUweights->Print("PUweights.pdf");

    //------------------------------
    // file loop
    //------------------------------

    unsigned int nEventsChain = 0;
    unsigned int nEvents_channel_migrated = 0;
    unsigned int nEvents = chain->GetEntries();
    unsigned int ntotal[10] = {0};
    unsigned int ncorrect[10] = {0};
    unsigned int nwithsol[10] = {0};
    unsigned int ncorrectwithsol[10] = {0};
    unsigned int ninaccept[10] = {0};
    unsigned int ninacceptwithsol[10] = {0};
    nEventsChain = nEvents;
    ULong64_t nEventsTotal = 0;
    int nevt_check = 0;
    int nevt_nlep[10] = {0};

    bool isData = name.Contains("data") ? true : false;

    //Define jet multiplicity requirement
    min_njets = 0;
    printf("[StopTreeLooper::loop] N JET min. requirement %i \n", min_njets);

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
        //nEvents = 100;
        for (ULong64_t event = 0; event < nEvents; ++event)
        {

            stopt.GetEntry(event);

            if(doDYtautau && abs(stopt.mcid1()) != 15) continue;
            if(doDYeemm && abs(stopt.mcid1()) == 15) continue;

            //if (event % 100 != 0) continue; //to skip 99 of every 100 events

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
                        cout << "At " << i_permille / 10. << "% of " << origname.Data() << " time is " << stwatch.RealTime() << " cpu: " << stwatch.CpuTime() << endl;
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
            nevt_nlep[stopt.ngoodlep()]++;

            run = stopt.run();
            ls = stopt.lumi();
            evt = stopt.event();
            weight = -999;

            //----------------------------------------------------------------------------
            // determine event weight
            // make 2 example histograms of nvtx and corresponding weight
            //----------------------------------------------------------------------------

            // to reweight from the nvtx distribution
            // float evtweight = isData ? 1. :
            //    ( stopt.weight() * 19.5 * stopt.nvtxweight() );
            float puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt, isData );
            if (systupPUShape) puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt_systup, isData );
            if (systdownPUShape) puweight = vtxweight_n( stopt.ntruepu(), h_pu_wgt_systdown, isData );
            if (systnoPUReweighting) puweight = 1.;


            float evtweight = isData ? 1. :
                              ( stopt.weight() * 19.5 * puweight );

            //make corrections to the mcatnlo baby weights (same as mgcor in singleLeptonLooper.cc). stopt.nleps() does not include hadronic taus, so count leptons based on their gen-level id
            if (name.Contains("mcatnlo")) {
                  if( stopt.lepPlus_status3_id() == -9999 && stopt.lepMinus_status3_id() == -9999 ) evtweight *= 1.028;
                  else if( stopt.lepPlus_status3_id() != -9999 && stopt.lepMinus_status3_id() != -9999 ) evtweight *= 0.945;
                  else evtweight *= 0.986;
            }


            plot1DUnderOverFlow("h_vtx",       stopt.nvtx(), evtweight, h_1d, 40, 0, 40);
            plot1DUnderOverFlow("h_vtxweight",     puweight, evtweight, h_1d, 41, -4., 4.);

            //----------------------------------------------------------------------------
            // apply preselection:
            // rho 0-40 GeV, MET filters, >=1 good lepton, veto 2 leptons dR < 0.1
            //----------------------------------------------------------------------------

            //if ( !passEvtSelection(name) ) continue; //cut on pfcalo_metdphi is enabled by default (i.e. second argument = true), which does not make sense for events with MET~0 like we can have in the emu channel.
            if ( !passEvtSelection(name, false) ) continue;

            //----------------------------------------------------------------------------
            // Function to perform MET phi corrections on-the-fly
            // Default branches are: tree->t1metphicorr_ and tree->t1metphicorrmt_
            //----------------------------------------------------------------------------

            pair<float, float> p_t1metphicorr =
                getPhiCorrMET( stopt.t1met10(), stopt.t1met10phi(), stopt.nvtx(), !isData);
            t1metphicorr    = p_t1metphicorr.first;
            t1metphicorrphi = p_t1metphicorr.second;
            met_x = t1metphicorr * cos(t1metphicorrphi);
            met_y = t1metphicorr * sin(t1metphicorrphi);

            t1metphicorrsmeared = stopt.t1met10s();

            //cout<<t1metphicorrphi - stopt.t1metphicorrphi()<<endl;

            //pfcalo_metratio = t1metphicorr / stopt.calomet();
            pfcalo_metdphi  = getdphi(t1metphicorrphi, stopt.calometphi());
            pfcalo_deltamet = sqrt( pow( t1metphicorr * sin(t1metphicorrphi) - stopt.calomet() * sin(stopt.calometphi()) , 2 ) + pow( t1metphicorr * cos(t1metphicorrphi) - stopt.calomet() * cos(stopt.calometphi()) , 2 ) );
            pfcalo_metratio = pfcalo_deltamet / t1metphicorr;


            //----------------------------------------------------------------------------
            //evaluate weights for PDF systematic variations
            //----------------------------------------------------------------------------
            PDFsystweights.clear();

            if ( calculatePDFsystweights && name.Contains("tt") ){

                  float   x1          = 0.0;      // momentum fraction for parton1
                  float   x2          = 0.0;
                  int     id1         = 0;        // pdgid of parton1
                  int     id2         = 0;
                  float   Q           = 0.0;      // event momentum scale
                  x1          = stopt.pdfx1();
                  x2          = stopt.pdfx2();
                  id1         = stopt.pdfid1();
                  id2         = stopt.pdfid2();
                  Q           = stopt.pdfQ();

                for (unsigned int subset = 0; subset < 41; subset++)
                {

                  // std::cout << "doing set, subset: " << set_ << ", " << subset << std::endl;
                  //LHAPDF::initPDFM(2,0);

                  LHAPDF::initPDFM(2, subset);

                  // generated pdf values
                  double fx1Q0gen = LHAPDF::xfxM(1, x1, Q, id1) / x1;
                  double fx2Q0gen = LHAPDF::xfxM(1, x2, Q, id2) / x2;
                  // subset pdf values
                  double fx1Qi = LHAPDF::xfxM(2, x1, Q, id1) / x1;
                  double fx2Qi = LHAPDF::xfxM(2, x2, Q, id2) / x2;
                  // calculate weight and fill histogram
                  PDFsystweights.push_back( ((fx1Qi*fx2Qi)/(fx1Q0gen*fx2Q0gen)) );

                }// end of loop over subset of PDFs
            }



            //----------------------------------------------------------------------------------------------------------------------------------------------
            //evaluate weights to correct unpolarised tau decay in mcatnlo (only really correct for 3-body tau decay, so use for systematic estimate only)
            //----------------------------------------------------------------------------------------------------------------------------------------------

            if ( weighttaudecay && name.Contains("mcatnlo") && name.Contains("ttdl") ) 
            {
                double weight_taudecay = 1.;

                TLorentzVector lepPlus_gen_status3(0, 0, 0, 0), lepMinus_gen_status3(0, 0, 0, 0), lepPlus_gen_status1(0, 0, 0, 0), lepMinus_gen_status1(0, 0, 0, 0);
                bool lepPlusIsTau = false;
                bool lepMinusIsTau = false;

                if( abs(stopt.lepPlus_status3_id()) == 15 ) {
                    lepPlusIsTau = true;
                    lepPlus_gen_status1.SetXYZT( stopt.lepPlus_status1().x(), stopt.lepPlus_status1().y(), stopt.lepPlus_status1().z(), stopt.lepPlus_status1().t() );
                    lepPlus_gen_status3.SetXYZT( stopt.lepPlus_status3_orig().x(), stopt.lepPlus_status3_orig().y(), stopt.lepPlus_status3_orig().z(), stopt.lepPlus_status3_orig().t() );
                }

                if( abs(stopt.lepMinus_status3_id()) == 15 ) {
                    lepMinusIsTau = true;
                    lepMinus_gen_status1.SetXYZT( stopt.lepMinus_status1().x(), stopt.lepMinus_status1().y(), stopt.lepMinus_status1().z(), stopt.lepMinus_status1().t() );
                    lepMinus_gen_status3.SetXYZT( stopt.lepMinus_status3_orig().x(), stopt.lepMinus_status3_orig().y(), stopt.lepMinus_status3_orig().z(), stopt.lepMinus_status3_orig().t() );
                }

                if (lepPlusIsTau)
                {
                    lepPlus_gen_status1.Boost(-lepPlus_gen_status3.BoostVector());
                    double cosTheta_lepPlus_status1 = lepPlus_gen_status1.Vect().Dot(lepPlus_gen_status3.Vect()) / (lepPlus_gen_status1.Vect().Mag() * lepPlus_gen_status3.Vect().Mag());
                    double EoverEmax_lepPlus = 2.*lepPlus_gen_status1.E() / lepPlus_gen_status3.M();
                    if (EoverEmax_lepPlus > 1.) EoverEmax_lepPlus = 1.;
                    //if(EoverEmax_lepPlus>1.) {cout<<"**********x>1********"<<endl; cout<<EoverEmax_lepPlus<<" "<<lepPlus_gen_status1.E()<<" "<<lepPlus_gen_status3.M()<<endl;}
                    //double weight_lepPlus = 1.+cosTheta_lepPlus_status1/3.;  //approximation: cosTheta dependence varies with x
                    double weight_lepPlus = 1. + cosTheta_lepPlus_status1 * ( 1. - 2.*EoverEmax_lepPlus ) / ( 2.*EoverEmax_lepPlus - 3. );
                    //cout<<"P "<<(weighttaudecay?weight_lepPlus:1.)<<endl;
                    weight_taudecay *= weight_lepPlus;
                    //cout<<weight_lepPlus<<" ";
                }
                if (lepMinusIsTau)
                {
                    lepMinus_gen_status1.Boost(-lepMinus_gen_status3.BoostVector());
                    double cosTheta_lepMinus_status1 = lepMinus_gen_status1.Vect().Dot(lepMinus_gen_status3.Vect()) / (lepMinus_gen_status1.Vect().Mag() * lepMinus_gen_status3.Vect().Mag());
                    double EoverEmax_lepMinus = 2.*lepMinus_gen_status1.E() / lepMinus_gen_status3.M();
                    if (EoverEmax_lepMinus > 1.) EoverEmax_lepMinus = 1.;
                    //if(EoverEmax_lepMinus>1.) {cout<<"**********x>1********"<<endl; cout<<EoverEmax_lepMinus<<" "<<lepMinus_gen_status1.E()<<" "<<lepMinus_gen_status3.M()<<endl;}
                    //double weight_lepMinus = 1.+cosTheta_lepMinus_status1/3.;  //approximation: cosTheta dependence varies with x
                    double weight_lepMinus = 1. + cosTheta_lepMinus_status1 * ( 1. - 2.*EoverEmax_lepMinus ) / ( 2.*EoverEmax_lepMinus - 3. );
                    //cout<<"M "<<(weighttaudecay?weight_lepMinus:1.)<<endl;
                    weight_taudecay *= weight_lepMinus;
                    //cout<<weight_lepMinus<<" ";
                }

                //cout<<weight_taudecay<<endl;
                evtweight *= weight_taudecay;

            } //weighttaudecay



            //----------------------------------------------------------------------------
            // get jet information
            //----------------------------------------------------------------------------

            jets.clear();
            bjets.clear();
            nonbjets.clear();
            bcandidates.clear();
            CSVdisc.clear();
            btag.clear();
            //sigma_jets.clear();
            mc.clear();

            n_jets  = 0;
            n_bjets = 0;
            n_ljets = 0;
            tobtecveto_ = false;

            //evaluate reshaped CSV values (for MC) and fill vector for sorting
            VofiCSV viCSV;
            viCSV.clear();

            for ( unsigned int ijet = 0 ; ijet < stopt.pfjets().size() ; ++ijet ) {
                float CSVdisc_temp = isData ? stopt.pfjets_csv().at(ijet)
                                    : nominalShape->reshape( stopt.pfjets().at(ijet).eta(),
                                            stopt.pfjets().at(ijet).pt(),
                                            stopt.pfjets_csv().at(ijet),
                                            stopt.pfjets_mcflavorAlgo().at(ijet) );

                if( systupBCShape && !isData ) CSVdisc_temp = upBCShape->reshape( stopt.pfjets().at(ijet).eta(), stopt.pfjets().at(ijet).pt(), stopt.pfjets_csv().at(ijet), stopt.pfjets_mcflavorAlgo().at(ijet) );
                if( systdownBCShape && !isData ) CSVdisc_temp = downBCShape->reshape( stopt.pfjets().at(ijet).eta(), stopt.pfjets().at(ijet).pt(), stopt.pfjets_csv().at(ijet), stopt.pfjets_mcflavorAlgo().at(ijet) );
                if( systupLShape && !isData ) CSVdisc_temp = upLShape->reshape( stopt.pfjets().at(ijet).eta(), stopt.pfjets().at(ijet).pt(), stopt.pfjets_csv().at(ijet), stopt.pfjets_mcflavorAlgo().at(ijet) );
                if( systdownLShape && !isData ) CSVdisc_temp = downLShape->reshape( stopt.pfjets().at(ijet).eta(), stopt.pfjets().at(ijet).pt(), stopt.pfjets_csv().at(ijet), stopt.pfjets_mcflavorAlgo().at(ijet) );

                indCSV iCSV = { CSVdisc_temp, ijet };
                CSVdisc.push_back( CSVdisc_temp );
                viCSV.push_back(iCSV);
            }

            sort(viCSV.begin(), viCSV.end(), sortICSVByCSV);

            //stopt.pfjets() are sorted by pT. Sort by CSV instead.
            //for ( unsigned int i = 0 ; i < stopt.pfjets().size() ; ++i )
            for ( unsigned int iii = 0 ; iii < stopt.pfjets().size() ; ++iii )
            {
                int i = viCSV.at(iii).ind;

                if ( stopt.pfjets().at(i).pt() < 30 )  continue;
                if ( fabs(stopt.pfjets().at(i).eta()) > 2.4 )  continue; //note, this was 2.5 for 7 TeV AFB
                //  if( stopt.pfjets_beta2_0p5().at(i)<0.2 )  continue;

                float sigma = stopt.pfjets_sigma().at(i);
                float unc   = stopt.pfjets_uncertainty().at(i);

                //smear the JER
                if (  !isData ) {

                    TRandom3 r3(  stopt.event()*1000 + i  );
                    double JER_smear_width = 0.;
                    if(scaleJERUp)            JER_smear_width = sqrt( pow( getDataMCRatioSystUp(stopt.pfjets().at(i).eta()) , 2 ) - 1. ) * sigma  / getDataMCRatioOld(stopt.pfjets().at(i).eta());  //pfjets_sigma() is already scaled by getDataMCRatio(Old) in the babies for MC, so must divide
                    else if(scaleJERDown)     JER_smear_width = sqrt( pow( getDataMCRatioSystDown(stopt.pfjets().at(i).eta()) , 2 ) - 1. ) * sigma  / getDataMCRatioOld(stopt.pfjets().at(i).eta());  //pfjets_sigma() is already scaled by getDataMCRatio(Old) in the babies for MC, so must divide
                    else                      JER_smear_width = sqrt( pow( getDataMCRatio(stopt.pfjets().at(i).eta()) , 2 ) - 1. ) * sigma  / getDataMCRatioOld(stopt.pfjets().at(i).eta());  //pfjets_sigma() is already scaled by getDataMCRatio(Old) in the babies for MC, so must divide

                    double JER_smear = r3.Gaus(0,JER_smear_width);

                    //cout<<stopt.pfjets().at(i)<<endl;

                    met_x -= JER_smear*stopt.pfjets().at(i).px();
                    met_y -= JER_smear*stopt.pfjets().at(i).py();
                    stopt.pfjets().at(i) *= (1. + JER_smear);

                    //cout<<stopt.pfjets().at(i)<<" "<<JER_smear_width<<" "<<JER_smear<<" "<<endl;
                }

                if ( (scaleJESMETDown || scaleJESMETUp) && isData ) {

                    if(scaleJESMETDown==scaleJESMETUp) {cout<<"error: can't vary JES up and down at the same time"<<endl; exit(0);}

                    float vardirection = scaleJESMETUp?1:-1;

                    //cout<<stopt.pfjets().at(i)<<endl;

                    //met_x -= vardirection*unc*stopt.pfjets().at(i).px(); //only jets with pt>30 GeV are stored, so can't recompute the deltaMET. Using the value from the stop babies instead.
                    //met_y -= vardirection*unc*stopt.pfjets().at(i).py(); //only jets with pt>30 GeV are stored, so can't recompute the deltaMET. Using the value from the stop babies instead.
                    stopt.pfjets().at(i) *= (1. + vardirection*unc);

                    //cout<<stopt.pfjets().at(i)<<" "<<vardirection<<" "<<unc<<" "<<endl;
                }

                TString currentfilename = currentFile->GetTitle();

                //  bool passMediumPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(),stopt.pfjets_mvaPUid().at(i),1);
                bool passTightPUid = true;
                //hack due to missing pfjets_mva5xPUid in the low mass DY baby
                if(! currentfilename.Contains("Mll10to50")  ) passTightPUid = passMVAJetId(stopt.pfjets().at(i).pt(), stopt.pfjets().at(i).eta(), stopt.pfjets_mva5xPUid().at(i), 0);

                if (!passTightPUid) continue;

                n_jets++;
                n_ljets++;

                jets.push_back( stopt.pfjets().at(i) );

                // check for tob/tec tracking problem
                if ((fabs(stopt.pfjets().at(i).eta()) > 0.9) && (fabs(stopt.pfjets().at(i).eta()) < 1.9))
                {
                    if ((stopt.pfjets_chm().at(i) - stopt.pfjets_neu().at(i)) > 40) tobtecveto_ = true;
                }

                float csv_nominal = CSVdisc.at(i);

                if (csv_nominal > 0.679)
                {
                    bjets.push_back( stopt.pfjets().at(i) );
                    //if more than 2 b-jets, use the ones with highest pT as b candidates
                    if (bcandidates.size() < 2) bcandidates.push_back( stopt.pfjets().at(i) );
                    n_bjets++;
                }
                else
                {
                    nonbjets.push_back( stopt.pfjets().at(i) );
                }
                btag.push_back( csv_nominal );

                if ( !isData ) mc.push_back  ( stopt.pfjets_mc3().at(i) );
                else mc.push_back  ( 0 );

                //sigma_jets.push_back(sigma);

                //count jets that are not overlapping with second lepton
                if (isData) continue;
                if (stopt.nleps() < 2) continue;
                if (stopt.mclep2().pt() < 30.) continue;
                if (ROOT::Math::VectorUtil::DeltaR(stopt.mclep2(), stopt.pfjets().at(i)) > 0.4 ) continue;
                n_ljets--;

            }

            if (n_jets < min_njets) continue;


            if ( (scaleLeptonEnergyDown || scaleLeptonEnergyUp) && isData ) {

                if(scaleLeptonEnergyDown==scaleLeptonEnergyUp) {cout<<"error: can't vary LES up and down at the same time"<<endl; exit(0);}

                float vardirection = scaleLeptonEnergyUp?1:-1;

                if( abs(stopt.id1())==11 ) {
                    float Lunc = abs(stopt.lep1().eta()) > 1.4442 ? 0.015 : 0.006;
                    met_x -= vardirection*Lunc*stopt.lep1().px();
                    met_y -= vardirection*Lunc*stopt.lep1().py();
                    stopt.lep1() *= (1. + vardirection*Lunc);
                    if ( stopt.id1() > 0 ) stopt.lepp() *= (1. + vardirection*Lunc);
                    else stopt.lepm() *= (1. + vardirection*Lunc);
                    //cout<<(1. + vardirection*Lunc)<<endl;
                }

                if( abs(stopt.id2())==11 ) {
                    float Lunc = abs(stopt.lep2().eta()) > 1.4442 ? 0.015 : 0.006;
                    met_x -= vardirection*Lunc*stopt.lep2().px();
                    met_y -= vardirection*Lunc*stopt.lep2().py();
                    stopt.lep2() *= (1. + vardirection*Lunc);
                    if ( stopt.id2() > 0 ) stopt.lepp() *= (1. + vardirection*Lunc);
                    else stopt.lepm() *= (1. + vardirection*Lunc);
                    //cout<<(1. + vardirection*Lunc)<<endl;
                }


            }

/*
//insufficient information in babies to do unclustered energy variation separately. It's included in t1metphicorrup and t1metphicorrdn.
            if ( scaleUnclusteredEnergyUp || scaleUnclusteredEnergyDown ) {
                float pfmetx = stopt.t1metphicorr() * cos( stopt.t1metphicorrphi() );
                float pfmety = stopt.t1metphicorr() * sin( stopt.t1metphicorrphi() );
                float pfmetxup = stopt.t1metphicorrup() * cos( stopt.t1metphicorrphiup() );
                float pfmetyup = stopt.t1metphicorrup() * sin( stopt.t1metphicorrphiup() );
                float pfmetxdn = stopt.t1metphicorrdn() * cos( stopt.t1metphicorrphidn() );
                float pfmetydn = stopt.t1metphicorrdn() * sin( stopt.t1metphicorrphidn() );

                cout<<pfmetxup-pfmetx<<" "<<pfmetxdn-pfmetx<<" "<<met_x<<" "<<pfmetxup<<" "<<pfmetx<<" "<<pfmetxdn<<endl;
                //met_x - pfmetxup
            }
*/

            //if making systematic variations, recalculate t1metphicorr using the updated met_x, met_y  (for MC we are always varying the JER to the smeared central value)
            if( ( (scaleLeptonEnergyDown || scaleLeptonEnergyUp) && isData ) || (!isData) ) {
                t1metphicorr = sqrt( met_x*met_x + met_y*met_y );
                t1metphicorrphi = atan2( met_y , met_x );
            }

            if( scaleJESMETDown ) {
                t1metphicorr = stopt.t1metphicorrdn();
                t1metphicorrphi = stopt.t1metphicorrphidn();
                met_x = t1metphicorr * cos(t1metphicorrphi);
                met_y = t1metphicorr * sin(t1metphicorrphi);
            }

            if( scaleJESMETUp ) {
                t1metphicorr = stopt.t1metphicorrup();
                t1metphicorrphi = stopt.t1metphicorrphiup();
                met_x = t1metphicorr * cos(t1metphicorrphi);
                met_y = t1metphicorr * sin(t1metphicorrphi);
            }

            //if <2 btagged jets, take the highest pT light jets as b candidates
            if (bcandidates.size() < 2 && nonbjets.size() > 0)
            {
                bcandidates.push_back( nonbjets.at(0) );
            }
            if (bcandidates.size() < 2 && nonbjets.size() > 1)
            {
                bcandidates.push_back( nonbjets.at(1) );
            }

             lep1b_mindR = 9999.;
             lep1b_mindPhi = 9999.;
             lep2b_mindR = 9999.;
             lep2b_mindPhi = 9999.;

            for (int i = 0; i < bcandidates.size(); ++i)
            {
                if (ROOT::Math::VectorUtil::DeltaR(stopt.lep1(), bcandidates.at(i)) < lep1b_mindR ) lep1b_mindR = ROOT::Math::VectorUtil::DeltaR(stopt.lep1(), bcandidates.at(i));
                if (ROOT::Math::VectorUtil::DeltaR(stopt.lep2(), bcandidates.at(i)) < lep2b_mindR ) lep2b_mindR = ROOT::Math::VectorUtil::DeltaR(stopt.lep2(), bcandidates.at(i));
                if (fabs(ROOT::Math::VectorUtil::DeltaPhi(stopt.lep1(), bcandidates.at(i))) < lep1b_mindPhi ) lep1b_mindPhi = fabs(ROOT::Math::VectorUtil::DeltaPhi(stopt.lep1(), bcandidates.at(i)));
                if (fabs(ROOT::Math::VectorUtil::DeltaPhi(stopt.lep2(), bcandidates.at(i))) < lep2b_mindPhi ) lep2b_mindPhi = fabs(ROOT::Math::VectorUtil::DeltaPhi(stopt.lep2(), bcandidates.at(i)));
            }
            //if(lep1b_mindPhi<0.6 || lep2b_mindPhi<0.6 ) continue;

            //----------------------------------------------------------------------------
            // Veto events with jet(s) with suspected TOB/TEC seeded tracking problem. http://www.t2.ucsd.edu/tastwiki/pub/CMS/20130426OsChats/wh_tobtecprob_olivito_260413.pdf
            //----------------------------------------------------------------------------

            if (doTobTecVeto && tobtecveto_) continue;

            //store dilepton mass to fill in baby tree (region 20<Mll<30 is worst affected by deltaPhi bump)
            dilmass = stopt.dilmass();

            //mlb variables sensitive to top mass without needing solver
            mlb_1 = bcandidates.size()<1 ? 0 : (stopt.lep1() + bcandidates.at(0)).M();
            mlb_2 = bcandidates.size()<2 ? 0 : (stopt.lep2() + bcandidates.at(1)).M();
            mlb_3 = bcandidates.size()<2 ? 0 : (stopt.lep1() + bcandidates.at(1)).M();
            mlb_4 = bcandidates.size()<1 ? 0 : (stopt.lep2() + bcandidates.at(0)).M();
            mlb_min = bcandidates.size()<1 ? 0 : ( bcandidates.size()<2 ? std::min(mlb_1, mlb_4) : std::min( std::min(mlb_1, mlb_2), std::min(mlb_3, mlb_4) ) );

            //----------------------------------------------------------------------------
            // Require event to pass full selection prior to ttbar solver unless we want to make control region plots.
            //----------------------------------------------------------------------------

            if ( !makeCRplots && !passFullSelection(isData) ) continue;

            //----------------------------------------------------------------------------
            // set up variables used by ttbar solver
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

            if (stopt.id1() > 0 && stopt.ngoodlep() > 1)
            {
                lepPlus.SetPtEtaPhiE(stopt.lep1().Pt(), stopt.lep1().Eta(), stopt.lep1().Phi(), stopt.lep1().E());
                lepMinus.SetPtEtaPhiE(stopt.lep2().Pt(), stopt.lep2().Eta(), stopt.lep2().Phi(), stopt.lep2().E());
            }
            else if (stopt.ngoodlep() > 1)
            {
                lepMinus.SetPtEtaPhiE(stopt.lep1().Pt(), stopt.lep1().Eta(), stopt.lep1().Phi(), stopt.lep1().E());
                lepPlus.SetPtEtaPhiE(stopt.lep2().Pt(), stopt.lep2().Eta(), stopt.lep2().Phi(), stopt.lep2().E());
            }

            //if ( stopt.id1() * stopt.id2() > 0 ) printf("[StopTreeLooper::loop] Same-sign leptons \n");

            if (stopt.ngoodlep()>1) if ( ( fabs(lepPlus.E() - stopt.lepp().E()) > 0.01 || fabs(lepMinus.E() - stopt.lepm().E()) > 0.01 ) && stopt.id1() * stopt.id2() < 0 ) printf("[StopTreeLooper::loop] Something went wrong with lepton assignments. %d %f %f %d %f %f \n", stopt.id1(), lepPlus.E(), stopt.lepp().E(), stopt.id2(), lepMinus.E(), stopt.lepm().E() );

            if(bcandidates.size()>0) jet1.SetPtEtaPhiE(bcandidates.at(0).Pt(), bcandidates.at(0).Eta(), bcandidates.at(0).Phi(), bcandidates.at(0).E());
            else jet1.SetPtEtaPhiM(4, 4., 0., 4.8);  //dummy jet for 0-jet CRs to avoid crash
            if(bcandidates.size()>1) jet2.SetPtEtaPhiE(bcandidates.at(1).Pt(), bcandidates.at(1).Eta(), bcandidates.at(1).Phi(), bcandidates.at(1).E());
            else jet2.SetPtEtaPhiM(4, 4.*jet1.Eta()/fabs(jet1.Eta()), 0., 4.8);  //for 1-jet events assume 2nd jet was soft to allow ttbar solution in 1j CRs

            top1_p4.SetPtEtaPhiE(0, 0, 0, 0);
            top2_p4.SetPtEtaPhiE(0, 0, 0, 0);
            nusum.SetPtEtaPhiE(0, 0, 0, 0);
            cms.SetPtEtaPhiE(0, 0, 0, 0);

            m_top = -999;
            m_top_B = -999;
            closestDeltaMET_maxwcombo = -999;
            closestDeltaMET_othercombo = -999;
            closestDeltaMET_bestcombo = -999;
            imaxweight = -1;
            closestApproach = false;

            //----------------------------------------------------------------------------
            // ttbar solver
            //----------------------------------------------------------------------------

            //if( stopt.ngoodlep() > 1 ) solvettbar(); //this takes a long time to run on DY MC, so apply the Z veto
            if( stopt.ngoodlep() > 1 && n_jets > 0 && (abs(stopt.id1()) != abs(stopt.id2()) || (n_jets > 1 && n_bjets > 0) || fabs( stopt.dilmass() - 91.) > 15. ) ) solvettbar();

            //cout<<"gtl: "<<__LINE__<<endl;

            m_top = m_top_B;

            //----------------------------------------------------------------------------
            // asymmetry calculations
            //----------------------------------------------------------------------------

            lep_charge_asymmetry = -999.0;
            lep_azimuthal_asymmetry = -999.0;
            lep_azimuthal_asymmetry2 = -999.0;
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

            //more variables for plots
            float lep_pseudorap_diff = -999.0;
            float lep_cosalpha = -999.0;
            float jet_azimuthal_asymmetry = -999.0;
            float jet_pseudorap_diff = -999.0;
            float jet_cosalpha = -999.0;

            if(stopt.ngoodlep()>1){
                
                //fully leptonic asymmetries
                lep_charge_asymmetry = abs(lepPlus.Eta()) - abs(lepMinus.Eta());
                lep_azimuthal_asymmetry = lepPlus.DeltaPhi(lepMinus);  //lep_azimuthal_asymmetry is same as lep_azimuthal_asymmetry2 but from -pi to pi instead of folding it over from 0 to pi
                lep_azimuthal_asymmetry2 = fabs(lepPlus.DeltaPhi(lepMinus));
                lep_dR = lepPlus.DeltaR(lepMinus);
                lep_dEta = lepPlus.Eta() - lepMinus.Eta();

                Mll = (lepPlus+lepMinus).M();
                Ptll = (lepPlus+lepMinus).Pt();
                Etall = (lepPlus+lepMinus).Eta();
                Phill = (lepPlus+lepMinus).Phi();

                //more variables for plots
                lep_pseudorap_diff = (lepPlus.Eta()) - (lepMinus.Eta());
                lep_cosalpha =  lepPlus.Vect().Dot( lepMinus.Vect() ) / (lepPlus.Vect().Mag() * lepMinus.Vect().Mag());
            }
            if(bcandidates.size()>1) {
                jet_azimuthal_asymmetry = jet1.DeltaPhi(jet2);
                jet_pseudorap_diff = jet1.Eta() - jet2.Eta();
                jet_cosalpha =  jet1.Vect().Dot( jet2.Vect() ) / (jet1.Vect().Mag() * jet2.Vect().Mag());
            }

            //variables that require ttbar solution
            if ( m_top > 0 )
            {

                if ( fabs(top1_p4.M() - mt_solver) > 1.  || fabs(top2_p4.M() - mt_solver) > 1. ) printf("[StopTreeLooper::loop] Something went wrong with ttbar solver. \n");

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
            // gen-level asymmetry calculations
            //----------------------------------------------------------------------------

            lepPlus_gen.SetPtEtaPhiE(0, 0, 0, 0);
            lepMinus_gen.SetPtEtaPhiE(0, 0, 0, 0);
            topplus_genp_p4.SetPtEtaPhiE(0, 0, 0, 0);
            topminus_genp_p4.SetPtEtaPhiE(0, 0, 0, 0);
            cms_gen.SetPtEtaPhiE(0, 0, 0, 0);

            lep_charge_asymmetry_gen = -999.0;
            lep_azimuthal_asymmetry_gen = -999.0;
            lep_azimuthal_asymmetry2_gen = -999.0;
            top_rapiditydiff_cms_gen = -999.0;
            top_pseudorapiditydiff_cms_gen = -999.0;
            top_rapiditydiff_Marco_gen = -999.0;
            top_costheta_cms_gen = -999.0;
            lepPlus_costheta_cms_gen = -999.0;
            lepMinus_costheta_cms_gen = -999.0;
            top_spin_correlation_gen = -999.0;
            lep_cos_opening_angle_gen = -999.0;
            tt_mass_gen = -999.0;
            ttRapidity2_gen = -999.0;
            tt_pT_gen = -999.0;
            top1_pt_gen = -999.0;
            top2_pt_gen = -999.0;

            //if ( name.Contains("ttdl") )
            if ( name.Contains("ttdl") &&  tree->GetListOfBranches()->FindObject("topPlus_status3") ) //to make it work with the old ttdl stop babies (reco-level only)
            {

                //use original status=3 tops
                topplus_genp_p4.SetPtEtaPhiE(stopt.topPlus_status3().Pt(), stopt.topPlus_status3().Eta(), stopt.topPlus_status3().Phi(), stopt.topPlus_status3().E());
                topminus_genp_p4.SetPtEtaPhiE(stopt.topMinus_status3().Pt(), stopt.topMinus_status3().Eta(), stopt.topMinus_status3().Phi(), stopt.topMinus_status3().E());

                //use leptons corrected to level before FSR from b
                lepPlus_gen.SetPtEtaPhiE(stopt.lepPlus_status3().Pt(), stopt.lepPlus_status3().Eta(), stopt.lepPlus_status3().Phi(), stopt.lepPlus_status3().E());
                lepMinus_gen.SetPtEtaPhiE(stopt.lepMinus_status3().Pt(), stopt.lepMinus_status3().Eta(), stopt.lepMinus_status3().Phi(), stopt.lepMinus_status3().E());

                if ( (fabs(topplus_genp_p4.M() - mt_solver) > 100.  || fabs(topminus_genp_p4.M() - mt_solver) > 100.) ) printf("[StopTreeLooper::loop] Something went wrong with gen-level tops. %f %f %d \n", topplus_genp_p4.M() , topminus_genp_p4.M() , stopt.nleps() );

                /*
                lep_charge_asymmetry_gen = abs(lepPlus_gen.Eta()) - abs(lepMinus_gen.Eta());
                lep_azimuthal_asymmetry_gen = lepPlus_gen.DeltaPhi(lepMinus_gen);
                lep_azimuthal_asymmetry2_gen = acos(cos(lep_azimuthal_asymmetry_gen));

                tt_mass_gen = (topplus_genp_p4 + topminus_genp_p4).M();
                if ( (fabs(tt_mass_gen - stopt.mttbar()) > 1.0) ) printf("[StopTreeLooper::loop] Something went wrong with gen-level ttbar. %f %f \n", tt_mass_gen, stopt.mttbar() );
                ttRapidity2_gen = (topplus_genp_p4 + topminus_genp_p4).Rapidity();

                top1_pt_gen =  topplus_genp_p4.Pt();
                top2_pt_gen =  topminus_genp_p4.Pt();

                top_rapiditydiff_cms_gen = (topplus_genp_p4.Rapidity() - topminus_genp_p4.Rapidity()) * (topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity());
                top_pseudorapiditydiff_cms_gen = abs(topplus_genp_p4.Eta()) - abs(topminus_genp_p4.Eta());
                top_rapiditydiff_Marco_gen = abs(topplus_genp_p4.Rapidity()) - abs(topminus_genp_p4.Rapidity());


                cms_gen = topplus_genp_p4 + topminus_genp_p4;
                tt_pT_gen = cms_gen.Pt();
                topplus_genp_p4.Boost(-cms_gen.BoostVector());
                topminus_genp_p4.Boost(-cms_gen.BoostVector());
                top_costheta_cms_gen = topplus_genp_p4.Vect().Dot(cms_gen.Vect()) / (topplus_genp_p4.Vect().Mag() * cms_gen.Vect().Mag());

                lepPlus_gen.Boost(-cms_gen.BoostVector());
                lepPlus_gen.Boost(-topplus_genp_p4.BoostVector());
                lepMinus_gen.Boost(-cms_gen.BoostVector());
                lepMinus_gen.Boost(-topminus_genp_p4.BoostVector());

                lepPlus_costheta_cms_gen = lepPlus_gen.Vect().Dot(topplus_genp_p4.Vect()) / (lepPlus_gen.Vect().Mag() * topplus_genp_p4.Vect().Mag());
                lepMinus_costheta_cms_gen = lepMinus_gen.Vect().Dot(topminus_genp_p4.Vect()) / (lepMinus_gen.Vect().Mag() * topminus_genp_p4.Vect().Mag());

                top_spin_correlation_gen = lepPlus_costheta_cms_gen * lepMinus_costheta_cms_gen;
                lep_cos_opening_angle_gen = lepPlus_gen.Vect().Dot(lepMinus_gen.Vect()) / (lepPlus_gen.Vect().Mag() * lepMinus_gen.Vect().Mag());

                //reset 4-vectors to lab frame
                topplus_genp_p4.SetPtEtaPhiE(stopt.t().Pt(), stopt.t().Eta(), stopt.t().Phi(), stopt.t().E());
                topminus_genp_p4.SetPtEtaPhiE(stopt.tbar().Pt(), stopt.tbar().Eta(), stopt.tbar().Phi(), stopt.tbar().E());
                lepPlus_gen.SetPtEtaPhiE(stopt.lep_t().Pt(), stopt.lep_t().Eta(), stopt.lep_t().Phi(), stopt.lep_t().E());
                lepMinus_gen.SetPtEtaPhiE(stopt.lep_tbar().Pt(), stopt.lep_tbar().Eta(), stopt.lep_tbar().Phi(), stopt.lep_tbar().E());



                //test values already stored in new babies
                if(fabs(lep_charge_asymmetry_gen  -  stopt.lep_charge_asymmetry_gen())>1e-2 && fabs(1.-lep_charge_asymmetry_gen/stopt.lep_charge_asymmetry_gen())>1e-2) cout<<"problem with stored lep_charge_asymmetry_gen "<<fabs(1.-lep_charge_asymmetry_gen/stopt.lep_charge_asymmetry_gen())<<" "<<stopt.lep_charge_asymmetry_gen()<<" "<<lep_charge_asymmetry_gen<<endl;
                if(fabs(lep_azimuthal_asymmetry_gen  -  stopt.lep_azimuthal_asymmetry_gen())>1e-2 && fabs(1.-lep_azimuthal_asymmetry_gen/stopt.lep_azimuthal_asymmetry_gen())>1e-2) cout<<"problem with stored lep_azimuthal_asymmetry_gen "<<fabs(1.-lep_azimuthal_asymmetry_gen/stopt.lep_azimuthal_asymmetry_gen())<<" "<<stopt.lep_azimuthal_asymmetry_gen()<<" "<<lep_azimuthal_asymmetry_gen<<endl;
                if(fabs(lep_azimuthal_asymmetry2_gen  -  stopt.lep_azimuthal_asymmetry2_gen())>1e-2 && fabs(1.-lep_azimuthal_asymmetry2_gen/stopt.lep_azimuthal_asymmetry2_gen())>1e-2) cout<<"problem with stored lep_azimuthal_asymmetry2_gen "<<fabs(1.-lep_azimuthal_asymmetry2_gen/stopt.lep_azimuthal_asymmetry2_gen())<<" "<<stopt.lep_azimuthal_asymmetry2_gen()<<" "<<lep_azimuthal_asymmetry2_gen<<endl;
                if(fabs(top_rapiditydiff_cms_gen  -  stopt.top_rapiditydiff_cms_gen())>1e-2 && fabs(1.-top_rapiditydiff_cms_gen/stopt.top_rapiditydiff_cms_gen())>1e-2) cout<<"problem with stored top_rapiditydiff_cms_gen "<<fabs(1.-top_rapiditydiff_cms_gen/stopt.top_rapiditydiff_cms_gen())<<" "<<stopt.top_rapiditydiff_cms_gen()<<" "<<top_rapiditydiff_cms_gen<<endl;
                if(fabs(top_pseudorapiditydiff_cms_gen  -  stopt.top_pseudorapiditydiff_cms_gen())>1e-2 && fabs(1.-top_pseudorapiditydiff_cms_gen/stopt.top_pseudorapiditydiff_cms_gen())>1e-2) cout<<"problem with stored top_pseudorapiditydiff_cms_gen "<<fabs(1.-top_pseudorapiditydiff_cms_gen/stopt.top_pseudorapiditydiff_cms_gen())<<" "<<stopt.top_pseudorapiditydiff_cms_gen()<<" "<<top_pseudorapiditydiff_cms_gen<<endl;
                if(fabs(top_rapiditydiff_Marco_gen  -  stopt.top_rapiditydiff_Marco_gen())>1e-2 && fabs(1.-top_rapiditydiff_Marco_gen/stopt.top_rapiditydiff_Marco_gen())>1e-2) cout<<"problem with stored top_rapiditydiff_Marco_gen "<<fabs(1.-top_rapiditydiff_Marco_gen/stopt.top_rapiditydiff_Marco_gen())<<" "<<stopt.top_rapiditydiff_Marco_gen()<<" "<<top_rapiditydiff_Marco_gen<<endl;
                if(fabs(top_costheta_cms_gen  -  stopt.top_costheta_cms_gen())>1e-2 && fabs(1.-top_costheta_cms_gen/stopt.top_costheta_cms_gen())>1e-2) cout<<"problem with stored top_costheta_cms_gen "<<fabs(1.-top_costheta_cms_gen/stopt.top_costheta_cms_gen())<<" "<<stopt.top_costheta_cms_gen()<<" "<<top_costheta_cms_gen<<endl;
                if(fabs(lepPlus_costheta_cms_gen  -  stopt.lepPlus_costheta_cms_gen())>1e-2 && fabs(1.-lepPlus_costheta_cms_gen/stopt.lepPlus_costheta_cms_gen())>1e-2) cout<<"problem with stored lepPlus_costheta_cms_gen "<<fabs(1.-lepPlus_costheta_cms_gen/stopt.lepPlus_costheta_cms_gen())<<" "<<stopt.lepPlus_costheta_cms_gen()<<" "<<lepPlus_costheta_cms_gen<<endl;
                if(fabs(lepMinus_costheta_cms_gen  -  stopt.lepMinus_costheta_cms_gen())>1e-2 && fabs(1.-lepMinus_costheta_cms_gen/stopt.lepMinus_costheta_cms_gen())>1e-2) cout<<"problem with stored lepMinus_costheta_cms_gen "<<fabs(1.-lepMinus_costheta_cms_gen/stopt.lepMinus_costheta_cms_gen())<<" "<<stopt.lepMinus_costheta_cms_gen()<<" "<<lepMinus_costheta_cms_gen<<endl;
                if(fabs(top_spin_correlation_gen  -  stopt.top_spin_correlation_gen())>1e-2 && fabs(1.-top_spin_correlation_gen/stopt.top_spin_correlation_gen())>1e-2) cout<<"problem with stored top_spin_correlation_gen "<<fabs(1.-top_spin_correlation_gen/stopt.top_spin_correlation_gen())<<" "<<stopt.top_spin_correlation_gen()<<" "<<top_spin_correlation_gen<<endl;
                if(fabs(lep_cos_opening_angle_gen  -  stopt.lep_cos_opening_angle_gen())>1e-2 && fabs(1.-lep_cos_opening_angle_gen/stopt.lep_cos_opening_angle_gen())>1e-2) cout<<"problem with stored lep_cos_opening_angle_gen "<<fabs(1.-lep_cos_opening_angle_gen/stopt.lep_cos_opening_angle_gen())<<" "<<stopt.lep_cos_opening_angle_gen()<<" "<<lep_cos_opening_angle_gen<<endl;
                if(fabs(tt_mass_gen  -  stopt.tt_mass_gen())>1e-2 && fabs(1.-tt_mass_gen/stopt.tt_mass_gen())>1e-2) cout<<"problem with stored tt_mass_gen "<<fabs(1.-tt_mass_gen/stopt.tt_mass_gen())<<" "<<stopt.tt_mass_gen()<<" "<<tt_mass_gen<<endl;
                if(fabs(ttRapidity2_gen  -  stopt.ttRapidity2_gen())>1e-2 && fabs(1.-ttRapidity2_gen/stopt.ttRapidity2_gen())>1e-2) cout<<"problem with stored ttRapidity2_gen "<<fabs(1.-ttRapidity2_gen/stopt.ttRapidity2_gen())<<" "<<stopt.ttRapidity2_gen()<<" "<<ttRapidity2_gen<<endl;
                if(fabs(tt_pT_gen  -  stopt.tt_pT_gen())>1e-2 && fabs(1.-tt_pT_gen/stopt.tt_pT_gen())>1e-2) cout<<"problem with stored tt_pT_gen "<<fabs(1.-tt_pT_gen/stopt.tt_pT_gen())<<" "<<stopt.tt_pT_gen()<<" "<<tt_pT_gen<<endl;
                //if(fabs(top1_pt_gen  -  stopt.topPlus_status3().Pt())>1e-2 && fabs(1.-top1_pt_gen/stopt.topPlus_status3().Pt())>1e-2) cout<<"problem with stored top1_pt_gen "<<fabs(1.-top1_pt_gen/stopt.topPlus_status3().Pt())<<" "<<stopt.topPlus_status3().Pt()<<" "<<top1_pt_gen<<endl;
                //if(fabs(top2_pt_gen  -  stopt.topMinus_status3().Pt())>1e-2 && fabs(1.-top2_pt_gen/stopt.topMinus_status3().Pt())>1e-2) cout<<"problem with stored top2_pt_gen "<<fabs(1.-top2_pt_gen/stopt.topMinus_status3().Pt())<<" "<<stopt.topMinus_status3().Pt()<<" "<<top2_pt_gen<<endl;
                if(fabs(stopt.topPlus_status3().Pt()  -  stopt.ptt())>1e-2 && fabs(1.-stopt.topPlus_status3().Pt()/stopt.ptt())>1e-2) cout<<"problem with stored stopt.topPlus_status3().Pt() "<<fabs(1.-stopt.topPlus_status3().Pt()/stopt.ptt())<<" "<<stopt.ptt()<<" "<<stopt.topPlus_status3().Pt()<<endl;
                if(fabs(stopt.topMinus_status3().Pt()  -  stopt.pttbar())>1e-2 && fabs(1.-stopt.topMinus_status3().Pt()/stopt.pttbar())>1e-2) cout<<"problem with stored stopt.topMinus_status3().Pt() "<<fabs(1.-stopt.topMinus_status3().Pt()/stopt.pttbar())<<" "<<stopt.pttbar()<<" "<<stopt.topMinus_status3().Pt()<<endl;
                */


                lep_charge_asymmetry_gen = stopt.lep_charge_asymmetry_gen();
                lep_azimuthal_asymmetry_gen = stopt.lep_azimuthal_asymmetry_gen();
                lep_azimuthal_asymmetry2_gen = stopt.lep_azimuthal_asymmetry2_gen();
                tt_mass_gen = stopt.tt_mass_gen();
                ttRapidity2_gen = stopt.ttRapidity2_gen();
                top_rapiditydiff_cms_gen = stopt.top_rapiditydiff_cms_gen();
                top_pseudorapiditydiff_cms_gen = stopt.top_pseudorapiditydiff_cms_gen();
                top_rapiditydiff_Marco_gen = stopt.top_rapiditydiff_Marco_gen();
                tt_pT_gen = stopt.tt_pT_gen();
                top_costheta_cms_gen = stopt.top_costheta_cms_gen();
                lepPlus_costheta_cms_gen = stopt.lepPlus_costheta_cms_gen();
                lepMinus_costheta_cms_gen = stopt.lepMinus_costheta_cms_gen();
                top_spin_correlation_gen = stopt.top_spin_correlation_gen();
                lep_cos_opening_angle_gen = stopt.lep_cos_opening_angle_gen();

                top1_pt_gen =  topplus_genp_p4.Pt();
                top2_pt_gen =  topminus_genp_p4.Pt();


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
            channel = -999;
            genchannel = -999;
            string flav_tag_sl;
            if ( abs(stopt.id1()) == 13 ) flav_tag_sl = "_muo";
            else if ( abs(stopt.id1()) == 11 ) flav_tag_sl = "_ele";
            else flav_tag_sl = "_mysterysl";
            string flav_tag_dl;
            if      ( abs(stopt.id1()) == abs(stopt.id2()) && abs(stopt.id1()) == 13 ) {
                flav_tag_dl = "_dimu";
                channel = 1;
            }
            else if ( abs(stopt.id1()) == abs(stopt.id2()) && abs(stopt.id1()) == 11 ) {
                flav_tag_dl = "_diel";
                channel = 0;
            }
            else if ( abs(stopt.id1()) != abs(stopt.id2()) && abs(stopt.id1()) == 13 ) {
                flav_tag_dl = "_muel";
                channel = 2;
            }
            else if ( abs(stopt.id1()) != abs(stopt.id2()) && abs(stopt.id1()) == 11 ) {
                flav_tag_dl = "_elmu";
                channel = 2;
            }
            else flav_tag_dl = "_mysterydl";
            string basic_flav_tag_dl = flav_tag_dl;
            if ( abs(stopt.id1()) != abs(stopt.id2()) && flav_tag_dl != "_mysterydl" ) basic_flav_tag_dl = "_mueg";
            //if (stopt.id1()*stopt.id2() > 0) basic_flav_tag_dl += "_SS";


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
            //float trigweight_dl = isData ? 1. : getdltrigweight(stopt.id1(), stopt.id2());
            float trigweight_dl = isData ? 1. : getdltrigweight_pteta(stopt.id1(), stopt.lep1().Pt(), stopt.lep1().Eta(), stopt.id2(), stopt.lep2().Pt(), stopt.lep2().Eta());


            if ( applyTopPtWeighting && (name.Contains("ttdl") || name.Contains("ttsl") || name.Contains("ttfake")) )
            {

                float pT_topplus_gen = stopt.t().Pt();
                float pT_topminus_gen = stopt.tbar().Pt();
                evtweight *= sqrt( TopPtWeight(pT_topplus_gen) * TopPtWeight(pT_topminus_gen) );
            }

            if ( systidisoeffweighting ) {
                double idisoweight = 1.;
                float lep1eta_temp = fabs(stopt.lep1().Eta());
                if ( abs(stopt.id1())==11 )  lep1eta_temp = lep1eta_temp > 1.44 ? 1.44 : lep1eta_temp;
                if ( abs(stopt.id1())==13 )  lep1eta_temp = lep1eta_temp > 2.09 ? 2.09 : lep1eta_temp;

                float lep2eta_temp = fabs(stopt.lep2().Eta());
                if ( abs(stopt.id2())==11 )  lep2eta_temp = lep2eta_temp > 1.44 ? 1.44 : lep2eta_temp;
                if ( abs(stopt.id2())==13 )  lep2eta_temp = lep2eta_temp > 2.09 ? 2.09 : lep2eta_temp;

                idisoweight *= getisoeffweight( stopt.id1(), stopt.lep1().Pt(), lep1eta_temp );
                idisoweight *= getisoeffweight( stopt.id2(), stopt.lep2().Pt(), lep2eta_temp );
                idisoweight *= getideffweight( stopt.id1(), stopt.lep1().Pt(), lep1eta_temp );
                idisoweight *= getideffweight( stopt.id2(), stopt.lep2().Pt(), lep2eta_temp );

                evtweight*=idisoweight;
                //cout<<idisoweight<<endl;
            }

            //check performance of b candidate selection
            if (  dataset_2l && (passFullSelection(isData) || passFullSelection_bveto(isData)) ) {
                //test if selected b candidates match the true b quarks
                if(name.Contains("ttdl")) {

                    int nmatched = 0;

                    if ( ROOT::Math::VectorUtil::DeltaR(stopt.bPlus_status3(), bcandidates.at(0)) < 0.4 ) nmatched++;
                    else if ( ROOT::Math::VectorUtil::DeltaR(stopt.bMinus_status3(), bcandidates.at(0)) < 0.4 ) nmatched++;
                    if ( ROOT::Math::VectorUtil::DeltaR(stopt.bPlus_status3(), bcandidates.at(0)) >= 0.4 ) if ( ROOT::Math::VectorUtil::DeltaR(stopt.bPlus_status3(), bcandidates.at(1)) < 0.4 ) nmatched++;
                    if ( ROOT::Math::VectorUtil::DeltaR(stopt.bMinus_status3(), bcandidates.at(0)) >= 0.4 ) if ( ROOT::Math::VectorUtil::DeltaR(stopt.bPlus_status3(), bcandidates.at(1)) >= 0.4 ) if ( ROOT::Math::VectorUtil::DeltaR(stopt.bMinus_status3(), bcandidates.at(1)) < 0.4 ) nmatched++;

                    if(nmatched==2) ncorrect[n_bjets]++;
                    if(nmatched>2) cout<<"something went wrong with matching"<<endl;
                    if(m_top>0 && nmatched==2) ncorrectwithsol[n_bjets]++;

                    if( fabs(stopt.bPlus_status3().eta()) < 2.4 && fabs(stopt.bPlus_status3().pt())>30. && fabs(stopt.bMinus_status3().eta()) < 2.4 && fabs(stopt.bMinus_status3().pt())>30. ) {
                        ninaccept[n_bjets]++;
                        if(m_top>0) ninacceptwithsol[n_bjets]++;
                    }
                }

                ntotal[n_bjets]++;
                if(m_top>0) nwithsol[n_bjets]++;

                //if (event % 100 == 0 || n_bjets>2) cout<<n_bjets<<" b-jets, within acceptance: "<<double(ninaccept[n_bjets])/double(ntotal[n_bjets])<<", correct permutation: "<<double(ncorrect[n_bjets])/double(ntotal[n_bjets])<<", with solution: "<<double(nwithsol[n_bjets])/double(ntotal[n_bjets])<<", in acceptance when with solution: "<<double(ninacceptwithsol[n_bjets])/double(nwithsol[n_bjets])<<", correct permutation when with solution: "<<double(ncorrectwithsol[n_bjets])/double(nwithsol[n_bjets])<<", ntotal: "<<ntotal[n_bjets]<<endl;
            }

            //veto on 3 or more b-tags because the chance of choosing the correct permutation is worse than in 1-tag events (66% [3tag] vs 68% [1tag] vs 96%[2tag] in events with a solution) and because the purity is questionable (data/MC agreement is bad in >2 b-tag bins).
            if(n_bjets>2) continue;

            //
            // SIGNAL REGION - dilepton + b-tag
            //

            if ( dataset_2l && passFullSelection(isData) )
            {
                //check cuts
                if ( (stopt.lep1().Pt() < 20.0 || stopt.lep2().Pt() < 20.0 ) ) cout << stopt.lep1().Pt() << " " << stopt.lep2().Pt() << endl;
                if ( (stopt.lep1() + stopt.lep2()).M() < 20.0 )  cout << (stopt.lep1() + stopt.lep2()).M() << endl;
                if ( t1metphicorr < 40. && basic_flav_tag_dl != "_mueg" )  cout << "MET: " << t1metphicorr << " " + basic_flav_tag_dl << endl;

                weight = evtweight * trigweight_dl;

                makeNJPlots( weight, h_1d_nj, "", basic_flav_tag_dl);
                makeSIGPlots( weight, h_1d_sig,  tag_btag  , basic_flav_tag_dl , "h_sig");
                makeSIGPlots( weight, h_1d_sig,  tag_btag  , "_all" , "h_sig");

                if ( name.Contains("ttdl") )
                {
                    //gen-level flavour
                    string basic_flav_tag_dl_gen = "_genunknown";
                    if(tree->GetListOfBranches()->FindObject("topPlus_status3")){ //to make it work with the old ttdl stop babies (reco-level only)
                    if ( stopt.lepPlus_status1_id() == -11 && stopt.lepMinus_status1_id() == 11 ) {genchannel = 0; basic_flav_tag_dl_gen = "_gendiel";}
                    else if ( stopt.lepPlus_status1_id() == -13 && stopt.lepMinus_status1_id() == 11 ) {genchannel = 2; basic_flav_tag_dl_gen = "_genmueg";}
                    else if ( stopt.lepPlus_status1_id() == -11 && stopt.lepMinus_status1_id() == 13 ) {genchannel = 2; basic_flav_tag_dl_gen = "_genmueg";}
                    else if ( stopt.lepPlus_status1_id() == -13 && stopt.lepMinus_status1_id() == 13 ) {genchannel = 1; basic_flav_tag_dl_gen = "_gendimu";}
                    else cout<<"indeterminate dilepton type: "<<stopt.lepPlus_status1_id()<<" "<<stopt.lepMinus_status1_id()<<endl;

                    if(channel!=genchannel) {
                        nEvents_channel_migrated++;
                        //cout<<nEvents_channel_migrated<<" "<<100.*double(nEvents_channel_migrated)/double(nEventsTotal)<<"%"<<endl;
                    }
                    }

                    makettPlots( weight, h_1d_sig, h_2d_sig, tag_btag  , basic_flav_tag_dl );
                    makettPlots( weight, h_1d_sig, h_2d_sig, tag_btag  , "_all" );
                    makeAccPlots( weight, h_1d_sig, h_2d_sig, tag_btag  , basic_flav_tag_dl );
                    if(basic_flav_tag_dl == "_mueg") makeAccPlots( weight, h_1d_sig, h_2d_sig, tag_btag  , flav_tag_dl );
                    makeAccPlots( weight, h_1d_sig, h_2d_sig, tag_btag  , basic_flav_tag_dl_gen );
                    makeAccPlots( weight, h_1d_sig, h_2d_sig, tag_btag  , "_all" );
                }

                FillBabyNtuple();

            }

            // Control regions
/*
                  //
                  // CR1 - single lepton + b-veto
                  //

                  // selection - 1 lepton + iso track veto
                  // Add b-tag veto
                  if ( dataset_1l && passSingleLeptonSelection(isData)
                   //&& passisotrk
                   && n_jets>=2 )
                {
                    //pre b-tag veto
                    makeSIGPlots( evtweight*trigweight, h_1d_cr1, "_prebveto", flav_tag_sl , "h_cr1");
                    makeSIGPlots( evtweight*trigweight, h_1d_cr1, "_prebveto", "_all" , "h_cr1");
                    //b-veto
                    if ( n_bjets==0 ) makeSIGPlots( evtweight*trigweight, h_1d_cr1, "" , flav_tag_sl, "h_cr1" );
                    if ( n_bjets==0 ) makeSIGPlots( evtweight*trigweight, h_1d_cr1, "" , "_all", "h_cr1" );

                }//end CR1 selection
*/

            // CR1 - dilepton control REGION - dilepton + b-tag + Z peak
            if ( dataset_2l && passFullSelection_Zpeak_inclusiveb(isData) && n_bjets > 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr1, "" , basic_flav_tag_dl , "h_cr1");
                makeSIGPlots( weight, h_1d_cr1, "" , "_all" , "h_cr1");

            }

            // CR1v - dilepton control REGION - dilepton + b-veto + Z peak
            if ( dataset_2l && passFullSelection_Zpeak_inclusiveb(isData) && n_bjets == 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr1v, "" , basic_flav_tag_dl , "h_cr1v");
                makeSIGPlots( weight, h_1d_cr1v, "" , "_all" , "h_cr1v");

            }

            //
            // CR2(v) - Z-peak for yields and mT resolution studies
            //

            // selection - SF dilepton, veto on isolated track in addition to 2 leptons, in z-peak
            if ( dataset_2l && passDileptonSelectionWithEndcapEls(isData) )
            {

              //invariant mass - basic check of inclusive distribution
              plot1DUnderOverFlow("h_z_dilmass"+flav_tag_dl, stopt.dilmass(), evtweight*trigweight_dl, h_1d_z,  30 , 76 , 106);

              if ( fabs( stopt.dilmass() - 91.) <= 15. )
                {

                  // if (n_jets>8)
                  //    cout<<"NJETS: "<<n_jets<<" * dataset: "<<stopt.dataset()
                  //        <<" run: "<<stopt.run()<<" lumi: "<<stopt.lumi()<<" event: "<<stopt.event()<<endl;

                  //z peak plots
                  plot1DUnderOverFlow("h_z_njets"    +flav_tag_dl, min(n_jets,4),  evtweight*trigweight_dl, h_1d_z, 5,0,5);
                  plot1DUnderOverFlow("h_z_njets_all"+flav_tag_dl, min(n_jets,9),  evtweight*trigweight_dl, h_1d_z, 10, 0, 10);
                  plot1DUnderOverFlow("h_z_nbjets"   +flav_tag_dl, min(n_bjets,3), evtweight*trigweight_dl, h_1d_z, 4, 0, 4);
                  makeZPlots( evtweight*trigweight_dl, h_1d_z, "", flav_tag_dl );

                  // Add stricter 3rd lepton veto
                  // require at least 2 jets
                  //CR2
                  if ( n_jets>=2 && n_bjets > 0 ) {

                  makeSIGPlots( evtweight*trigweight_dl, h_1d_cr2, "", basic_flav_tag_dl , "h_cr2");
                  makeSIGPlots( evtweight*trigweight_dl, h_1d_cr2, "", "_all" , "h_cr2");

                  }
                  //CR2v
                  if ( n_jets>=2 && n_bjets == 0 ) {

                  makeSIGPlots( evtweight*trigweight_dl, h_1d_cr2v, "", basic_flav_tag_dl , "h_cr2v");
                  makeSIGPlots( evtweight*trigweight_dl, h_1d_cr2v, "", "_all" , "h_cr2v");

                  }
                }
            }//end CR2 selection


  
            // CR0 - dilepton control REGION - dilepton + 0 b-tag
            if ( dataset_2l && passFullSelection_bveto(isData) )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr0, "" , basic_flav_tag_dl , "h_cr0");
                makeSIGPlots( weight, h_1d_cr0, "" , "_all" , "h_cr0");

            }

            // CR3 - dilepton control REGION - dilepton + no MET cut + b-tag
            if ( dataset_2l && passFullSelection_noMETcut_inclusiveb(isData) && n_bjets > 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr3, "" , basic_flav_tag_dl , "h_cr3");
                makeSIGPlots( weight, h_1d_cr3, "" , "_all" , "h_cr3");

            }

            // CR3v - dilepton control REGION - dilepton + no MET cut + b-veto
            if ( dataset_2l && passFullSelection_noMETcut_inclusiveb(isData) && n_bjets == 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr3v, "" , basic_flav_tag_dl , "h_cr3v");
                makeSIGPlots( weight, h_1d_cr3v, "" , "_all" , "h_cr3v");

            }

            // CR4 - dilepton control REGION - dilepton + 1 b-tagged jet
            if ( dataset_2l && passFullSelection_1jet_inclusiveb(isData) && n_bjets > 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr4, "" , basic_flav_tag_dl , "h_cr4");
                makeSIGPlots( weight, h_1d_cr4, "" , "_all" , "h_cr4");

            }

            // CR4v - dilepton control REGION - dilepton + 1 untagged jet
            if ( dataset_2l && passFullSelection_1jet_inclusiveb(isData) && n_bjets == 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr4v, "" , basic_flav_tag_dl , "h_cr4v");
                makeSIGPlots( weight, h_1d_cr4v, "" , "_all" , "h_cr4v");

            }

            // CR40 - dilepton control REGION - dilepton + 0 jets
            if ( dataset_2l && passFullSelection_0jets(isData) )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr40, "" , basic_flav_tag_dl , "h_cr40");
                makeSIGPlots( weight, h_1d_cr40, "" , "_all" , "h_cr40");

            }

            // CR5 - dilepton control REGION - SS dilepton + b-tag
            if ( dataset_2l && passFullSelection_SS_inclusiveb(isData) && n_bjets > 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr5, "" , basic_flav_tag_dl , "h_cr5");
                makeSIGPlots( weight, h_1d_cr5, "" , "_all" , "h_cr5");

            }

            // CR5v - dilepton control REGION - SS dilepton + b-veto
            if ( dataset_2l && passFullSelection_SS_inclusiveb(isData) && n_bjets == 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr5v, "" , basic_flav_tag_dl , "h_cr5v");
                makeSIGPlots( weight, h_1d_cr5v, "" , "_all" , "h_cr5v");

            }

            // CR6 - dilepton control REGION - SS dilepton + no MET cut + b-tag
            if ( dataset_2l && passFullSelection_SS_noMETcut_inclusiveb(isData) && n_bjets > 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr6, "" , basic_flav_tag_dl , "h_cr6");
                makeSIGPlots( weight, h_1d_cr6, "" , "_all" , "h_cr6");

            }

            // CR6v - dilepton control REGION - SS dilepton + no MET cut + b-veto
            if ( dataset_2l && passFullSelection_SS_noMETcut_inclusiveb(isData) && n_bjets == 0 )
            {
                weight = evtweight * trigweight_dl;

                makeSIGPlots( weight, h_1d_cr6v, "" , basic_flav_tag_dl , "h_cr6v");
                makeSIGPlots( weight, h_1d_cr6v, "" , "_all" , "h_cr6v");

            }


        } // end event loop

        // delete tree;

        //stwatch.Stop();
        //cout << "time is " << stwatch.CpuTime() << endl;
        //stwatch.Start();

    } // end file loop

    for (int nb = 0; nb < 6; ++nb)
    {
        cout<<nb<<" b-jets, within acceptance: "<<double(ninaccept[nb])/double(ntotal[nb])<<", correct permutation: "<<double(ncorrect[nb])/double(ntotal[nb])<<", with solution: "<<double(nwithsol[nb])/double(ntotal[nb])<<", in acceptance when with solution: "<<double(ninacceptwithsol[nb])/double(nwithsol[nb])<<", correct permutation when with solution: "<<double(ncorrectwithsol[nb])/double(nwithsol[nb])<<", ntotal: "<<ntotal[nb]<<endl;
    }

    //
    // finish
    //
    printf("[StopTreeLooper::loop] Closing baby ntuple \n");
    CloseBabyNtuple();

    cout << "N EVENT CHECK " << nevt_check << endl;
    TFile outfile(m_outfilename_.c_str(), "RECREATE") ;
    printf("[StopTreeLooper::loop] Saving histograms to %s\n", m_outfilename_.c_str());

    std::map<std::string, TH1D *>::iterator it1d;
    for (it1d = h_1d.begin(); it1d != h_1d.end(); it1d++)
    {
        it1d->second->Write();
        delete it1d->second;
    }

    outfile.Write();
    outfile.Close();

    TFile outfile_sig(Form("SIG%s", m_outfilename_.c_str()), "RECREATE") ;
    printf("[StopTreeLooper::loop] Saving SIG histograms to %s\n", m_outfilename_.c_str());

    std::map<std::string, TH1D *>::iterator it1d_sig;
    for (it1d_sig = h_1d_sig.begin(); it1d_sig != h_1d_sig.end(); it1d_sig++)
    {
        it1d_sig->second->Write();
        delete it1d_sig->second;
    }

    std::map<std::string, TH2D *>::iterator it2d_sig;
    for (it2d_sig = h_2d_sig.begin(); it2d_sig != h_2d_sig.end(); it2d_sig++)
    {
        it2d_sig->second->Write();
        delete it2d_sig->second;
    }

    outfile_sig.Write();
    outfile_sig.Close();

      //control regions
      //h_1d_cr*

      TFile outfile_cr1(Form("CR1%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR1 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr1;
      for(it1d_cr1=h_1d_cr1.begin(); it1d_cr1!=h_1d_cr1.end(); it1d_cr1++) {
        it1d_cr1->second->Write();
        delete it1d_cr1->second;
      }

      outfile_cr1.Write();
      outfile_cr1.Close();

      TFile outfile_cr1v(Form("CR1v%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR1v histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr1v;
      for(it1d_cr1v=h_1d_cr1v.begin(); it1d_cr1v!=h_1d_cr1v.end(); it1d_cr1v++) {
        it1d_cr1v->second->Write();
        delete it1d_cr1v->second;
      }

      outfile_cr1v.Write();
      outfile_cr1v.Close();


      TFile outfile_cr2(Form("CR2%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR2 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr2;
      for(it1d_cr2=h_1d_cr2.begin(); it1d_cr2!=h_1d_cr2.end(); it1d_cr2++) {
        it1d_cr2->second->Write();
        delete it1d_cr2->second;
      }

      outfile_cr2.Write();
      outfile_cr2.Close();

      TFile outfile_cr2v(Form("CR2v%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR2v histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr2v;
      for(it1d_cr2v=h_1d_cr2v.begin(); it1d_cr2v!=h_1d_cr2v.end(); it1d_cr2v++) {
        it1d_cr2v->second->Write();
        delete it1d_cr2v->second;
      }

      outfile_cr2v.Write();
      outfile_cr2v.Close();


      TFile outfile_cr0(Form("CR0%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR0 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr0;
      for(it1d_cr0=h_1d_cr0.begin(); it1d_cr0!=h_1d_cr0.end(); it1d_cr0++) {
        it1d_cr0->second->Write();
        delete it1d_cr0->second;
      }

      outfile_cr0.Write();
      outfile_cr0.Close();


      TFile outfile_cr3(Form("CR3%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR3 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr3;
      for(it1d_cr3=h_1d_cr3.begin(); it1d_cr3!=h_1d_cr3.end(); it1d_cr3++) {
        it1d_cr3->second->Write();
        delete it1d_cr3->second;
      }

      outfile_cr3.Write();
      outfile_cr3.Close();

      TFile outfile_cr3v(Form("CR3v%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR3v histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr3v;
      for(it1d_cr3v=h_1d_cr3v.begin(); it1d_cr3v!=h_1d_cr3v.end(); it1d_cr3v++) {
        it1d_cr3v->second->Write();
        delete it1d_cr3v->second;
      }

      outfile_cr3v.Write();
      outfile_cr3v.Close();


      TFile outfile_cr4(Form("CR4%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR4 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr4;
      for(it1d_cr4=h_1d_cr4.begin(); it1d_cr4!=h_1d_cr4.end(); it1d_cr4++) {
        it1d_cr4->second->Write();
        delete it1d_cr4->second;
      }

      outfile_cr4.Write();
      outfile_cr4.Close();

      TFile outfile_cr4v(Form("CR4v%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR4v histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr4v;
      for(it1d_cr4v=h_1d_cr4v.begin(); it1d_cr4v!=h_1d_cr4v.end(); it1d_cr4v++) {
        it1d_cr4v->second->Write();
        delete it1d_cr4v->second;
      }

      outfile_cr4v.Write();
      outfile_cr4v.Close();

      TFile outfile_cr40(Form("CR40%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR40 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr40;
      for(it1d_cr40=h_1d_cr40.begin(); it1d_cr40!=h_1d_cr40.end(); it1d_cr40++) {
        it1d_cr40->second->Write();
        delete it1d_cr40->second;
      }

      outfile_cr40.Write();
      outfile_cr40.Close();


      TFile outfile_cr5(Form("CR5%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR5 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr5;
      for(it1d_cr5=h_1d_cr5.begin(); it1d_cr5!=h_1d_cr5.end(); it1d_cr5++) {
        it1d_cr5->second->Write();
        delete it1d_cr5->second;
      }

      outfile_cr5.Write();
      outfile_cr5.Close();

      TFile outfile_cr5v(Form("CR5v%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR5v histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr5v;
      for(it1d_cr5v=h_1d_cr5v.begin(); it1d_cr5v!=h_1d_cr5v.end(); it1d_cr5v++) {
        it1d_cr5v->second->Write();
        delete it1d_cr5v->second;
      }

      outfile_cr5v.Write();
      outfile_cr5v.Close();


      TFile outfile_cr6(Form("CR6%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR6 histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr6;
      for(it1d_cr6=h_1d_cr6.begin(); it1d_cr6!=h_1d_cr6.end(); it1d_cr6++) {
        it1d_cr6->second->Write();
        delete it1d_cr6->second;
      }

      outfile_cr6.Write();
      outfile_cr6.Close();

      TFile outfile_cr6v(Form("CR6v%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving CR6v histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_cr6v;
      for(it1d_cr6v=h_1d_cr6v.begin(); it1d_cr6v!=h_1d_cr6v.end(); it1d_cr6v++) {
        it1d_cr6v->second->Write();
        delete it1d_cr6v->second;
      }

      outfile_cr6v.Write();
      outfile_cr6v.Close();


    TFile outfile_nj(Form("NJ%s", m_outfilename_.c_str()), "RECREATE") ;
    printf("[StopTreeLooper::loop] Saving NJ histograms to %s\n", m_outfilename_.c_str());

    std::map<std::string, TH1D *>::iterator it1d_nj;
    for (it1d_nj = h_1d_nj.begin(); it1d_nj != h_1d_nj.end(); it1d_nj++)
    {
        it1d_nj->second->Write();
        delete it1d_nj->second;
    }

    outfile_nj.Write();
    outfile_nj.Close();

      TFile outfile_z(Form("Z%s",m_outfilename_.c_str()),"RECREATE") ;
      printf("[StopTreeLooper::loop] Saving Z histograms to %s\n", m_outfilename_.c_str());

      std::map<std::string, TH1D*>::iterator it1d_z;
      for(it1d_z=h_1d_z.begin(); it1d_z!=h_1d_z.end(); it1d_z++) {
        it1d_z->second->Write();
        delete it1d_z->second;
      }

      outfile_z.Write();
      outfile_z.Close();


    printf("[StopTreeLooper::loop] nevt[nlep] %i %i %i %i %i %i %i\n", nevt_nlep[0], nevt_nlep[1], nevt_nlep[2], nevt_nlep[3], nevt_nlep[4], nevt_nlep[5], nevt_nlep[6]) ;
    if( name.Contains("ttdl") ) cout<<"Number of events that migrated between channels: "<<nEvents_channel_migrated<<", "<<100.*double(nEvents_channel_migrated)/double(nEventsTotal)<<"%"<<endl;

    already_seen.clear();

    gROOT->cd();

    printf("[StopTreeLooper::loop] Finished %s\n", origname.Data());

}



void StopTreeLooper::makeSIGPlots( float evtweight, std::map<std::string, TH1D *> &h_1d,
                                   string tag_selection, string flav_tag, string hist_tag)
{

    int nbins = 80;

    plot1DUnderOverFlow(hist_tag+"_lep_charge_asymmetry" + tag_selection + flav_tag, lep_charge_asymmetry , evtweight, h_1d, nbins, -4, 4);
    plot1DUnderOverFlow(hist_tag+"_lep_azimuthal_asymmetry" + tag_selection + flav_tag, lep_azimuthal_asymmetry , evtweight, h_1d, nbins, -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_lep_azimuthal_asymmetry2" + tag_selection + flav_tag, lep_azimuthal_asymmetry2 , evtweight, h_1d, nbins, 0, TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_m_top" + tag_selection + flav_tag, m_top , evtweight, h_1d, nbins - 1, 171.4, 173.6);

    if (m_top > 0)
    {
        plot1DUnderOverFlow(hist_tag+"_top_rapiditydiff_cms" + tag_selection + flav_tag, top_rapiditydiff_cms , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow(hist_tag+"_top_pseudorapiditydiff_cms" + tag_selection + flav_tag, top_pseudorapiditydiff_cms , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow(hist_tag+"_top_rapiditydiff_Marco" + tag_selection + flav_tag, top_rapiditydiff_Marco , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow(hist_tag+"_top_costheta_cms" + tag_selection + flav_tag, top_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow(hist_tag+"_lepPlus_costheta_cms" + tag_selection + flav_tag, lepPlus_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow(hist_tag+"_lepMinus_costheta_cms" + tag_selection + flav_tag, lepMinus_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow(hist_tag+"_lep_costheta_cms" + tag_selection + flav_tag, lepPlus_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow(hist_tag+"_lep_costheta_cms" + tag_selection + flav_tag, lepMinus_costheta_cms , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow(hist_tag+"_top_spin_correlation" + tag_selection + flav_tag, top_spin_correlation , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow(hist_tag+"_lep_cos_opening_angle" + tag_selection + flav_tag, lep_cos_opening_angle , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow(hist_tag+"_tt_mass" + tag_selection + flav_tag, tt_mass , evtweight, h_1d, nbins, 0, 1600);
        plot1DUnderOverFlow(hist_tag+"_ttRapidity2" + tag_selection + flav_tag, ttRapidity2 , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow(hist_tag+"_tt_pT" + tag_selection + flav_tag, tt_pT , evtweight, h_1d, nbins, 0, 400);
        plot1DUnderOverFlow(hist_tag+"_top1_pt" + tag_selection + flav_tag, top1_pt , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow(hist_tag+"_top2_pt" + tag_selection + flav_tag, top2_pt , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow(hist_tag+"_top1_p_CM" + tag_selection + flav_tag, top1_p_CM , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow(hist_tag+"_top2_p_CM" + tag_selection + flav_tag, top2_p_CM , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow(hist_tag+"_top_rapiditydiffsigned_cms" + tag_selection + flav_tag, top_rapiditydiffsigned_cms , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow(hist_tag+"_maxAMWTweight" + tag_selection + flav_tag, AMWT_weights.at(imaxweight) , evtweight, h_1d, nbins, 0, 1);
        for (int i = 0; i < int(AMWT_weights.size()); ++i)
        {
            if (i != imaxweight) plot1DUnderOverFlow(hist_tag+"_otherAMWTweights" + tag_selection + flav_tag, AMWT_weights.at(i) , evtweight, h_1d, nbins, 0, 1);
            if (i != imaxweight) plot1DUnderOverFlow(hist_tag+"_otherAMWTweightsratio" + tag_selection + flav_tag, AMWT_weights.at(i) / AMWT_weights.at(imaxweight) , evtweight, h_1d, nbins, 0, 1);
        }
        if (closestApproach)
        {
            plot1DUnderOverFlow(hist_tag+"_closestDeltaMET_bestcombo" + tag_selection + flav_tag, closestDeltaMET_bestcombo , evtweight, h_1d, nbins, 0, 200);
            plot1DUnderOverFlow(hist_tag+"_closestDeltaMET_maxwcombo" + tag_selection + flav_tag, closestDeltaMET_maxwcombo , evtweight, h_1d, nbins, 0, 200);
            if (closestDeltaMET_othercombo > 0) plot1DUnderOverFlow(hist_tag+"_closestDeltaMET_othercombo" + tag_selection + flav_tag, closestDeltaMET_othercombo , evtweight, h_1d, nbins, 0, 200);
            plot1DUnderOverFlow(hist_tag+"_maxAMWTweight_closestApproach" + tag_selection + flav_tag, AMWT_weights.at(imaxweight) , evtweight, h_1d, nbins, 0, 1);
        }
    }

    //channel histograms to be used for yields tables
    plot1DUnderOverFlow(hist_tag+"_channel" + tag_selection + flav_tag, channel, evtweight, h_1d, 4 , 0, 4);
    plot1DUnderOverFlow(hist_tag+"_channel" + tag_selection + flav_tag, 3, evtweight, h_1d, 4 , 0, 4);
    if ( m_top > 0)
    {
        plot1DUnderOverFlow(hist_tag+"_channel_withttbarsol" + tag_selection + flav_tag, channel, evtweight, h_1d, 4 , 0, 4);
        plot1DUnderOverFlow(hist_tag+"_channel_withttbarsol" + tag_selection + flav_tag, 3, evtweight, h_1d, 4 , 0, 4);
    }

    plot1DUnderOverFlow(hist_tag+"_n_jets" + tag_selection + flav_tag, n_jets, evtweight, h_1d, 8 , 0, 8);
    plot1DUnderOverFlow(hist_tag+"_n_bjets" + tag_selection + flav_tag, n_bjets, evtweight, h_1d, 8 , 0, 8);
    plot1DUnderOverFlow(hist_tag+"_nvtx" + tag_selection + flav_tag,  stopt.nvtx(), evtweight, h_1d, 40, 0, 40);

    plot1DUnderOverFlow(hist_tag+"_met" + tag_selection + flav_tag, t1metphicorr, evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_metphi" + tag_selection + flav_tag, t1metphicorrphi, evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_met_smeared" + tag_selection + flav_tag, t1metphicorrsmeared, evtweight, h_1d, nbins , 0, 500);

    plot1DUnderOverFlow(hist_tag+"_lep1b_mindR" + tag_selection + flav_tag, lep1b_mindR, evtweight, h_1d, nbins , 0, 5);
    plot1DUnderOverFlow(hist_tag+"_lep2b_mindR" + tag_selection + flav_tag, lep2b_mindR, evtweight, h_1d, nbins , 0, 5);
    plot1DUnderOverFlow(hist_tag+"_lep1b_mindPhi" + tag_selection + flav_tag, lep1b_mindPhi, evtweight, h_1d, nbins , 0, TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_lep2b_mindPhi" + tag_selection + flav_tag, lep2b_mindPhi, evtweight, h_1d, nbins , 0, TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_lep_dR" + tag_selection + flav_tag, lep_dR, evtweight, h_1d, nbins , 0, 5);
    plot1DUnderOverFlow(hist_tag+"_lep_dEta" + tag_selection + flav_tag, lep_dEta, evtweight, h_1d, nbins , -5, 5);

    plot1DUnderOverFlow(hist_tag+"_Mll" + tag_selection + flav_tag, Mll, evtweight, h_1d, nbins , 0, 800);
    plot1DUnderOverFlow(hist_tag+"_Ptll" + tag_selection + flav_tag, Ptll, evtweight, h_1d, nbins , 0, 400);
    plot1DUnderOverFlow(hist_tag+"_Etall" + tag_selection + flav_tag, Etall, evtweight, h_1d, nbins , -5, 5);
    plot1DUnderOverFlow(hist_tag+"_Phill" + tag_selection + flav_tag, Phill, evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());


    //pT
    plot1DUnderOverFlow(hist_tag+"_lepPlus_Pt" + tag_selection + flav_tag, stopt.lepp().Pt(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_lepMinus_Pt" + tag_selection + flav_tag, stopt.lepm().Pt(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_lepPt" + tag_selection + flav_tag, stopt.lepp().Pt(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_lepPt" + tag_selection + flav_tag, stopt.lepm().Pt(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id1()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepPt_ele" + tag_selection + flav_tag, stopt.lep1().Pt(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id2()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepPt_ele" + tag_selection + flav_tag, stopt.lep2().Pt(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id1()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepPt_muo" + tag_selection + flav_tag, stopt.lep1().Pt(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id2()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepPt_muo" + tag_selection + flav_tag, stopt.lep2().Pt(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet0_Pt" + tag_selection + flav_tag, jet1.Pt(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet1_Pt" + tag_selection + flav_tag, jet2.Pt(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet_Pt" + tag_selection + flav_tag, jet1.Pt(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet_Pt" + tag_selection + flav_tag, jet2.Pt(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b0_Pt" + tag_selection + flav_tag, bjets.at(0).Pt(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b1_Pt" + tag_selection + flav_tag, bjets.at(1).Pt(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b_Pt" + tag_selection + flav_tag, bjets.at(0).Pt(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b_Pt" + tag_selection + flav_tag, bjets.at(1).Pt(), evtweight, h_1d, nbins , 0, 500);
    else if (nonbjets.size()>0) plot1DUnderOverFlow(hist_tag+"_nonb_Pt" + tag_selection + flav_tag, nonbjets.at(0).Pt(), evtweight, h_1d, nbins , 0, 500);

    //eta
    plot1DUnderOverFlow(hist_tag+"_lepPlus_Eta" + tag_selection + flav_tag, stopt.lepp().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    plot1DUnderOverFlow(hist_tag+"_lepMinus_Eta" + tag_selection + flav_tag, stopt.lepm().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    plot1DUnderOverFlow(hist_tag+"_lepEta" + tag_selection + flav_tag, stopt.lepp().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    plot1DUnderOverFlow(hist_tag+"_lepEta" + tag_selection + flav_tag, stopt.lepm().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if ( abs(stopt.id1()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepEta_ele" + tag_selection + flav_tag, stopt.lep1().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if ( abs(stopt.id2()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepEta_ele" + tag_selection + flav_tag, stopt.lep2().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if ( abs(stopt.id1()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepEta_muo" + tag_selection + flav_tag, stopt.lep1().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if ( abs(stopt.id2()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepEta_muo" + tag_selection + flav_tag, stopt.lep2().Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    plot1DUnderOverFlow(hist_tag+"_jet0_Eta" + tag_selection + flav_tag, jet1.Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    plot1DUnderOverFlow(hist_tag+"_jet1_Eta" + tag_selection + flav_tag, jet2.Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    plot1DUnderOverFlow(hist_tag+"_jet_Eta" + tag_selection + flav_tag, jet1.Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    plot1DUnderOverFlow(hist_tag+"_jet_Eta" + tag_selection + flav_tag, jet2.Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b0_Eta" + tag_selection + flav_tag, bjets.at(0).Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b1_Eta" + tag_selection + flav_tag, bjets.at(1).Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b_Eta" + tag_selection + flav_tag, bjets.at(0).Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b_Eta" + tag_selection + flav_tag, bjets.at(1).Eta(), evtweight, h_1d, 52 , -2.6, 2.6);
    else if (nonbjets.size()>0) plot1DUnderOverFlow(hist_tag+"_nonb_Eta" + tag_selection + flav_tag, nonbjets.at(0).Eta(), evtweight, h_1d, 52 , -2.6, 2.6);

    //phi
    plot1DUnderOverFlow(hist_tag+"_lepPlus_Phi" + tag_selection + flav_tag, stopt.lepp().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_lepMinus_Phi" + tag_selection + flav_tag, stopt.lepm().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_lepPhi" + tag_selection + flav_tag, stopt.lepp().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_lepPhi" + tag_selection + flav_tag, stopt.lepm().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if ( abs(stopt.id1()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepPhi_ele" + tag_selection + flav_tag, stopt.lep1().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if ( abs(stopt.id2()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepPhi_ele" + tag_selection + flav_tag, stopt.lep2().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if ( abs(stopt.id1()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepPhi_muo" + tag_selection + flav_tag, stopt.lep1().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if ( abs(stopt.id2()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepPhi_muo" + tag_selection + flav_tag, stopt.lep2().Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_jet0_Phi" + tag_selection + flav_tag, jet1.Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_jet1_Phi" + tag_selection + flav_tag, jet2.Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_jet_Phi" + tag_selection + flav_tag, jet1.Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_jet_Phi" + tag_selection + flav_tag, jet2.Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b0_Phi" + tag_selection + flav_tag, bjets.at(0).Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b1_Phi" + tag_selection + flav_tag, bjets.at(1).Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b_Phi" + tag_selection + flav_tag, bjets.at(0).Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b_Phi" + tag_selection + flav_tag, bjets.at(1).Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    else if (nonbjets.size()>0) plot1DUnderOverFlow(hist_tag+"_nonb_Phi" + tag_selection + flav_tag, nonbjets.at(0).Phi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());

    //E
    plot1DUnderOverFlow(hist_tag+"_lepPlus_E" + tag_selection + flav_tag, stopt.lepp().E(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_lepMinus_E" + tag_selection + flav_tag, stopt.lepm().E(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_lepE" + tag_selection + flav_tag, stopt.lepp().E(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_lepE" + tag_selection + flav_tag, stopt.lepm().E(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id1()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepE_ele" + tag_selection + flav_tag, stopt.lep1().E(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id2()) == 11 ) plot1DUnderOverFlow(hist_tag+"_lepE_ele" + tag_selection + flav_tag, stopt.lep2().E(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id1()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepE_muo" + tag_selection + flav_tag, stopt.lep1().E(), evtweight, h_1d, nbins , 0, 500);
    if ( abs(stopt.id2()) == 13 ) plot1DUnderOverFlow(hist_tag+"_lepE_muo" + tag_selection + flav_tag, stopt.lep2().E(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet0_E" + tag_selection + flav_tag, jet1.E(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet1_E" + tag_selection + flav_tag, jet2.E(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet_E" + tag_selection + flav_tag, jet1.E(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_jet_E" + tag_selection + flav_tag, jet2.E(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b0_E" + tag_selection + flav_tag, bjets.at(0).E(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b1_E" + tag_selection + flav_tag, bjets.at(1).E(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 0) plot1DUnderOverFlow(hist_tag+"_b_E" + tag_selection + flav_tag, bjets.at(0).E(), evtweight, h_1d, nbins , 0, 500);
    if (n_bjets > 1) plot1DUnderOverFlow(hist_tag+"_b_E" + tag_selection + flav_tag, bjets.at(1).E(), evtweight, h_1d, nbins , 0, 500);
    else if (nonbjets.size()>0) plot1DUnderOverFlow(hist_tag+"_nonb_E" + tag_selection + flav_tag, nonbjets.at(0).E(), evtweight, h_1d, nbins , 0, 500);

    //mlb variables
    plot1DUnderOverFlow(hist_tag+"_mlb" + tag_selection + flav_tag, mlb_1 , evtweight, h_1d, nbins, 0., 800.);
    plot1DUnderOverFlow(hist_tag+"_mlb" + tag_selection + flav_tag, mlb_2 , evtweight, h_1d, nbins, 0., 800.);
    plot1DUnderOverFlow(hist_tag+"_mlb" + tag_selection + flav_tag, mlb_3 , evtweight, h_1d, nbins, 0., 800.);
    plot1DUnderOverFlow(hist_tag+"_mlb" + tag_selection + flav_tag, mlb_4 , evtweight, h_1d, nbins, 0., 800.);
    plot1DUnderOverFlow(hist_tag+"_mlb_min" + tag_selection + flav_tag, mlb_min , evtweight, h_1d, nbins, 0., 400.);

    //check HO and TOB/TEC cleanup cut variables
    plot1DUnderOverFlow(hist_tag+"_calomet" + tag_selection + flav_tag, stopt.calomet(), evtweight, h_1d, nbins , 0, 500);
    plot1DUnderOverFlow(hist_tag+"_calometphi" + tag_selection + flav_tag, stopt.calometphi(), evtweight, h_1d, nbins , -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_pfcalo_metratio" + tag_selection + flav_tag, pfcalo_metratio , evtweight, h_1d, nbins, 0, 4.);
    plot1DUnderOverFlow(hist_tag+"_pfcalodPhi" + tag_selection + flav_tag, pfcalo_metdphi , evtweight, h_1d, nbins, 0, TMath::Pi());
    plot1DUnderOverFlow(hist_tag+"_pfcalo_deltamet" + tag_selection + flav_tag, pfcalo_deltamet , evtweight, h_1d, nbins, 0, 200.);

    float pfcalo_deltametx = t1metphicorr * sin(t1metphicorrphi) - stopt.calomet() * sin(stopt.calometphi());
    float pfcalo_deltamety = t1metphicorr * cos(t1metphicorrphi) - stopt.calomet() * cos(stopt.calometphi());

    plot1DUnderOverFlow(hist_tag+"_pfcalo_deltametx" + tag_selection + flav_tag, pfcalo_deltametx , evtweight, h_1d, nbins, -100, 100.);
    plot1DUnderOverFlow(hist_tag+"_pfcalo_deltamety" + tag_selection + flav_tag, pfcalo_deltamety , evtweight, h_1d, nbins, -100, 100.);

}



void StopTreeLooper::makeAccPlots( float evtweight, std::map<std::string, TH1D *> &h_1d, std::map<std::string, TH2D *> &h_2d,
                                   string tag_selection, string flav_tag)
{

    int nbins = 240; //240 bins in the asymmetry variable makes more sense than 80 since it can be divided into 12 and we will be rebinning to 6 or 12 bins for the acceptance histograms

    int nbinsmtt = 120;
    int nbinsttpt = 300;
    int nbinsttrapidity2 = 300;

    double minmtt = 0.0;
    double minttpt = 0.0;
    double minttrapidity2 = 0.0;

    double maxmtt = 1200.0;
    double maxttpt = 300.0;
    double maxttrapidity2 = 3.0;

    //1D distributions
    plot1DUnderOverFlow("h_numerator_lep_charge_asymmetry_gen" + tag_selection + flav_tag, lep_charge_asymmetry_gen , evtweight, h_1d, nbins, -2, 2);
    plot1DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry_gen , evtweight, h_1d, nbins, -TMath::Pi(), TMath::Pi());
    plot1DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry2_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry2_gen , evtweight, h_1d, nbins, 0, TMath::Pi());

    if (m_top > 0)
    {
        plot1DUnderOverFlow("h_numerator_top_rapiditydiff_cms_gen" + tag_selection + flav_tag, top_rapiditydiff_cms_gen , evtweight, h_1d, nbins, -2, 2);
        plot1DUnderOverFlow("h_numerator_top_pseudorapiditydiff_cms_gen" + tag_selection + flav_tag, top_pseudorapiditydiff_cms_gen , evtweight, h_1d, nbins, -2, 2);
        plot1DUnderOverFlow("h_numerator_top_rapiditydiff_Marco_gen" + tag_selection + flav_tag, top_rapiditydiff_Marco_gen , evtweight, h_1d, nbins, -2, 2);
        plot1DUnderOverFlow("h_numerator_top_costheta_cms_gen" + tag_selection + flav_tag, top_costheta_cms_gen , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_numerator_lepPlus_costheta_cms_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_numerator_lepMinus_costheta_cms_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_numerator_lep_costheta_cms_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_numerator_lep_costheta_cms_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_numerator_top_spin_correlation_gen" + tag_selection + flav_tag, top_spin_correlation_gen , evtweight, h_1d, nbins, -1, 1);
        plot1DUnderOverFlow("h_numerator_lep_cos_opening_angle_gen" + tag_selection + flav_tag, lep_cos_opening_angle_gen , evtweight, h_1d, nbins, -1, 1);
    }

    //2D distributions vs Mttbar
    if (m_top > 0)
    {
        plot2DUnderOverFlow("h_numerator_lep_charge_asymmetry_vs_mtt_gen" + tag_selection + flav_tag, lep_charge_asymmetry_gen , tt_mass_gen , evtweight, h_2d, nbins, -2, 2, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry_vs_mtt_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry_gen , tt_mass_gen , evtweight, h_2d, nbins, -TMath::Pi(), TMath::Pi(), nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry2_vs_mtt_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry2_gen , tt_mass_gen , evtweight, h_2d, nbins, 0, TMath::Pi(), nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_top_rapiditydiff_cms_vs_mtt_gen" + tag_selection + flav_tag, top_rapiditydiff_cms_gen , tt_mass_gen , evtweight, h_2d, nbins, -2, 2, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_top_pseudorapiditydiff_cms_vs_mtt_gen" + tag_selection + flav_tag, top_pseudorapiditydiff_cms_gen , tt_mass_gen , evtweight, h_2d, nbins, -2, 2, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_top_rapiditydiff_Marco_vs_mtt_gen" + tag_selection + flav_tag, top_rapiditydiff_Marco_gen , tt_mass_gen , evtweight, h_2d, nbins, -2, 2, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_top_costheta_cms_vs_mtt_gen" + tag_selection + flav_tag, top_costheta_cms_gen , tt_mass_gen , evtweight, h_2d, nbins, -1, 1, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_lepPlus_costheta_cms_vs_mtt_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , tt_mass_gen , evtweight, h_2d, nbins, -1, 1, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_lepMinus_costheta_cms_vs_mtt_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , tt_mass_gen , evtweight, h_2d, nbins, -1, 1, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_lep_costheta_cms_vs_mtt_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , tt_mass_gen , evtweight, h_2d, nbins, -1, 1, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_lep_costheta_cms_vs_mtt_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , tt_mass_gen , evtweight, h_2d, nbins, -1, 1, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_top_spin_correlation_vs_mtt_gen" + tag_selection + flav_tag, top_spin_correlation_gen , tt_mass_gen , evtweight, h_2d, nbins, -1, 1, nbinsmtt, minmtt, maxmtt);
        plot2DUnderOverFlow("h_numerator_lep_cos_opening_angle_vs_mtt_gen" + tag_selection + flav_tag, lep_cos_opening_angle_gen , tt_mass_gen , evtweight, h_2d, nbins, -1, 1, nbinsmtt, minmtt, maxmtt);
    }

    //2D distributions vs ttbar pT
    if (m_top > 0)
    {
        plot2DUnderOverFlow("h_numerator_lep_charge_asymmetry_vs_ttpt_gen" + tag_selection + flav_tag, lep_charge_asymmetry_gen , tt_pT_gen , evtweight, h_2d, nbins, -2, 2, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry_vs_ttpt_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry_gen , tt_pT_gen , evtweight, h_2d, nbins, -TMath::Pi(), TMath::Pi(), nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry2_vs_ttpt_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry2_gen , tt_pT_gen , evtweight, h_2d, nbins, 0, TMath::Pi(), nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_top_rapiditydiff_cms_vs_ttpt_gen" + tag_selection + flav_tag, top_rapiditydiff_cms_gen , tt_pT_gen , evtweight, h_2d, nbins, -2, 2, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_top_pseudorapiditydiff_cms_vs_ttpt_gen" + tag_selection + flav_tag, top_pseudorapiditydiff_cms_gen , tt_pT_gen , evtweight, h_2d, nbins, -2, 2, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_top_rapiditydiff_Marco_vs_ttpt_gen" + tag_selection + flav_tag, top_rapiditydiff_Marco_gen , tt_pT_gen , evtweight, h_2d, nbins, -2, 2, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_top_costheta_cms_vs_ttpt_gen" + tag_selection + flav_tag, top_costheta_cms_gen , tt_pT_gen , evtweight, h_2d, nbins, -1, 1, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_lepPlus_costheta_cms_vs_ttpt_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , tt_pT_gen , evtweight, h_2d, nbins, -1, 1, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_lepMinus_costheta_cms_vs_ttpt_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , tt_pT_gen , evtweight, h_2d, nbins, -1, 1, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_lep_costheta_cms_vs_ttpt_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , tt_pT_gen , evtweight, h_2d, nbins, -1, 1, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_lep_costheta_cms_vs_ttpt_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , tt_pT_gen , evtweight, h_2d, nbins, -1, 1, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_top_spin_correlation_vs_ttpt_gen" + tag_selection + flav_tag, top_spin_correlation_gen , tt_pT_gen , evtweight, h_2d, nbins, -1, 1, nbinsttpt, minttpt, maxttpt);
        plot2DUnderOverFlow("h_numerator_lep_cos_opening_angle_vs_ttpt_gen" + tag_selection + flav_tag, lep_cos_opening_angle_gen , tt_pT_gen , evtweight, h_2d, nbins, -1, 1, nbinsttpt, minttpt, maxttpt);
    }

    //2D distributions vs |ttbar rapidity|
    if (m_top > 0)
    {
        plot2DUnderOverFlow("h_numerator_lep_charge_asymmetry_vs_ttrapidity2_gen" + tag_selection + flav_tag, lep_charge_asymmetry_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -2, 2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry_vs_ttrapidity2_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -TMath::Pi(), TMath::Pi(), nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_lep_azimuthal_asymmetry2_vs_ttrapidity2_gen" + tag_selection + flav_tag, lep_azimuthal_asymmetry2_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, 0, TMath::Pi(), nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_top_rapiditydiff_cms_vs_ttrapidity2_gen" + tag_selection + flav_tag, top_rapiditydiff_cms_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -2, 2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_top_pseudorapiditydiff_cms_vs_ttrapidity2_gen" + tag_selection + flav_tag, top_pseudorapiditydiff_cms_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -2, 2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_top_rapiditydiff_Marco_vs_ttrapidity2_gen" + tag_selection + flav_tag, top_rapiditydiff_Marco_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -2, 2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_top_costheta_cms_vs_ttrapidity2_gen" + tag_selection + flav_tag, top_costheta_cms_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -1, 1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_lepPlus_costheta_cms_vs_ttrapidity2_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -1, 1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_lepMinus_costheta_cms_vs_ttrapidity2_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -1, 1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_lep_costheta_cms_vs_ttrapidity2_gen" + tag_selection + flav_tag, lepPlus_costheta_cms_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -1, 1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_lep_costheta_cms_vs_ttrapidity2_gen" + tag_selection + flav_tag, lepMinus_costheta_cms_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -1, 1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_top_spin_correlation_vs_ttrapidity2_gen" + tag_selection + flav_tag, top_spin_correlation_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -1, 1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
        plot2DUnderOverFlow("h_numerator_lep_cos_opening_angle_vs_ttrapidity2_gen" + tag_selection + flav_tag, lep_cos_opening_angle_gen , fabs(ttRapidity2_gen) , evtweight, h_2d, nbins, -1, 1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
    }

}







void StopTreeLooper::makettPlots( float evtweight, std::map<std::string, TH1D *> &h_1d, std::map<std::string, TH2D *> &h_2d,
                                  string tag_selection, string flav_tag)
{

    int nbins = 80;


    //gen-level distributions

    if (m_top > 0)
    {
        plot1DUnderOverFlow("h_sig_tt_mass_gen" + tag_selection + flav_tag, tt_mass_gen , evtweight, h_1d, nbins, 0, 1600);
        plot1DUnderOverFlow("h_sig_ttRapidity2_gen" + tag_selection + flav_tag, ttRapidity2_gen , evtweight, h_1d, nbins, -4, 4);
        plot1DUnderOverFlow("h_sig_tt_pT_gen" + tag_selection + flav_tag, tt_pT_gen , evtweight, h_1d, nbins, 0, 400);
        plot1DUnderOverFlow("h_sig_top1_pt_gen" + tag_selection + flav_tag, top1_pt_gen , evtweight, h_1d, nbins, 0, 800);
        plot1DUnderOverFlow("h_sig_top2_pt_gen" + tag_selection + flav_tag, top2_pt_gen , evtweight, h_1d, nbins, 0, 800);
    }


    //ttbar solution resolution plots

    plot1DUnderOverFlow("h_sig_lep_charge_asymmetry_Delta" + tag_selection + flav_tag, lep_charge_asymmetry - lep_charge_asymmetry_gen , evtweight, h_1d, nbins, -3.5, 3.5);
    plot2DUnderOverFlow("h_sig_lep_charge_asymmetry_Delta_2D" + tag_selection + flav_tag, lep_charge_asymmetry - lep_charge_asymmetry_gen ,              lep_charge_asymmetry_gen  , evtweight, h_2d, nbins, -3.5, 3.5,    12, bins_lepChargeAsym);
    plot1DUnderOverFlow("h_sig_lep_azimuthal_asymmetry_Delta" + tag_selection + flav_tag, lep_azimuthal_asymmetry - lep_azimuthal_asymmetry_gen , evtweight, h_1d, nbins, -2.*TMath::Pi(), 2.*TMath::Pi());
    plot2DUnderOverFlow("h_sig_lep_azimuthal_asymmetry_Delta_2D" + tag_selection + flav_tag, lep_azimuthal_asymmetry - lep_azimuthal_asymmetry_gen ,              lep_azimuthal_asymmetry_gen  , evtweight, h_2d, nbins, -2.*TMath::Pi(), 2.*TMath::Pi(),    12, bins_lepAzimAsym);
    plot1DUnderOverFlow("h_sig_lep_azimuthal_asymmetry2_Delta" + tag_selection + flav_tag, lep_azimuthal_asymmetry2 - lep_azimuthal_asymmetry2_gen , evtweight, h_1d, nbins, -TMath::Pi(), TMath::Pi());
    plot2DUnderOverFlow("h_sig_lep_azimuthal_asymmetry2_Delta_2D" + tag_selection + flav_tag, lep_azimuthal_asymmetry2 - lep_azimuthal_asymmetry2_gen ,              lep_azimuthal_asymmetry2_gen  , evtweight, h_2d, nbins, -TMath::Pi(), TMath::Pi(),    12, bins_lepAzimAsym2);

    if (m_top > 0)
    {
        plot1DUnderOverFlow("h_sig_top_rapiditydiff_cms_Delta" + tag_selection + flav_tag, top_rapiditydiff_cms - top_rapiditydiff_cms_gen , evtweight, h_1d, nbins, -3.5, 3.5);
        plot2DUnderOverFlow("h_sig_top_rapiditydiff_cms_Delta_2D" + tag_selection + flav_tag, top_rapiditydiff_cms - top_rapiditydiff_cms_gen ,              top_rapiditydiff_cms_gen  , evtweight, h_2d, nbins, -3.5, 3.5,    6, bins_rapiditydiff);
        plot1DUnderOverFlow("h_sig_top_pseudorapiditydiff_cms_Delta" + tag_selection + flav_tag, top_pseudorapiditydiff_cms - top_pseudorapiditydiff_cms_gen , evtweight, h_1d, nbins, -3.5, 3.5);
        plot2DUnderOverFlow("h_sig_top_pseudorapiditydiff_cms_Delta_2D" + tag_selection + flav_tag, top_pseudorapiditydiff_cms - top_pseudorapiditydiff_cms_gen ,              top_pseudorapiditydiff_cms_gen  , evtweight, h_2d, nbins, -3.5, 3.5,    6, bins_pseudorapiditydiff);
        plot1DUnderOverFlow("h_sig_top_rapiditydiff_Marco_Delta" + tag_selection + flav_tag, top_rapiditydiff_Marco - top_rapiditydiff_Marco_gen , evtweight, h_1d, nbins, -3.5, 3.5);
        plot2DUnderOverFlow("h_sig_top_rapiditydiff_Marco_Delta_2D" + tag_selection + flav_tag, top_rapiditydiff_Marco - top_rapiditydiff_Marco_gen ,              top_rapiditydiff_Marco_gen  , evtweight, h_2d, nbins, -3.5, 3.5,    6, bins_rapiditydiffMarco);
        plot1DUnderOverFlow("h_sig_top_costheta_cms_Delta" + tag_selection + flav_tag, top_costheta_cms - top_costheta_cms_gen , evtweight, h_1d, nbins, -2, 2);
        plot2DUnderOverFlow("h_sig_top_costheta_cms_Delta_2D" + tag_selection + flav_tag, top_costheta_cms - top_costheta_cms_gen ,              top_costheta_cms_gen  , evtweight, h_2d, nbins, -2, 2,    6, bins_topCosTheta);
        plot1DUnderOverFlow("h_sig_lepPlus_costheta_cms_Delta" + tag_selection + flav_tag, lepPlus_costheta_cms - lepPlus_costheta_cms_gen , evtweight, h_1d, nbins, -2, 2);
        plot2DUnderOverFlow("h_sig_lepPlus_costheta_cms_Delta_2D" + tag_selection + flav_tag, lepPlus_costheta_cms - lepPlus_costheta_cms_gen ,              lepPlus_costheta_cms_gen  , evtweight, h_2d, nbins, -2, 2,    6, bins_lepCosTheta);
        plot1DUnderOverFlow("h_sig_lepMinus_costheta_cms_Delta" + tag_selection + flav_tag, lepMinus_costheta_cms - lepMinus_costheta_cms_gen , evtweight, h_1d, nbins, -2, 2);
        plot2DUnderOverFlow("h_sig_lepMinus_costheta_cms_Delta_2D" + tag_selection + flav_tag, lepMinus_costheta_cms - lepMinus_costheta_cms_gen ,              lepMinus_costheta_cms_gen  , evtweight, h_2d, nbins, -2, 2,    6, bins_lepCosTheta);
        plot1DUnderOverFlow("h_sig_lep_costheta_cms_Delta" + tag_selection + flav_tag, lepPlus_costheta_cms - lepPlus_costheta_cms_gen , evtweight, h_1d, nbins, -2, 2);
        plot2DUnderOverFlow("h_sig_lep_costheta_cms_Delta_2D" + tag_selection + flav_tag, lepPlus_costheta_cms - lepPlus_costheta_cms_gen ,              lepPlus_costheta_cms_gen  , evtweight, h_2d, nbins, -2, 2,    6, bins_lepCosTheta);
        plot1DUnderOverFlow("h_sig_lep_costheta_cms_Delta" + tag_selection + flav_tag, lepMinus_costheta_cms - lepMinus_costheta_cms_gen , evtweight, h_1d, nbins, -2, 2);
        plot2DUnderOverFlow("h_sig_lep_costheta_cms_Delta_2D" + tag_selection + flav_tag, lepMinus_costheta_cms - lepMinus_costheta_cms_gen ,              lepMinus_costheta_cms_gen  , evtweight, h_2d, nbins, -2, 2,    6, bins_lepCosTheta);
        plot1DUnderOverFlow("h_sig_top_spin_correlation_Delta" + tag_selection + flav_tag, top_spin_correlation - top_spin_correlation_gen , evtweight, h_1d, nbins, -2, 2);
        plot2DUnderOverFlow("h_sig_top_spin_correlation_Delta_2D" + tag_selection + flav_tag, top_spin_correlation - top_spin_correlation_gen ,              top_spin_correlation_gen  , evtweight, h_2d, nbins, -2, 2,    6, bins_topSpinCorr);
        plot1DUnderOverFlow("h_sig_lep_cos_opening_angle_Delta" + tag_selection + flav_tag, lep_cos_opening_angle - lep_cos_opening_angle_gen , evtweight, h_1d, nbins, -2, 2);
        plot2DUnderOverFlow("h_sig_lep_cos_opening_angle_Delta_2D" + tag_selection + flav_tag, lep_cos_opening_angle - lep_cos_opening_angle_gen ,              lep_cos_opening_angle_gen  , evtweight, h_2d, nbins, -2, 2,    6, bins_lepCosOpeningAngle);

        plot1DUnderOverFlow("h_sig_tt_mass_Delta" + tag_selection + flav_tag, tt_mass - tt_mass_gen , evtweight, h_1d, nbins, -500, 500);
        plot2DUnderOverFlow("h_sig_tt_mass_Delta_2D" + tag_selection + flav_tag, tt_mass - tt_mass_gen ,              tt_mass_gen  , evtweight, h_2d, nbins, -500, 500,    3, ybinsmtt);
        plot1DUnderOverFlow("h_sig_ttRapidity2_Delta" + tag_selection + flav_tag, ttRapidity2 - ttRapidity2_gen , evtweight, h_1d, nbins, -3.5, 3.5);
        plot2DUnderOverFlow("h_sig_ttRapidity2_Delta_2D" + tag_selection + flav_tag, ttRapidity2 - ttRapidity2_gen ,              fabs(ttRapidity2_gen)  , evtweight, h_2d, nbins, -3.5, 3.5,    3, ybinsttrapidity2);
        plot1DUnderOverFlow("h_sig_absttRapidity2_Delta" + tag_selection + flav_tag, fabs(ttRapidity2) - fabs(ttRapidity2_gen) , evtweight, h_1d, nbins, -3.5, 3.5);
        plot2DUnderOverFlow("h_sig_absttRapidity2_Delta_2D" + tag_selection + flav_tag, fabs(ttRapidity2) - fabs(ttRapidity2_gen) ,              fabs(ttRapidity2_gen)  , evtweight, h_2d, nbins, -3.5, 3.5,    3, ybinsttrapidity2);
        plot1DUnderOverFlow("h_sig_tt_pT_Delta" + tag_selection + flav_tag, tt_pT - tt_pT_gen , evtweight, h_1d, nbins, -200, 200);
        plot2DUnderOverFlow("h_sig_tt_pT_Delta_2D" + tag_selection + flav_tag, tt_pT - tt_pT_gen ,              tt_pT_gen  , evtweight, h_2d, nbins, -200, 200,    3, ybinsttpt);
        plot1DUnderOverFlow("h_sig_top1_pt_Delta" + tag_selection + flav_tag, top1_pt - top1_pt_gen , evtweight, h_1d, nbins, -300, 300);
        plot2DUnderOverFlow("h_sig_top1_pt_Delta_2D" + tag_selection + flav_tag, top1_pt - top1_pt_gen ,              top1_pt_gen  , evtweight, h_2d, nbins, -300, 300,    3, 0, 300);
        plot1DUnderOverFlow("h_sig_top2_pt_Delta" + tag_selection + flav_tag, top2_pt - top2_pt_gen , evtweight, h_1d, nbins, -300, 300);
        plot2DUnderOverFlow("h_sig_top2_pt_Delta_2D" + tag_selection + flav_tag, top2_pt - top2_pt_gen ,              top2_pt_gen  , evtweight, h_2d, nbins, -300, 300,    3, 0, 300);
    }

    if (m_top > 0)
    {
        double top1dotgen = top1_vecs.at(imaxweight).Vect().Dot( topplus_genp_p4.Vect() ) / top1_vecs.at(imaxweight).Vect().Mag() / topplus_genp_p4.Vect().Mag();
        double top1dotgent2 = top1_vecs.at(imaxweight).Vect().Dot( topminus_genp_p4.Vect() ) / top1_vecs.at(imaxweight).Vect().Mag() / topminus_genp_p4.Vect().Mag();
        double top2dotgent1 = top2_vecs.at(imaxweight).Vect().Dot( topplus_genp_p4.Vect() ) / top2_vecs.at(imaxweight).Vect().Mag() / topplus_genp_p4.Vect().Mag();
        double top2dotgen = top2_vecs.at(imaxweight).Vect().Dot( topminus_genp_p4.Vect() ) / top2_vecs.at(imaxweight).Vect().Mag() / topminus_genp_p4.Vect().Mag();
        double top1Pratio = ( top1_vecs.at(imaxweight).Vect().Mag() - topplus_genp_p4.Vect().Mag() ) / ( top1_vecs.at(imaxweight).Vect().Mag() + topplus_genp_p4.Vect().Mag() );
        double top2Pratio = ( top2_vecs.at(imaxweight).Vect().Mag() - topminus_genp_p4.Vect().Mag() ) / ( top2_vecs.at(imaxweight).Vect().Mag() + topminus_genp_p4.Vect().Mag() );

        /*
            double nu1dotgen = nu1_vecs.at(imaxweight).Vect().Dot( nuPlus_gen.Vect() ) / nu1_vecs.at(imaxweight).Vect().Mag() / nuPlus_gen.Vect().Mag();
            double nu1dotgennu2 = nu1_vecs.at(imaxweight).Vect().Dot( nuMinus_gen.Vect() ) / nu1_vecs.at(imaxweight).Vect().Mag() / nuMinus_gen.Vect().Mag();
            double nu2dotgennu1 = nu2_vecs.at(imaxweight).Vect().Dot( nuPlus_gen.Vect() ) / nu2_vecs.at(imaxweight).Vect().Mag() / nuPlus_gen.Vect().Mag();
            double nu2dotgen = nu2_vecs.at(imaxweight).Vect().Dot( nuMinus_gen.Vect() ) / nu2_vecs.at(imaxweight).Vect().Mag() / nuMinus_gen.Vect().Mag();
            double nu1Pratio = ( nu1_vecs.at(imaxweight).Vect().Mag() - nuPlus_gen.Vect().Mag() ) / ( nu1_vecs.at(imaxweight).Vect().Mag() + nuPlus_gen.Vect().Mag() );
            double nu2Pratio = ( nu2_vecs.at(imaxweight).Vect().Mag() - nuMinus_gen.Vect().Mag() ) / ( nu2_vecs.at(imaxweight).Vect().Mag() + nuMinus_gen.Vect().Mag() );


            TLorentzVector nusum_gen = nuPlus_gen + nuMinus_gen;
        */

        //double met_x = t1metphicorr * cos(t1metphicorrphi);
        //double met_y = t1metphicorr * sin(t1metphicorrphi);

        //double DeltaMETsol_gen = sqrt( pow( nusum.Px() - nusum_gen.Px() , 2 ) + pow( nusum.Py() - nusum_gen.Py() , 2 ) );
        //double DeltaMETmeas_gen = sqrt( pow( met_x - nusum_gen.Px() , 2 ) + pow( met_y - nusum_gen.Py() , 2 ) );

        //cout<<"about to plot2D"<<endl;



        plot2DUnderOverFlow("h_tt_topdotgen_vs_MET" + tag_selection + flav_tag, acos(top1dotgen), t1metphicorr, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
        plot2DUnderOverFlow("h_tt_topdotgen_vs_MET" + tag_selection + flav_tag, acos(top2dotgen), t1metphicorr, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);

        plot2DUnderOverFlow("h_tt_topdotgen_vs_closestDeltaMET" + tag_selection + flav_tag, acos(top1dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);
        plot2DUnderOverFlow("h_tt_topdotgen_vs_closestDeltaMET" + tag_selection + flav_tag, acos(top2dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);

        plot2DUnderOverFlow("h_tt_topdotgen_vs_njets" + tag_selection + flav_tag, acos(top1dotgen), n_jets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
        plot2DUnderOverFlow("h_tt_topdotgen_vs_njets" + tag_selection + flav_tag, acos(top2dotgen), n_jets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
        plot2DUnderOverFlow("h_tt_topdotgen_vs_nbtag" + tag_selection + flav_tag, acos(top1dotgen), n_bjets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);
        plot2DUnderOverFlow("h_tt_topdotgen_vs_nbtag" + tag_selection + flav_tag, acos(top2dotgen), n_bjets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);

        //plot2DUnderOverFlow("h_tt_topdotgen_vs_ntaugen" + tag_selection + flav_tag, acos(top1dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);
        //plot2DUnderOverFlow("h_tt_topdotgen_vs_ntaugen" + tag_selection + flav_tag, acos(top2dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);

        plot2DUnderOverFlow("h_tt_topdotgen_vs_mttgen" + tag_selection + flav_tag, acos(top1dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);
        plot2DUnderOverFlow("h_tt_topdotgen_vs_mttgen" + tag_selection + flav_tag, acos(top2dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);



        plot2DUnderOverFlow("h_tt_topPratio_vs_MET" + tag_selection + flav_tag, fabs(top1Pratio), t1metphicorr, 1., h_2d, 40, 0., 1., 40, 0., 200);
        plot2DUnderOverFlow("h_tt_topPratio_vs_MET" + tag_selection + flav_tag, fabs(top2Pratio), t1metphicorr, 1., h_2d, 40, 0., 1., 40, 0., 200);

        plot2DUnderOverFlow("h_tt_topPratio_vs_closestDeltaMET" + tag_selection + flav_tag, fabs(top1Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);
        plot2DUnderOverFlow("h_tt_topPratio_vs_closestDeltaMET" + tag_selection + flav_tag, fabs(top2Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);

        plot2DUnderOverFlow("h_tt_topPratio_vs_njets" + tag_selection + flav_tag, fabs(top1Pratio), n_jets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
        plot2DUnderOverFlow("h_tt_topPratio_vs_njets" + tag_selection + flav_tag, fabs(top2Pratio), n_jets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
        plot2DUnderOverFlow("h_tt_topPratio_vs_nbtag" + tag_selection + flav_tag, fabs(top1Pratio), n_bjets, 1., h_2d, 40, 0., 1., 4, 0., 4.);
        plot2DUnderOverFlow("h_tt_topPratio_vs_nbtag" + tag_selection + flav_tag, fabs(top2Pratio), n_bjets, 1., h_2d, 40, 0., 1., 4, 0., 4.);

        //plot2DUnderOverFlow("h_tt_topPratio_vs_ntaugen" + tag_selection + flav_tag, fabs(top1Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);
        //plot2DUnderOverFlow("h_tt_topPratio_vs_ntaugen" + tag_selection + flav_tag, fabs(top2Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);

        plot2DUnderOverFlow("h_tt_topPratio_vs_mttgen" + tag_selection + flav_tag, fabs(top1Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);
        plot2DUnderOverFlow("h_tt_topPratio_vs_mttgen" + tag_selection + flav_tag, fabs(top2Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);


        /*


            plot2DUnderOverFlow("h_tt_nudotgen_vs_MET" + tag_selection + flav_tag, acos(nu1dotgen), t1metphicorr, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);
            plot2DUnderOverFlow("h_tt_nudotgen_vs_MET" + tag_selection + flav_tag, acos(nu2dotgen), t1metphicorr, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 200);

            plot2DUnderOverFlow("h_tt_nudotgen_vs_closestDeltaMET" + tag_selection + flav_tag, acos(nu1dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);
            plot2DUnderOverFlow("h_tt_nudotgen_vs_closestDeltaMET" + tag_selection + flav_tag, acos(nu2dotgen), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 3.1415926536, 40, 0., 100);

            plot2DUnderOverFlow("h_tt_nudotgen_vs_njets" + tag_selection + flav_tag, acos(nu1dotgen), n_jets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
            plot2DUnderOverFlow("h_tt_nudotgen_vs_njets" + tag_selection + flav_tag, acos(nu2dotgen), n_jets, 1., h_2d, 40, 0., 3.1415926536, 8, 0., 8.);
            plot2DUnderOverFlow("h_tt_nudotgen_vs_nbtag" + tag_selection + flav_tag, acos(nu1dotgen), n_bjets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);
            plot2DUnderOverFlow("h_tt_nudotgen_vs_nbtag" + tag_selection + flav_tag, acos(nu2dotgen), n_bjets, 1., h_2d, 40, 0., 3.1415926536, 4, 0., 4.);

            //plot2DUnderOverFlow("h_tt_nudotgen_vs_ntaugen" + tag_selection + flav_tag, acos(nu1dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);
            //plot2DUnderOverFlow("h_tt_nudotgen_vs_ntaugen" + tag_selection + flav_tag, acos(nu2dotgen), ntaus, 1., h_2d, 40, 0., 3.1415926536, 3, 0., 3.);

            plot2DUnderOverFlow("h_tt_nudotgen_vs_mttgen" + tag_selection + flav_tag, acos(nu1dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);
            plot2DUnderOverFlow("h_tt_nudotgen_vs_mttgen" + tag_selection + flav_tag, acos(nu2dotgen), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 3.1415926536, 40, 285., 1485.);



            plot2DUnderOverFlow("h_tt_nuPratio_vs_MET" + tag_selection + flav_tag, fabs(nu1Pratio), t1metphicorr, 1., h_2d, 40, 0., 1., 40, 0., 200);
            plot2DUnderOverFlow("h_tt_nuPratio_vs_MET" + tag_selection + flav_tag, fabs(nu2Pratio), t1metphicorr, 1., h_2d, 40, 0., 1., 40, 0., 200);

            plot2DUnderOverFlow("h_tt_nuPratio_vs_closestDeltaMET" + tag_selection + flav_tag, fabs(nu1Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);
            plot2DUnderOverFlow("h_tt_nuPratio_vs_closestDeltaMET" + tag_selection + flav_tag, fabs(nu2Pratio), closestDeltaMET_bestcombo, 1., h_2d, 40, 0., 1., 40, 0., 100);

            plot2DUnderOverFlow("h_tt_nuPratio_vs_njets" + tag_selection + flav_tag, fabs(nu1Pratio), n_jets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
            plot2DUnderOverFlow("h_tt_nuPratio_vs_njets" + tag_selection + flav_tag, fabs(nu2Pratio), n_jets, 1., h_2d, 40, 0., 1., 8, 0., 8.);
            plot2DUnderOverFlow("h_tt_nuPratio_vs_nbtag" + tag_selection + flav_tag, fabs(nu1Pratio), n_bjets, 1., h_2d, 40, 0., 1., 4, 0., 4.);
            plot2DUnderOverFlow("h_tt_nuPratio_vs_nbtag" + tag_selection + flav_tag, fabs(nu2Pratio), n_bjets, 1., h_2d, 40, 0., 1., 4, 0., 4.);

            //plot2DUnderOverFlow("h_tt_nuPratio_vs_ntaugen" + tag_selection + flav_tag, fabs(nu1Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);
            //plot2DUnderOverFlow("h_tt_nuPratio_vs_ntaugen" + tag_selection + flav_tag, fabs(nu2Pratio), ntaus, 1., h_2d, 40, 0., 1., 3, 0., 3.);

            plot2DUnderOverFlow("h_tt_nuPratio_vs_mttgen" + tag_selection + flav_tag, fabs(nu1Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);
            plot2DUnderOverFlow("h_tt_nuPratio_vs_mttgen" + tag_selection + flav_tag, fabs(nu2Pratio), (topplus_genp_p4 + topminus_genp_p4).M(), 1., h_2d, 40, 0., 1., 40, 285., 1485.);
        */




        //if (top1sdp > -998) plot1DUnderOverFlow("h_tt_top1sdp" + tag_selection + flav_tag, acos(top1sdp), 1., h_1d, 40, 0., 3.1415926536);
        //if (top2sdp > -998) plot1DUnderOverFlow("h_tt_top2sdp" + tag_selection + flav_tag, acos(top2sdp), 1., h_1d, 40, 0., 3.1415926536);

        plot1DUnderOverFlow("h_tt_top1dotgen" + tag_selection + flav_tag, acos(top1dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow("h_tt_top2dotgen" + tag_selection + flav_tag, acos(top2dotgen), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow("h_tt_top2dotgent1" + tag_selection + flav_tag, acos(top2dotgent1), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow("h_tt_top1dotgent2" + tag_selection + flav_tag, acos(top1dotgent2), 1., h_1d, 40, 0., 3.1415926536);
        plot1DUnderOverFlow("h_tt_top1Pratio" + tag_selection + flav_tag, top1Pratio, 1., h_1d, 40, -1., 1.);
        plot1DUnderOverFlow("h_tt_top2Pratio" + tag_selection + flav_tag, top2Pratio, 1., h_1d, 40, -1., 1.);
        /*
            plot1DUnderOverFlow("h_tt_nu1dotgen" + tag_selection + flav_tag, acos(nu1dotgen), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_nu2dotgen" + tag_selection + flav_tag, acos(nu2dotgen), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_nu2dotgennu1" + tag_selection + flav_tag, acos(nu2dotgennu1), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_nu1dotgennu2" + tag_selection + flav_tag, acos(nu1dotgennu2), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_nu1Pratio" + tag_selection + flav_tag, nu1Pratio, 1., h_1d, 40, -1., 1.);
            plot1DUnderOverFlow("h_tt_nu2Pratio" + tag_selection + flav_tag, nu2Pratio, 1., h_1d, 40, -1., 1.);
        */
        //plot1DUnderOverFlow("h_tt_DeltaMETsol_gen" + tag_selection + flav_tag,  DeltaMETsol_gen , 1., h_1d, 40, 0., 200.);
        //plot1DUnderOverFlow("h_tt_DeltaMETmeas_gen" + tag_selection + flav_tag,  DeltaMETmeas_gen , 1., h_1d, 40, 0., 200.);

        if (closestApproach)
        {
            plot1DUnderOverFlow("h_tt_top1dotgen_closest_bestsol" + tag_selection + flav_tag, acos(top1dotgen), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_top2dotgen_closest_bestsol" + tag_selection + flav_tag, acos(top2dotgen), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_top1Pratio_closest_bestsol" + tag_selection + flav_tag, top1Pratio, 1., h_1d, 40, -1., 1.);
            plot1DUnderOverFlow("h_tt_top2Pratio_closest_bestsol" + tag_selection + flav_tag, top2Pratio, 1., h_1d, 40, -1., 1.);
            /*
                    plot1DUnderOverFlow("h_tt_nu1dotgen_closest_bestsol" + tag_selection + flav_tag, acos(nu1dotgen), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_nu2dotgen_closest_bestsol" + tag_selection + flav_tag, acos(nu2dotgen), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_nu1Pratio_closest_bestsol" + tag_selection + flav_tag, nu1Pratio, 1., h_1d, 40, -1., 1.);
                    plot1DUnderOverFlow("h_tt_nu2Pratio_closest_bestsol" + tag_selection + flav_tag, nu2Pratio, 1., h_1d, 40, -1., 1.);
            */
            plot1DUnderOverFlow("h_tt_closestDeltaMET_closest_maxwsol" + tag_selection + flav_tag,      closestDeltaMET_maxwcombo , 1., h_1d, 40, 0., 200.);
            plot1DUnderOverFlow("h_tt_closestDeltaMET_closest_bestsol" + tag_selection + flav_tag,      closestDeltaMET_bestcombo , 1., h_1d, 40, 0., 200.);
            if (closestDeltaMET_othercombo > 0) plot1DUnderOverFlow("h_tt_closestDeltaMET_closest_othersol" + tag_selection + flav_tag,      closestDeltaMET_othercombo , 1., h_1d, 40, 0., 200.);

            //plot1DUnderOverFlow("h_tt_DeltaMETsol_gen_closest" + tag_selection + flav_tag,  DeltaMETsol_gen , 1., h_1d, 40, 0., 200.);
            //plot1DUnderOverFlow("h_tt_DeltaMETmeas_gen_closest" + tag_selection + flav_tag,  DeltaMETmeas_gen , 1., h_1d, 40, 0., 200.);

            for (int i = 0; i < int(top1_vecs.size()); ++i)
            {
                if (i != imaxweight)
                {
                    double top1dotgen_i = top1_vecs.at(i).Vect().Dot( topplus_genp_p4.Vect() ) / top1_vecs.at(i).Vect().Mag() / topplus_genp_p4.Vect().Mag();
                    double top2dotgen_i = top2_vecs.at(i).Vect().Dot( topminus_genp_p4.Vect() ) / top2_vecs.at(i).Vect().Mag() / topminus_genp_p4.Vect().Mag();
                    double top1Pratio_i = ( top1_vecs.at(i).Vect().Mag() - topplus_genp_p4.Vect().Mag() ) / ( top1_vecs.at(i).Vect().Mag() + topplus_genp_p4.Vect().Mag() );
                    double top2Pratio_i = ( top2_vecs.at(i).Vect().Mag() - topminus_genp_p4.Vect().Mag() ) / ( top2_vecs.at(i).Vect().Mag() + topminus_genp_p4.Vect().Mag() );
                    plot1DUnderOverFlow("h_tt_top1dotgen_closest_othersol" + tag_selection + flav_tag, acos(top1dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_top2dotgen_closest_othersol" + tag_selection + flav_tag, acos(top2dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_top1Pratio_closest_othersol" + tag_selection + flav_tag, top1Pratio_i, 1., h_1d, 40, -1., 1.);
                    plot1DUnderOverFlow("h_tt_top2Pratio_closest_othersol" + tag_selection + flav_tag, top2Pratio_i, 1., h_1d, 40, -1., 1.);
                }
            }


        }
        else
        {
            plot1DUnderOverFlow("h_tt_top1dotgen_max" + tag_selection + flav_tag, acos(top1dotgen), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_top2dotgen_max" + tag_selection + flav_tag, acos(top2dotgen), 1., h_1d, 40, 0., 3.1415926536);
            plot1DUnderOverFlow("h_tt_top1Pratio_max" + tag_selection + flav_tag, top1Pratio, 1., h_1d, 40, -1., 1.);
            plot1DUnderOverFlow("h_tt_top2Pratio_max" + tag_selection + flav_tag, top2Pratio, 1., h_1d, 40, -1., 1.);
            /*
                    plot1DUnderOverFlow("h_tt_nu1dotgen_max" + tag_selection + flav_tag, acos(nu1dotgen), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_nu2dotgen_max" + tag_selection + flav_tag, acos(nu2dotgen), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_nu1Pratio_max" + tag_selection + flav_tag, nu1Pratio, 1., h_1d, 40, -1., 1.);
                    plot1DUnderOverFlow("h_tt_nu2Pratio_max" + tag_selection + flav_tag, nu2Pratio, 1., h_1d, 40, -1., 1.);
            */
            plot1DUnderOverFlow("h_tt_closestDeltaMET_max" + tag_selection + flav_tag,      closestDeltaMET_maxwcombo , 1., h_1d, 40, 0., 200.);

            //plot1DUnderOverFlow("h_tt_DeltaMETsol_gen_max" + tag_selection + flav_tag,  DeltaMETsol_gen , 1., h_1d, 40, 0., 200.);
            //plot1DUnderOverFlow("h_tt_DeltaMETmeas_gen_max" + tag_selection + flav_tag,  DeltaMETmeas_gen , 1., h_1d, 40, 0., 200.);

            for (int i = 0; i < int(top1_vecs.size()); ++i)
            {
                if (i != imaxweight)
                {
                    double top1dotgen_i = top1_vecs.at(i).Vect().Dot( topplus_genp_p4.Vect() ) / top1_vecs.at(i).Vect().Mag() / topplus_genp_p4.Vect().Mag();
                    double top2dotgen_i = top2_vecs.at(i).Vect().Dot( topminus_genp_p4.Vect() ) / top2_vecs.at(i).Vect().Mag() / topminus_genp_p4.Vect().Mag();
                    double top1Pratio_i = ( top1_vecs.at(i).Vect().Mag() - topplus_genp_p4.Vect().Mag() ) / ( top1_vecs.at(i).Vect().Mag() + topplus_genp_p4.Vect().Mag() );
                    double top2Pratio_i = ( top2_vecs.at(i).Vect().Mag() - topminus_genp_p4.Vect().Mag() ) / ( top2_vecs.at(i).Vect().Mag() + topminus_genp_p4.Vect().Mag() );
                    plot1DUnderOverFlow("h_tt_top1dotgen_othersols" + tag_selection + flav_tag, acos(top1dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_top2dotgen_othersols" + tag_selection + flav_tag, acos(top2dotgen_i), 1., h_1d, 40, 0., 3.1415926536);
                    plot1DUnderOverFlow("h_tt_top1Pratio_othersols" + tag_selection + flav_tag, top1Pratio_i, 1., h_1d, 40, -1., 1.);
                    plot1DUnderOverFlow("h_tt_top2Pratio_othersols" + tag_selection + flav_tag, top2Pratio_i, 1., h_1d, 40, -1., 1.);
                }
            }



        }
    }// m_top>0


}









//template function, not used
void StopTreeLooper::makeCRPlots( float evtweight, std::map<std::string, TH1D*> &h_1d, 
                   string tag_selection, string flav_tag, string hist_tag ) 
{
  int nbins = 50;
  float h_xmin = 0.;
  float h_xmax = 500.;

  plot1DUnderOverFlow(hist_tag+"_njets"    +tag_selection+flav_tag, n_jets,  evtweight, h_1d, 5,0,5);
  plot1DUnderOverFlow(hist_tag+"_njets_all"+tag_selection+flav_tag, n_jets,  evtweight, h_1d, 10, 0, 10);
  //default met
  plot1DUnderOverFlow(hist_tag+"_met"+tag_selection+flav_tag, t1metphicorr, evtweight, h_1d, nbins-5, 50., h_xmax);
  //lepton pt - enters mT calculation
  plot1DUnderOverFlow(hist_tag+"_leppt"+tag_selection+flav_tag, stopt.lep1().Pt(), evtweight, h_1d, nbins, h_xmin, h_xmax);
  //lepton phi
  plot1DUnderOverFlow(hist_tag+"_lepphi"+tag_selection+flav_tag, stopt.lep1().Phi(), evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());
  //met phi
  plot1DUnderOverFlow(hist_tag+"_metphi"+tag_selection+flav_tag, t1metphicorrphi, evtweight, h_1d, 30, -1.*TMath::Pi(), TMath::Pi());

  //check HO and TOB/TEC cleanup cut variables
  plot1DUnderOverFlow(hist_tag+"_pfcaloMET"+tag_selection+flav_tag, pfcalo_metratio , evtweight, h_1d, 100, 0, 4.);
  plot1DUnderOverFlow(hist_tag+"_pfcalodPhi"+tag_selection+flav_tag, pfcalo_metdphi , evtweight, h_1d, 100, 0, TMath::Pi());  
}







void StopTreeLooper::makeNJPlots( float evtweight, std::map<std::string, TH1D *> &h_1d,
                                  string tag_selection, string flav_tag )
{

    plot1DUnderOverFlow("h_njets"    + tag_selection,          min(n_jets, 4), evtweight, h_1d, 4, 1, 5);
    plot1DUnderOverFlow("h_njets"    + tag_selection + flav_tag, min(n_jets, 4), evtweight, h_1d, 4, 1, 5);
    plot1DUnderOverFlow("h_njets_all" + tag_selection,          min(n_jets, 8), evtweight, h_1d, 7, 1, 8);
    plot1DUnderOverFlow("h_njets_all" + tag_selection + flav_tag, min(n_jets, 8), evtweight, h_1d, 7, 1, 8);

}

void StopTreeLooper::makeZPlots( float evtweight, std::map<std::string, TH1D *> &h_1d,
                                 string tag_selection, string flav_tag )
{

    int nbins = 30;
    float h_xmin = 0.;
    float h_xmax = 300.;

    plot1DUnderOverFlow("h_z_met" + tag_selection + flav_tag, t1metphicorr, evtweight, h_1d, nbins, h_xmin, h_xmax);

    string lep1type =  abs(stopt.id1()) == 13 ? "h_muo" : "h_ele";
    string lep2type =  abs(stopt.id2()) == 13 ? "h_muo" : "h_ele";
    plot1DUnderOverFlow(lep1type + "pt" + tag_selection + flav_tag, min(stopt.lep1().Pt(), (float)199.99), evtweight, h_1d, 40, 20, 200);
    plot1DUnderOverFlow(lep2type + "pt" + tag_selection + flav_tag, min(stopt.lep2().Pt(), (float)199.99), evtweight, h_1d, 40, 20, 200);

    plot1DUnderOverFlow("h_z_leppt"  + tag_selection + flav_tag, min(stopt.lep1().Pt(), (float)299.99), evtweight, h_1d, 50, 20., 300.);
    plot1DUnderOverFlow("h_z_lepeta" + tag_selection + flav_tag, stopt.lep1().Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    plot1DUnderOverFlow("h_z_lep2pt" + tag_selection + flav_tag, min(stopt.lep2().Pt(), (float)199.99), evtweight, h_1d, 50, 20., 200.);
    plot1DUnderOverFlow("h_z_lep2eta" + tag_selection + flav_tag, stopt.lep2().Eta(), evtweight, h_1d, 24, -2.4, 2.4);

    float dphi_metlep = getdphi(stopt.lep1().Phi(), t1metphicorrphi);
    plot1DUnderOverFlow("h_z_dphi_metl" + tag_selection + flav_tag, dphi_metlep, evtweight, h_1d, 15, 0., 3.14159);

    if ( n_jets < 1 ) return;
    plot1DUnderOverFlow("h_z_j1pt" + tag_selection + flav_tag, min(jets.at(0).Pt(), (float)399.99), evtweight, h_1d, 20, 30., 400.);
    plot1DUnderOverFlow("h_z_j1eta" + tag_selection + flav_tag, jets.at(0).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    if ( n_jets < 2 ) return;
    plot1DUnderOverFlow("h_z_j2pt" + tag_selection + flav_tag, min(jets.at(1).Pt(), (float)299.99), evtweight, h_1d, 20, 30., 300.);
    plot1DUnderOverFlow("h_z_j2eta" + tag_selection + flav_tag, jets.at(1).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    if ( n_jets < 3 ) return;
    plot1DUnderOverFlow("h_z_j3pt" + tag_selection + flav_tag, min(jets.at(2).Pt(), (float)199.99), evtweight, h_1d, 20, 30., 200.);
    plot1DUnderOverFlow("h_z_j3eta" + tag_selection + flav_tag, jets.at(2).Eta(), evtweight, h_1d, 24, -2.4, 2.4);
    if ( n_jets < 4 ) return;
    plot1DUnderOverFlow("h_z_j4pt" + tag_selection + flav_tag, min(jets.at(3).Pt(), (float)119.99), evtweight, h_1d, 20, 30., 120.);
    plot1DUnderOverFlow("h_z_j4eta" + tag_selection + flav_tag, jets.at(3).Eta(), evtweight, h_1d, 24, -2.4, 2.4);

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
    int imaxweightcombos[2] = { -1, -1};
    double maxweightcombos[2] = { -1, -1};
    double avgweightcombos[2] = {0, 0};

    //now repeat using the Betchart solver
    //double met_x = t1metphicorr * cos(t1metphicorrphi);
    //double met_y = t1metphicorr * sin(t1metphicorrphi);

    //create lines of python code to transmit the input values
    TString l0 = Form("l0 = lv(%0.8f,%0.8f,%0.8f,%0.8f)", lepPlus.Pt(), lepPlus.Eta(), lepPlus.Phi(), lepPlus.E());
    TString l1 = Form("l1 = lv(%0.8f,%0.8f,%0.8f,%0.8f)", lepMinus.Pt(), lepMinus.Eta(), lepMinus.Phi(), lepMinus.E());
    TString j0 = Form("j0 = lv(%0.8f,%0.8f,%0.8f,%0.8f)", jet1.Pt(), jet1.Eta(), jet1.Phi(), jet1.E());
    TString j1 = Form("j1 = lv(%0.8f,%0.8f,%0.8f,%0.8f)", jet2.Pt(), jet2.Eta(), jet2.Phi(), jet2.E());
    TString metxy = Form("metx, mety = %0.8f, %0.8f", met_x, met_y);

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
    int num_err_sols[2] = {0, 0};

    for (int icombo = 0; icombo < 2; ++icombo)
    {
        if (icombo == 0) TPython::Exec("dnsC = doubleNeutrinoSolutionsCheckLinAlg((j0, j1), (l0, l1), (metx, mety))");
        if (icombo == 1) TPython::Exec("dnsC = doubleNeutrinoSolutionsCheckLinAlg((j1, j0), (l0, l1), (metx, mety))");
        TPython::Exec("dns = dnsC.dns");
        //TPython::Exec("dns = doubleNeutrinoSolutions((j0, j1), (l0, l1), (metx, mety))");
        TPython::Exec("soltest = 1");
        TPython::Exec("if dns==0: soltest = 0");
        //TPython::Exec("print soltest");
        int soltest  = TPython::Eval("soltest");


        if (soltest)
        {
            TPython::Exec("solutions = dns.nunu_s");
            //TPython::Exec("print solutions");
            TPython::Exec("nSolB = len(solutions)");

            const int nSolB  = TPython::Eval("nSolB");
            //cout<<"nSolB: "<<nSolB<<endl;
            double nusols[nSolB][2][3];

            for (int is = 0; is < nSolB; ++is)
            {
                for (int inu = 0; inu < 2; ++inu)
                {
                    for (int ix = 0; ix < 3; ++ix)
                    {
                        TString sols = Form("solutions[%0d][%0d][%0d]", is, inu, ix);
                        //cout<<sols<<endl;
                        nusols[is][inu][ix]  = TPython::Eval(sols);
                        //if(nusols[is][inu][ix] == -1) cout<<nusols[is][inu][ix]<<endl;
                    }
                }
                TLorentzVector nu1_vec , nu2_vec, lvTop1, lvTop2;
                nu1_vec.SetXYZM( nusols[is][0][0] , nusols[is][0][1] , nusols[is][0][2] , 0 );
                nu2_vec.SetXYZM( nusols[is][1][0] , nusols[is][1][1] , nusols[is][1][2] , 0 );

                //calculate t and tbar solutions and weights. The dalitz weight only distinguishes between the two combos, while the PDF weight is different for each solution.
                double sol_weight = -1;
                if (icombo == 0)
                {
                    lvTop1 = lepPlus + nu1_vec + jet1;
                    lvTop2 = lepMinus + nu2_vec + jet2;
                    sol_weight = get_pdf_weight(lvTop1, lvTop2) * get_dalitz_prob(lepPlus, lvTop1) * get_dalitz_prob(lepMinus, lvTop2);
                    //cout<<get_pdf_weight(lvTop1, lvTop2) << " " << get_dalitz_prob(lepPlus, lvTop1) << " " << get_dalitz_prob(lepMinus, lvTop2)<<endl;
                    ncombo0++;
                }
                if (icombo == 1)
                {
                    lvTop1 = lepPlus + nu1_vec + jet2;
                    lvTop2 = lepMinus + nu2_vec + jet1;
                    sol_weight = get_pdf_weight(lvTop1, lvTop2) * get_dalitz_prob(lepPlus, lvTop1) * get_dalitz_prob(lepMinus, lvTop2);
                    ncombo1++;
                }

                TLorentzVector lvW1 = lepPlus + nu1_vec;
                TLorentzVector lvW2 = lepMinus + nu2_vec;
                //cout<<"combo "<<icombo<<" solution "<<is<<" weight "<<sol_weight<<" masses: "<<lvTop1.M()<<" "<<lvTop2.M()<<" "<<lvW1.M()<<" "<<lvW2.M()<<endl;

                //don't use solutions with numerical error in solution (output masses don't match input). The input masses used are hard-coded in nuSolutions.py.
                if (  (fabs(mt_solver - lvTop1.M()) > 1.0 || fabs(mt_solver - lvTop2.M()) > 1.0) )
                {
                    num_err_sols[icombo]++;
                    continue;
                }

                if (  (fabs(mW_solver - lvW1.M()) > 1.0 || fabs(mW_solver - lvW2.M()) > 1.0) )
                {
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
        //cout<<"w: "<<AMWT_weights.at(is)<<endl;
        if (AMWT_weights.at(is) > maxweight)
        {
            imaxweight = is;
            maxweight = AMWT_weights.at(is);
        }
        if ( is < ncombo0_filled )
        {
            avgweightcombos[0] += AMWT_weights.at(is);
            if (AMWT_weights.at(is) > maxweightcombos[0])
            {
                imaxweightcombos[0] = is;
                maxweightcombos[0] = AMWT_weights.at(is);
            }
        }
        else
        {
            avgweightcombos[1] += AMWT_weights.at(is);
            if (AMWT_weights.at(is) > maxweightcombos[1])
            {
                imaxweightcombos[1] = is;
                maxweightcombos[1] = AMWT_weights.at(is);
            }
        }
    }

    //using average instead of sum with useMaxCombo gives slightly worse resolution, so commented out
    //if(ncombo0_filled>0) avgweightcombos[0] /= ncombo0_filled;
    //if(ncombo1_filled>0) avgweightcombos[1] /= ncombo1_filled;

    //cout<<AMWT_weights.size() <<" "<<ncombo0_filled <<" "<<ncombo1_filled <<" "<<imaxweight<<" "<<maxweight<<" "<<imaxweightcombos[0]<<" "<<maxweightcombos[0]<<" "<<avgweightcombos[0]<<" "<<imaxweightcombos[1]<<" "<<maxweightcombos[1]<<" "<<avgweightcombos[1]<<endl;

    if (useMaxCombo && ncombo0_filled > 0 && ncombo1_filled > 0 ) imaxweight = (avgweightcombos[0] > avgweightcombos[1]) ? imaxweightcombos[0] : imaxweightcombos[1];

    //if( (ncombo0 == 1 || ncombo1 == 1) && (ncombo0 > 1 || ncombo1 > 1) && (ncombo0 - num_err_sols[0] == 0 || ncombo1 - num_err_sols[1] == 0)  )  cout<<"all exact solutions have numerr: "<<ncombo0<<" "<<num_err_sols[0]<<" "<<ncombo1<<" "<<num_err_sols[1]<<" imax: "<<imaxweight<<" "<<imaxweightcombos[0]<<" "<<imaxweightcombos[1]<<endl;

    //don't take "closest approach" solution if exact solutions are available
    if (ncombo0 == 1)
    {
        if (ncombo1 > 1 && ncombo1_filled > 0) imaxweight = imaxweightcombos[1];
        else closestApproach = true;
    }
    if (ncombo1 == 1)
    {
        if (ncombo0 > 1 && ncombo0_filled > 0) imaxweight = imaxweightcombos[0];
        else closestApproach = true;
    }

    //cout<<imaxweight<<endl;
    //cout<<AMWT_weights.size()<<" "<<ncombo0<<" "<<ncombo1<<" "<<ncombo0_filled<<" "<<ncombo1_filled<<" "<<endl;

    if ( AMWT_weights.size() > 0)
    {
        nusum = nu1_vecs.at(imaxweight) + nu2_vecs.at(imaxweight);
        closestDeltaMET_maxwcombo = sqrt( pow( nusum.Px() - met_x , 2 ) + pow( nusum.Py() - met_y , 2 ) );

        if (closestApproach && ncombo0 == 1 && ncombo1 == 1 && ncombo0_filled == 1 && ncombo1_filled == 1 )
        {
            TLorentzVector nusum_othercombo = nu1_vecs.at(1 - imaxweight) + nu2_vecs.at(1 - imaxweight);
            closestDeltaMET_othercombo = sqrt( pow( nusum_othercombo.Px() - met_x , 2 ) + pow( nusum_othercombo.Py() - met_y , 2 ) );
            if (useClosestDeltaMET && closestDeltaMET_othercombo < closestDeltaMET_maxwcombo ) imaxweight = 1 - imaxweight;
            if (doDeltaMETcut && !useClosestDeltaMET && closestDeltaMET_maxwcombo > deltaMETcut && closestDeltaMET_othercombo < closestDeltaMET_maxwcombo ) imaxweight = 1 - imaxweight;
        }

        m_top_B = (top1_vecs.at(imaxweight).M() + top2_vecs.at(imaxweight).M()) / 2.;
        nusum = nu1_vecs.at(imaxweight) + nu2_vecs.at(imaxweight);
        closestDeltaMET_bestcombo = sqrt( pow( nusum.Px() - met_x , 2 ) + pow( nusum.Py() - met_y , 2 ) );

        //cut on deltaMET (measured vs closest solution)
        if (doDeltaMETcut && closestDeltaMET_bestcombo > deltaMETcut) m_top_B = -999.;

        if (imaxweight < 0 || imaxweight >= int(AMWT_weights.size()) ) cout << "something went wrong choosing the best solution" << endl;

    }

    if (m_top_B > 0)
    {
        top1_p4 = top1_vecs.at(imaxweight);
        top2_p4 = top2_vecs.at(imaxweight);
    }

}

double StopTreeLooper::get_pdf_weight( TLorentzVector &t1, TLorentzVector &t2 )
{
    //CM energy
    double e_com = 8000;

    //Determine x1 and x2
    double x1 = ( t1.E() + t2.E() + t1.Pz() + t2.Pz() ) / e_com;
    double x2 = ( t1.E() + t2.E() - t1.Pz() - t2.Pz() ) / e_com;

    vector <double> f1, f2;

    f1 = LHAPDF::xfxM(1, x1, mt_solver);
    f2 = LHAPDF::xfxM(1, x2, mt_solver);

    // The order of f:
    //    -t  -b  -c  -s  -u  -d   g   d   u   s   c   b   t
    //    -6  -5  -4  -3  -2  -1   0   1   2   3   4   5   6
    //     0   1   2   3   4   5   6   7   8   9   10  11  12

    double sbar1 = f1[3], sbar2 = f2[3];
    double ubar1 = f1[4], ubar2 = f2[4];
    double dbar1 = f1[5], dbar2 = f2[5];
    double g1    = f1[6], g2    = f2[6];
    double d1    = f1[7], d2    = f2[7];
    double u1    = f1[8], u2    = f2[8];
    double s1    = f1[9], s2    = f2[9];

    //Should glue-glue be doubled? Probably not, but plot histo later
    double pdf_prob = (u1 * ubar2 + u2 * ubar1 +
                       d1 * dbar2 + d2 * dbar1 +
                       s1 * sbar2 + s2 * sbar1 +
                       g1 * g2);

    return pdf_prob;
}

double StopTreeLooper::get_dalitz_prob( TLorentzVector &lep, TLorentzVector &top )
{
    double mb = mb_solver;
    double mw = mW_solver;
    double mte = lep.Dot( top );
    double mt = top.M();
    double mt2 = mt * mt;
    double mb2 = mb * mb;
    double mw2 = mw * mw;
    double mt2_mb2 = mt2 - mb2;

    return 4. * mte * ( mt2 - mb2 - 2. * mte ) / ( mt2_mb2 * mt2_mb2 + mw2 * ( mt2 + mb2 ) - 2. * mw2 * mw2 );
}


bool StopTreeLooper::passFullSelection(bool isData)
{
    //isolated track veto. Not used because it doesn't help signal/bkg because the main background has exactly 2 leptons, and because it's difficult for MC to model it well.
    //if ( stopt.trkpt10loose() > 0. && stopt.trkreliso10loose() < 0.1 ) return false;

    //tau veto
    //if (!passTauVeto()) return false;

    //if ( pfcalo_deltamet > 80. ) continue; //can't really cut on this because the data/MC agreement is very poor

    bool passFull = false;
    if ( passDileptonSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
            && n_bjets > 0
            && n_jets > 1
            && stopt.dilmass() >= 20.0
            && (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}



bool StopTreeLooper::passFullSelection_bveto(bool isData)
{
    bool passFull = false;
    if ( passDileptonSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
            && n_bjets == 0
            && n_jets > 1
            && stopt.dilmass() >= 20.0
            && (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}


bool StopTreeLooper::passFullSelection_1jet_inclusiveb(bool isData)
{
    bool passFull = false;
    if ( passDileptonSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
            //&& n_bjets == 1
            && n_jets == 1
            && stopt.dilmass() >= 20.0
            && (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}

bool StopTreeLooper::passFullSelection_0jets(bool isData)
{
    bool passFull = false;
    if ( passDileptonSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
            //&& n_bjets == 0
            && n_jets == 0
            && stopt.dilmass() >= 20.0
            && (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}


bool StopTreeLooper::passFullSelection_Zpeak_inclusiveb(bool isData)
{
    bool passFull = false;
    if ( passDileptonSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && fabs( stopt.dilmass() - 91.) <= 15.
            //&& n_bjets > 0
            && n_jets > 1
            && stopt.dilmass() >= 20.0
            && (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}


bool StopTreeLooper::passFullSelection_noMETcut_inclusiveb(bool isData)
{
    bool passFull = false;
    if ( passDileptonSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
            //&& n_bjets > 0
            && n_jets > 1
            && stopt.dilmass() >= 20.0
            //&& (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}



bool StopTreeLooper::passFullSelection_SS_inclusiveb(bool isData)
{
    bool passFull = false;
    if ( passDileptonSSSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
            //&& n_bjets > 0
            && n_jets > 1
            && stopt.dilmass() >= 20.0
            && (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}


bool StopTreeLooper::passFullSelection_SS_noMETcut_inclusiveb(bool isData)
{
    bool passFull = false;
    if ( passDileptonSSSelectionWithEndcapEls(isData)
            //events with mlb_min beyond the kinematic edge cannot be solved and have lower signal:background ratio. Could also cut on second-lowest mlb, but this would remove some events where we got one of the bs wrong that are still OK for the purely leptonic variables.
            && mlb_min <= mlb_max
            && (abs(stopt.id1()) != abs(stopt.id2()) || fabs( stopt.dilmass() - 91.) > 15. )
            //&& n_bjets > 0
            && n_jets > 1
            && stopt.dilmass() >= 20.0
            //&& (abs(stopt.id1()) != abs(stopt.id2()) || t1metphicorr >= 40. )
       ) passFull = true;

    return passFull;
}








//-------------------------------------
// Book the baby ntuple
//-------------------------------------
void StopTreeLooper::MakeBabyNtuple(const char *babyFilename)
{
    TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
    rootdir->cd();
    babyFile_ = new TFile(Form("%s", babyFilename), "RECREATE");
    babyFile_->cd();
    babyTree_ = new TTree("tree", "A Baby Ntuple");

    babyTree_->Branch("run",                   &run,                 "run/I"                  );
    babyTree_->Branch("ls",                    &ls,                  "ls/I"                   );
    babyTree_->Branch("evt",                   &evt,                 "evt/I"                  );
    babyTree_->Branch("channel",                &channel,              "channel/I"               );
    babyTree_->Branch("genchannel",                &genchannel,              "genchannel/I"               );
    babyTree_->Branch("t_mass",                &m_top,              "t_mass/F"               );
    babyTree_->Branch("weight",                &weight,              "weight/D"               );
    babyTree_->Branch("PDFsystweights", "std::vector<double>", &PDFsystweights );
    //babyTree_->Branch("Nsolns",                &Nsolns_,              "Nsolns/I"               );
    //babyTree_->Branch("massltb",               &massltb_,             "massltb/F"              );
    //babyTree_->Branch("massllb",               &massllb_,             "massllb/F"              );
    //babyTree_->Branch("dr_ltjet_gen",          &dr_ltjet_gen_,        "dr_ltjet_gen/F"         );
    //babyTree_->Branch("dr_lljet_gen",          &dr_lljet_gen_,        "dr_lljet_gen/F"         );
    //babyTree_->Branch("ndavtx",                &ndavtx_,              "ndavtx/I"               );
    babyTree_->Branch("dilmass",               &dilmass,             "dilmass/F"              );
    babyTree_->Branch("tt_mass",               &tt_mass,             "tt_mass/F"              );
    //babyTree_->Branch("ttRapidity",            &ttRapidity_,          "ttRapidity/F"           );
    babyTree_->Branch("ttRapidity2",            &ttRapidity2,          "ttRapidity2/F"           );
    babyTree_->Branch("ttRapidity2_gen",            &ttRapidity2_gen,          "ttRapidity2_gen/F"           );
    babyTree_->Branch("ttPt",                   &tt_pT,          "ttPt/F"           );
    babyTree_->Branch("ttPt_gen",                   &tt_pT_gen,          "ttPt_gen/F"           );
    babyTree_->Branch("lep_charge_asymmetry",  &lep_charge_asymmetry, "lep_charge_asymmetry/F" );
    //babyTree_->Branch("lep_pseudorap_diff",    &lep_pseudorap_diff,  "lep_pseudorap_diff/F"   );
    babyTree_->Branch("lep_azimuthal_asymmetry", &lep_azimuthal_asymmetry,  "lep_azimuthal_asymmetry/F"   );
    babyTree_->Branch("lep_azimuthal_asymmetry2", &lep_azimuthal_asymmetry2,  "lep_azimuthal_asymmetry2/F"   );
    babyTree_->Branch("top_spin_correlation",  &top_spin_correlation, "top_spin_correlation/F" );
    babyTree_->Branch("lep_cos_opening_angle",  &lep_cos_opening_angle, "lep_cos_opening_angle/F" );
    babyTree_->Branch("top_costheta_cms",      &top_costheta_cms,    "top_costheta_cms/F"     );
    babyTree_->Branch("lepPlus_costheta_cms",      &lepPlus_costheta_cms, "lepPlus_costheta_cms/F"     );
    babyTree_->Branch("lepMinus_costheta_cms",      &lepMinus_costheta_cms, "lepMinus_costheta_cms/F"     );
    babyTree_->Branch("top_rapiditydiff_cms",      &top_rapiditydiff_cms, "top_rapiditydiff_cms/F"     );
    babyTree_->Branch("top_rapiditydiff_Marco",      &top_rapiditydiff_Marco, "top_rapiditydiff_Marco/F"     );
    babyTree_->Branch("top_pseudorapiditydiff_cms",      &top_pseudorapiditydiff_cms, "top_pseudorapiditydiff_cms/F"     );
    babyTree_->Branch("top_rapiditydiff_cms_gen",      &top_rapiditydiff_cms_gen, "top_rapiditydiff_cms_gen/F"     );
    babyTree_->Branch("top_rapiditydiff_Marco_gen",      &top_rapiditydiff_Marco_gen, "top_rapiditydiff_Marco_gen/F"     );
    babyTree_->Branch("top_pseudorapiditydiff_cms_gen",      &top_pseudorapiditydiff_cms_gen, "top_pseudorapiditydiff_cms_gen/F"     );
    babyTree_->Branch("tt_mass_gen",           &tt_mass_gen,          "tt_mass_gen/F"              );
    //babyTree_->Branch("ttRapidity_gen",            &ttRapidity_gen,          "ttRapidity_gen/F"            );
    babyTree_->Branch("lep_charge_asymmetry_gen",  &lep_charge_asymmetry_gen, "lep_charge_asymmetry_gen/F" );
    babyTree_->Branch("lep_azimuthal_asymmetry_gen",    &lep_azimuthal_asymmetry_gen,  "lep_azimuthal_asymmetry_gen/F"   );
    babyTree_->Branch("lep_azimuthal_asymmetry2_gen",    &lep_azimuthal_asymmetry2_gen,  "lep_azimuthal_asymmetry2_gen/F"   );
    babyTree_->Branch("top_spin_correlation_gen",  &top_spin_correlation_gen, "top_spin_correlation_gen/F"  );
    babyTree_->Branch("lep_cos_opening_angle_gen",  &lep_cos_opening_angle_gen, "lep_cos_opening_angle_gen/F"  );
    babyTree_->Branch("top_costheta_cms_gen",      &top_costheta_cms_gen,    "top_costheta_cms_gen/F"      );
    babyTree_->Branch("lepPlus_costheta_cms_gen",      &lepPlus_costheta_cms_gen, "lepPlus_costheta_cms_gen/F"      );
    babyTree_->Branch("lepMinus_costheta_cms_gen",      &lepMinus_costheta_cms_gen, "lepMinus_costheta_cms_gen/F"      );
}

//----------------------------------
// Fill the baby
//----------------------------------
void StopTreeLooper::FillBabyNtuple()
{
    babyTree_->Fill();
}

//--------------------------------
// Close the baby
//--------------------------------
void StopTreeLooper::CloseBabyNtuple()
{
    babyFile_->cd();
    babyTree_->Write();
    babyFile_->Close();
}


