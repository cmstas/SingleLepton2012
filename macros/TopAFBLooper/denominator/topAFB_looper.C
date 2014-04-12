#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <set>
#include <algorithm>
#include <map>


#include "TPython.h"
// ROOT includes
#include "TSystem.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TChainElement.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "TLorentzVector.h"
// TAS includes
#include "../../../CORE/CMS2.h"
#include "../../../CORE/mcSelections.h"
//#include "../CORE/utilities.cc"
#include "./topAFB_looper.h"
#include "Histograms.cc"

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;


double topAFB_looper::TopPtWeight(double topPt)
{
    if ( topPt < 0 ) return 1;
    if (topPt > 400) topPt = 400;

    //double result = (1.4 / 1000000.0) * topPt * topPt - (2.0 / 1000.0) * topPt + 1.2; //old 7TeV fit (l+j and 2l combined)
    double result = exp(0.199 - 0.00166 * topPt); //new 7TeV fit (l+j and 2l combined)
    //note this fit is for data/madgraph, and we are using MC@NLO



    return result;
}


topAFB_looper::topAFB_looper()
{
    applyTopPtWeighting = false;
    weighttaudecay = false;
}


topAFB_looper::~topAFB_looper()
{

}


void topAFB_looper::ScanChain(TChain *chain, vector<TString> v_Cuts, string prefix, float lumi )
{
    //deal with the cuts
    applyTopPtWeighting = find(v_Cuts.begin(), v_Cuts.end(), "applyTopPtWeighting" ) != v_Cuts.end();
    weighttaudecay = find(v_Cuts.begin(), v_Cuts.end(), "weighttaudecay") != v_Cuts.end();

    bookHistos(prefix.c_str(), 1, 1);
    bool isData = false;
    bool applyNoCuts = true;

    //--------------------------
    // File and Event Loop
    //---------------------------
    TObjArray *listOfFiles = chain->GetListOfFiles();
    unsigned int nEventsChain = 0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    double nEvents_noCuts = 0;
    double nEvents_noCuts_dil = 0;
    ULong64_t nEvents_noCuts_novtxweight = 0;
    ULong64_t nEvents_noCuts_novtxweight_dil = 0;
    double nSelectedEvents = 0;
    unsigned int npreSelectedEvents = 0;
    unsigned int npreSelectedEvents_genmatch1 = 0;
    unsigned int npreSelectedEvents_genmatch2 = 0;
    TIter fileIter(listOfFiles);
    map<int, int> m_events;
    while (TChainElement *currentFile = (TChainElement *)fileIter.Next() )
    {

        TString filename = currentFile->GetTitle();

        TFile f(filename.Data());
        TTree *tree = (TTree *)f.Get("Events");
        TTreeCache::SetLearnEntries(10);
        tree->SetCacheSize(128 * 1024 * 1024);
        cms2.Init(tree);
        unsigned int nEvents = tree->GetEntries();
        bool print_evt_weight = true;
        for (unsigned int event = 0; event < nEvents; ++event)
        {
            // Event Loop
            tree->LoadTree(event);
            cms2.GetEntry(event);


            if (print_evt_weight && !isData)
            {
                cout << "Event Weight = " << evt_scale1fb() * lumi << endl;
                print_evt_weight = false;
            }

            int myType = 2;
            unsigned int jetBin = 2;
            int ndavtx = 0;
            double weight = 1.0;
            float lepPlus_costheta_cms , lep_azimuthal_asymmetry , lep_azimuthal_asymmetry_2 , lep_charge_asymmetry , lep_pseudorap_diff , top_costheta_cms;
            float lepMinus_costheta_cms;
            float top_pseudorapiditydiff_cms , top_rapiditydiff_Marco , top_rapiditydiff_cms , top_spin_correlation , ttRapidity , ttRapidity2, tt_mass , tt_mass_nojetsmear , tt_pT;
            float lep_cos_opening_angle;
            double top1sdp = -999.;
            double top2sdp = -999.;
            float m_top = -999.0;
            float m_top_S = -999.0;
            float m_top_B = -999.0;
            float m_top_nojetsmear = -999.0;
            double mass_ltb, mass_llb;

            vector <float> AMWTweight_nojetsmear;
            vector <TLorentzVector> top1_p4, top2_p4, top1_nojetsmear_p4, top2_nojetsmear_p4;
            vector <TLorentzVector> nu1_vecs , nu2_vecs;
            vector <TLorentzVector> top1_vecs , top2_vecs;
            vector <double> AMWT_weights;
            int imaxweight = -1;
            double maxweight = -1;
            int imaxweightcombos[2] = { -1, -1};
            double maxweightcombos[2] = { -1, -1};
            double avgweightcombos[2] = {0, 0};
            bool closestApproach = false;
            double closestDeltaMET_maxwcombo = -999;
            double closestDeltaMET_othercombo = -999;
            double closestDeltaMET_bestcombo = -999;
            TLorentzVector nusum;
            TLorentzVector cms, lepPlus, lepMinus, jet1, jet2;
            int Nsolns = 1;

            float ndavtxweight = 1.;

            //get the channels correct
            int nels = 0;
            int nmus = 0;
            int ntaus = 0;
            int nleps = 0;
            if (!isData)
                nleps = leptonGenpCount_lepTauDecays(nels, nmus, ntaus);

            if (prefix == "ttdil"    &&  nleps == 2) ++nEventsTotal;
            // Progress feedback to the user
            if (nEventsTotal % 2000 == 0)
            {
                // xterm magic from L. Vacavant and A. Cerri
                if (isatty(1))
                {
                    printf("\015\033[32m ---> \033[1m\033[31m%4.1f%%"
                           "\033[0m\033[32m <---\033[0m\015", (float)nEventsTotal / (nEventsChain * 0.01));
                    fflush(stdout);
                }
            }//if(nEventsTotal%20000 == 0) {


            //daughter lepton angle in tau rest frame to check if MC is correctly using the tau polarisation
            double weight_taudecay = 1.;
            if (ntaus > 0)
            {

                TLorentzVector lepPlus_gen(0, 0, 0, 0), lepMinus_gen(0, 0, 0, 0), lepPlus_gen_status1(0, 0, 0, 0), lepMinus_gen_status1(0, 0, 0, 0);
                bool lepPlusIsTau = false;
                bool lepPlusIsTauEle = false;
                bool lepPlusIsTauMuo = false;
                bool lepMinusIsTau = false;
                bool lepMinusIsTauEle = false;
                bool lepMinusIsTauMuo = false;

                for (unsigned int i = 0; i < genps_p4().size(); i++)
                {
                    if (genps_status()[i] == 3 && fabs(genps_id_mother()[i]) == 24 )
                    {

                        if (genps_id()[i] == -15 )
                        {
                            //if(genps_id()[i]==-15 && genps_lepdaughter_id()[i].size()==3 ) {   //cosTheta* distribution is only flat for unpolarised tau decay in the simple case of 3 daughters
                            double lepPlus_genX = 0;
                            double lepPlus_genY = 0;
                            double lepPlus_genZ = 0;
                            double lepPlus_genT = 0;

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 12)
                                {
                                    lepPlusIsTauEle = true;
                                    lepPlusIsTau = true;
                                    break;
                                }
                                if ( daughterID == 14)
                                {
                                    lepPlusIsTauMuo = true;
                                    lepPlusIsTau = true;
                                    break;
                                }
                                // only neutrinos are definitely from leptonic tau decay (see comment for leptonGenpCount_lepTauDecays)
                            }

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 11 && lepPlusIsTauEle)
                                {
                                    lepPlus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepPlusIsTauEle = false;  //so that any secondary daughter electrons don't overwrite the first one
                                }
                                if ( daughterID == 13 && lepPlusIsTauMuo)
                                {
                                    lepPlus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepPlusIsTauMuo = false;  //so that any secondary daughter muons don't overwrite the first one
                                }
                                lepPlus_genX += genps_lepdaughter_p4()[i][kk].x();
                                lepPlus_genY += genps_lepdaughter_p4()[i][kk].y();
                                lepPlus_genZ += genps_lepdaughter_p4()[i][kk].z();
                                lepPlus_genT += genps_lepdaughter_p4()[i][kk].T();
                                lepPlus_gen.SetXYZT(lepPlus_genX, lepPlus_genY, lepPlus_genZ, lepPlus_genT);
                            }

                        }

                        if (genps_id()[i] == 15 )
                        {
                            //if(genps_id()[i]==15 && genps_lepdaughter_id()[i].size()==3 ) {   //cosTheta* distribution is only flat for unpolarised tau decay in the simple case of 3 daughters
                            double lepMinus_genX = 0;
                            double lepMinus_genY = 0;
                            double lepMinus_genZ = 0;
                            double lepMinus_genT = 0;

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 12)
                                {
                                    lepMinusIsTauEle = true;
                                    lepMinusIsTau = true;
                                    break;
                                }
                                if ( daughterID == 14)
                                {
                                    lepMinusIsTauMuo = true;
                                    lepMinusIsTau = true;
                                    break;
                                }
                                // only neutrinos are definitely from leptonic tau decay (see comment for leptonGenpCount_lepTauDecays)
                            }

                            for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                            {
                                int daughterID = abs(genps_lepdaughter_id()[i][kk]);
                                if ( daughterID == 11 && lepMinusIsTauEle)
                                {
                                    lepMinus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepMinusIsTauEle = false;  //so that any secondary daughter electrons don't overwrite the first one
                                }
                                if ( daughterID == 13 && lepMinusIsTauMuo)
                                {
                                    lepMinus_gen_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                    lepMinusIsTauMuo = false;  //so that any secondary daughter muons don't overwrite the first one
                                }
                                lepMinus_genX += genps_lepdaughter_p4()[i][kk].x();
                                lepMinus_genY += genps_lepdaughter_p4()[i][kk].y();
                                lepMinus_genZ += genps_lepdaughter_p4()[i][kk].z();
                                lepMinus_genT += genps_lepdaughter_p4()[i][kk].T();
                                lepMinus_gen.SetXYZT(lepMinus_genX, lepMinus_genY, lepMinus_genZ, lepMinus_genT);
                            }

                        }

                    }
                }

                bool fillEvent = true;
                if (prefix == "ttdil"    &&  nleps != 2) fillEvent = false;
                if (prefix == "ttotr"    &&  nleps == 2) fillEvent = false;
                if (prefix == "DYee"     &&  nels != 2) fillEvent = false;
                if (prefix == "DYmm"     &&  nmus != 2) fillEvent = false;
                if (prefix == "DYtautau" &&  ntaus != 2) fillEvent = false;

                if (lepPlusIsTau)
                {
                    lepPlus_gen_status1.Boost(-lepPlus_gen.BoostVector());
                    double cosTheta_lepPlus_status1 = lepPlus_gen_status1.Vect().Dot(lepPlus_gen.Vect()) / (lepPlus_gen_status1.Vect().Mag() * lepPlus_gen.Vect().Mag());
                    double EoverEmax_lepPlus = 2.*lepPlus_gen_status1.E() / lepPlus_gen.M();
                    if (EoverEmax_lepPlus > 1.) EoverEmax_lepPlus = 1.;
                    //if(EoverEmax_lepPlus>1.) {cout<<"**********x>1********"<<endl; cout<<EoverEmax_lepPlus<<" "<<lepPlus_gen_status1.E()<<" "<<lepPlus_gen.M()<<endl;}
                    //double weight_lepPlus = 1.+cosTheta_lepPlus_status1/3.;  //approximation: cosTheta dependence varies with x
                    double weight_lepPlus = 1. + cosTheta_lepPlus_status1 * ( 1. - 2.*EoverEmax_lepPlus ) / ( 2.*EoverEmax_lepPlus - 3. );
                    //cout<<"P "<<(weighttaudecay?weight_lepPlus:1.)<<endl;
                    if (fillEvent) fillHistos( hlepPlusCosThetaTau_gen,  cosTheta_lepPlus_status1, weighttaudecay ? weight_lepPlus : 1., myType, jetBin);
                    if (fillEvent) fillHistos( hlepPlusxTau_gen,  EoverEmax_lepPlus, weighttaudecay ? weight_lepPlus : 1., myType, jetBin);
                    weight_taudecay *= weight_lepPlus;
                }
                if (lepMinusIsTau)
                {
                    lepMinus_gen_status1.Boost(-lepMinus_gen.BoostVector());
                    double cosTheta_lepMinus_status1 = lepMinus_gen_status1.Vect().Dot(lepMinus_gen.Vect()) / (lepMinus_gen_status1.Vect().Mag() * lepMinus_gen.Vect().Mag());
                    double EoverEmax_lepMinus = 2.*lepMinus_gen_status1.E() / lepMinus_gen.M();
                    if (EoverEmax_lepMinus > 1.) EoverEmax_lepMinus = 1.;
                    //if(EoverEmax_lepMinus>1.) {cout<<"**********x>1********"<<endl; cout<<EoverEmax_lepMinus<<" "<<lepMinus_gen_status1.E()<<" "<<lepMinus_gen.M()<<endl;}
                    //double weight_lepMinus = 1.+cosTheta_lepMinus_status1/3.;  //approximation: cosTheta dependence varies with x
                    double weight_lepMinus = 1. + cosTheta_lepMinus_status1 * ( 1. - 2.*EoverEmax_lepMinus ) / ( 2.*EoverEmax_lepMinus - 3. );
                    //cout<<"M "<<(weighttaudecay?weight_lepMinus:1.)<<endl;
                    if (fillEvent) fillHistos( hlepMinusCosThetaTau_gen,  cosTheta_lepMinus_status1, weighttaudecay ? weight_lepMinus : 1., myType, jetBin);
                    if (fillEvent) fillHistos( hlepMinusxTau_gen,  EoverEmax_lepMinus, weighttaudecay ? weight_lepMinus : 1., myType, jetBin);
                    weight_taudecay *= weight_lepMinus;
                }

            } //ntaus>0
            //cout<<weight_taudecay<<endl;

            /*
            // Caculate PDF Systematics
            double pdf_weight = 1.0;
            if (applyPDFWeight && prefix == "ttdil"  ){

                for (unsigned int subset = 0; subset < nsets; subset++)
                {

                  // std::cout << "doing set, subset: " << set_ << ", " << subset << std::endl;
                  LHAPDF::initPDFM(2,0);

                  float   x1          = 0.0;      // momentum fraction for parton1
                  float   x2          = 0.0;
                  int     id1         = 0;        // pdgid of parton1
                  int     id2         = 0;
                  float   Q           = 0.0;      // event momentum scale
                  x1          = cms2.pdfinfo_x1();
                  x2          = cms2.pdfinfo_x2();
                  id1         = cms2.pdfinfo_id1();
                  id2         = cms2.pdfinfo_id2();
                  Q           = cms2.pdfinfo_scale();
                  // generated pdf values
                  double fx1Q0gen = LHAPDF::xfxM(genset_, x1, Q, id1) / x1;
                  double fx2Q0gen = LHAPDF::xfxM(genset_, x2, Q, id2) / x2;
                  // subset pdf values
                  double fx1Qi = LHAPDF::xfxM(set_, x1, Q, id1) / x1;
                  double fx2Qi = LHAPDF::xfxM(set_, x2, Q, id2) / x2;
                  // calculate weight and fill histogram
                  pdf_weight = ((fx1Qi*fx2Qi)/(fx1Q0gen*fx2Q0gen));
                  //cout << fx1Qi <<endl;
                  // cout << fx1Q0gen <<endl;
                  // cout << pdf_weight <<endl;


                }// end of loop over subset of PDFs
            }
            */



            if (!isData) weight = evt_scale1fb() * lumi;
            //negative weights for MC@NLO
            if (prefix == "ttdil" || prefix == "ttotr") weight = weight * (fabs(genps_weight()) / genps_weight());
            //tau decay cosTheta* weighting
            if (weighttaudecay && (prefix == "ttdil" || prefix == "ttotr")  && ntaus > 0) weight *= weight_taudecay;

            nEvents_noCuts_novtxweight += 1;
            if (isData) nEvents_noCuts += 1.;
            else  nEvents_noCuts += ndavtxweight;

            if (prefix == "ttdil"    &&  nleps != 2) continue;
            if (prefix == "ttotr"    &&  nleps == 2) continue;
            if (prefix == "DYee"     &&  nels != 2) continue;
            if (prefix == "DYmm"     &&  nmus != 2) continue;
            if (prefix == "DYtautau" &&  ntaus != 2) continue;

            nEvents_noCuts_novtxweight_dil += 1;
            if (isData) nEvents_noCuts_dil += 1.;
            else  nEvents_noCuts_dil += ndavtxweight;


            if (true)
            {

                unsigned int hypIdx = 999;
                vector<LorentzVector>  v_goodJets_cand_p4;
                vector<unsigned int> v_goodJets;
                unsigned int nJets = 0;
                int nBtagJets = -1;
                LorentzVector lt_p4;
                LorentzVector ll_p4;
                pair<float, float> p_met; //met and met phi
                float thefirstJet_pt = -999;
                float thesecondJet_pt = -999;
                float theSumBtagJetPt = -999;
                int i_ltbjet = -1;
                int i_llbjet = -1;


                if ( applyTopPtWeighting && (prefix == "ttdil" || prefix == "ttotr") )
                {

                    TLorentzVector topplus_genp_p4(0, 0, 0, 0), topminus_genp_p4(0, 0, 0, 0);
                    bool hastop = false;
                    bool hastbar = false;

                    for (unsigned int i = 0; i < genps_p4().size(); i++)
                    {
                        if (genps_status()[i] == 3)
                        {

                            if (genps_id()[i] == 6 )
                            {
                                topplus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                         genps_p4()[i].y(),
                                                         genps_p4()[i].z(),
                                                         genps_p4()[i].t()
                                                       );
                                hastop = true;
                            }
                            else if (genps_id()[i] == -6 )
                            {
                                topminus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                          genps_p4()[i].y(),
                                                          genps_p4()[i].z(),
                                                          genps_p4()[i].t()
                                                        );
                                hastbar = true;
                            }

                        }
                    }

                    float pT_topplus_gen = topplus_genp_p4.Pt();
                    float pT_topminus_gen = topminus_genp_p4.Pt();
                    if (hastop && hastbar)
                    {

                        weight = weight * sqrt( TopPtWeight(pT_topplus_gen) * TopPtWeight(pT_topminus_gen) );
                    }
                }


                int i_smear = 0;


                float dr_ltjet_gen = -999.0;
                float dr_lljet_gen = -999.0;
                float tt_mass_gen;
                float tt_pT_gen;
                float ttRapidity_gen;
                float ttRapidity2_gen;
                float top_costheta_cms_gen;
                float lep_charge_asymmetry_gen;
                float lep_azimuthal_asymmetry_gen;
                float lep_azimuthal_asymmetry2_gen;
                float top_spin_correlation_gen;
                float lep_cos_opening_angle_gen;
                float lepPlus_costheta_cms_gen;
                float lepMinus_costheta_cms_gen;
                float top_rapiditydiff_cms_gen;
                float top_rapiditydiff_Marco_gen;
                float top_pseudorapiditydiff_cms_gen;
                // generator level plots
                //if(!isData && (prefix == "ttdil"|| prefix == "wprime400"|| prefix == "wprime600" || prefix == "axigluonR")){
                if (!isData)
                {



                    TLorentzVector topplus_genp_p4(0, 0, 0, 0), topminus_genp_p4(0, 0, 0, 0), cms_gen(0, 0, 0, 0), lepPlus_gen(0, 0, 0, 0), lepMinus_gen(0, 0, 0, 0), bPlus_gen(0, 0, 0, 0), bMinus_gen(0, 0, 0, 0), nuPlus_gen(0, 0, 0, 0), nuMinus_gen(0, 0, 0, 0);
                    TLorentzVector topPlus_status1(0, 0, 0, 0), topMinus_status1(0, 0, 0, 0), WPlus_status1(0, 0, 0, 0), WMinus_status1(0, 0, 0, 0), lepPlus_status1(0, 0, 0, 0), lepMinus_status1(0, 0, 0, 0), bPlus_status1(0, 0, 0, 0), bMinus_status1(0, 0, 0, 0), nuPlus_status1(0, 0, 0, 0), nuMinus_status1(0, 0, 0, 0);

                    bool from_gluon = true;
                    for (unsigned int i = 0; i < genps_p4().size(); i++)
                    {
                        if (genps_status()[i] == 3)
                        {
                            if ((genps_id_mother()[i] == 6 || genps_id_mother()[i] == 24 ))
                            {
                                if ( (genps_id()[i] == -11 || genps_id()[i] == -13 ||  genps_id()[i] == -15) )
                                {
                                    lepPlus_gen.SetXYZT(genps_p4()[i].x(),
                                                        genps_p4()[i].y(),
                                                        genps_p4()[i].z(),
                                                        genps_p4()[i].t()
                                                       );

                                    //status = 1 lepton
                                    if ( genps_id()[i] != -15 )
                                    {
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[i][kk];
                                            if ( daughterID == genps_id()[i] )
                                            {
                                                lepPlus_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                                continue;
                                            }
                                            //need to add all status=1 photons in a DR<0.1 cone around the lepton.
                                        }
                                    }

                                }
                                else if ( (genps_id()[i] == 12 || genps_id()[i] == 14 ||  genps_id()[i] == 16) )
                                {
                                    nuPlus_gen.SetXYZT(genps_p4()[i].x(),
                                                       genps_p4()[i].y(),
                                                       genps_p4()[i].z(),
                                                       genps_p4()[i].t()
                                                      );

                                    //status = 1 neutrino
                                    if ( genps_id()[i] != 16 )
                                    {
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[i][kk];
                                            if ( daughterID == genps_id()[i] )
                                            {
                                                nuPlus_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                                continue;
                                            }
                                        }
                                    }

                                }
                                else if ( genps_id()[i] == 5)
                                {
                                    bPlus_gen.SetXYZT(genps_p4()[i].x(),
                                                      genps_p4()[i].y(),
                                                      genps_p4()[i].z(),
                                                      genps_p4()[i].t()
                                                     );
                                }
                            }
                            else if ( (genps_id_mother()[i] == -6 || genps_id_mother()[i] == -24 ))
                            {
                                if ( (genps_id()[i] == 11 || genps_id()[i] == 13 ||  genps_id()[i] == 15) )
                                {

                                    lepMinus_gen.SetXYZT( genps_p4()[i].x(),
                                                          genps_p4()[i].y(),
                                                          genps_p4()[i].z(),
                                                          genps_p4()[i].t()
                                                        );

                                    //status = 1 lepton
                                    if ( genps_id()[i] != 15 )
                                    {
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[i][kk];
                                            if ( daughterID == genps_id()[i] )
                                            {
                                                lepMinus_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                                continue;
                                            }
                                            //need to add all status=1 photons in a DR<0.1 cone around the lepton.
                                        }
                                    }

                                }
                                else if ( (genps_id()[i] == -12 || genps_id()[i] == -14 ||  genps_id()[i] == -16) )
                                {

                                    nuMinus_gen.SetXYZT( genps_p4()[i].x(),
                                                         genps_p4()[i].y(),
                                                         genps_p4()[i].z(),
                                                         genps_p4()[i].t()
                                                       );

                                    //status = 1 neutrino
                                    if ( genps_id()[i] != -16 )
                                    {
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[i].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[i][kk];
                                            if ( daughterID == genps_id()[i] )
                                            {
                                                nuMinus_status1.SetXYZT( genps_lepdaughter_p4()[i][kk].x(), genps_lepdaughter_p4()[i][kk].y(), genps_lepdaughter_p4()[i][kk].z(), genps_lepdaughter_p4()[i][kk].t() );
                                                continue;
                                            }
                                        }
                                    }

                                }
                                else if ( genps_id()[i] == -5)
                                {
                                    bMinus_gen.SetXYZT(genps_p4()[i].x(),
                                                       genps_p4()[i].y(),
                                                       genps_p4()[i].z(),
                                                       genps_p4()[i].t()
                                                      );
                                }
                            }

                            if (genps_id()[i] == 6 )
                            {
                                topplus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                         genps_p4()[i].y(),
                                                         genps_p4()[i].z(),
                                                         genps_p4()[i].t()
                                                       );
                                if (abs(genps_id_mother()[i]) == 21)
                                {
                                    from_gluon = true;
                                }
                                else
                                {
                                    from_gluon = false;
                                }
                            }
                            else if (genps_id()[i] == -6 )
                            {
                                topminus_genp_p4.SetXYZT( genps_p4()[i].x(),
                                                          genps_p4()[i].y(),
                                                          genps_p4()[i].z(),
                                                          genps_p4()[i].t()
                                                        );
                                if (abs(genps_id_mother()[i]) == 21)
                                {
                                    from_gluon = true;
                                }
                                else
                                {
                                    from_gluon = false;
                                }

                            }

                        } //status = 3

                    } //genparticles loop




                    float m_topminus_gen = topminus_genp_p4.M();
                    float m_topplus_gen = topplus_genp_p4.M();

                    tt_mass_gen = (topplus_genp_p4 + topminus_genp_p4).M();
                    ttRapidity_gen = topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity();
                    ttRapidity2_gen = (topplus_genp_p4 + topminus_genp_p4).Rapidity();
                    //ttRapidity_gen = topplus_genp_p4.Eta() + topminus_genp_p4.Eta();

                    top_rapiditydiff_cms_gen = (topplus_genp_p4.Rapidity() - topminus_genp_p4.Rapidity()) * (topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity());
                    top_pseudorapiditydiff_cms_gen = abs(topplus_genp_p4.Eta()) - abs(topminus_genp_p4.Eta());
                    top_rapiditydiff_Marco_gen = abs(topplus_genp_p4.Rapidity()) - abs(topminus_genp_p4.Rapidity());


                    cms_gen = topplus_genp_p4 + topminus_genp_p4;
                    tt_pT_gen = cms_gen.Pt();
                    topplus_genp_p4.Boost(-cms_gen.BoostVector());
                    topminus_genp_p4.Boost(-cms_gen.BoostVector());
                    top_costheta_cms_gen = topplus_genp_p4.Vect().Dot(cms_gen.Vect()) / (topplus_genp_p4.Vect().Mag() * cms_gen.Vect().Mag());


                    lep_charge_asymmetry_gen = abs(lepPlus_gen.Eta()) - abs(lepMinus_gen.Eta());
                    lep_azimuthal_asymmetry_gen = lepPlus_gen.DeltaPhi(lepMinus_gen);
                    lep_azimuthal_asymmetry2_gen = acos(cos(lep_azimuthal_asymmetry_gen));

                    lepPlus_gen.Boost(-cms_gen.BoostVector());
                    lepPlus_gen.Boost(-topplus_genp_p4.BoostVector());
                    lepMinus_gen.Boost(-cms_gen.BoostVector());
                    lepMinus_gen.Boost(-topminus_genp_p4.BoostVector());

                    lepPlus_costheta_cms_gen = lepPlus_gen.Vect().Dot(topplus_genp_p4.Vect()) / (lepPlus_gen.Vect().Mag() * topplus_genp_p4.Vect().Mag());
                    lepMinus_costheta_cms_gen = lepMinus_gen.Vect().Dot(topminus_genp_p4.Vect()) / (lepMinus_gen.Vect().Mag() * topminus_genp_p4.Vect().Mag());

                    top_spin_correlation_gen = lepPlus_costheta_cms_gen * lepMinus_costheta_cms_gen;
                    lep_cos_opening_angle_gen = lepPlus_gen.Vect().Dot(lepMinus_gen.Vect()) / (lepPlus_gen.Vect().Mag() * lepMinus_gen.Vect().Mag());

                    fillHistos( htopMass_plus_gen, m_topplus_gen ,  weight, myType, jetBin, Nsolns);
                    fillHistos( htopMass_minus_gen, m_topminus_gen ,  weight, myType, jetBin, Nsolns);
                    fillHistos( httMass_gen, tt_mass_gen ,  weight, myType, jetBin, Nsolns);
                    fillHistos( httpT_gen, tt_pT_gen ,  weight, myType, jetBin, Nsolns);
                    fillHistos( hlepChargeAsym_gen, lep_charge_asymmetry_gen ,  weight, myType, jetBin, Nsolns);
                    fillHistos( hlepAzimAsym_gen, lep_azimuthal_asymmetry_gen ,  weight, myType, jetBin, Nsolns);
                    fillHistos( hlepAzimAsym2_gen, lep_azimuthal_asymmetry2_gen ,  weight, myType, jetBin, Nsolns);
                    if (m_top > 0 || applyNoCuts)
                    {
                        fillHistos( htopSpinCorr_gen, top_spin_correlation_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosOpeningAngle_gen, lep_cos_opening_angle_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( htopCosTheta_gen, top_costheta_cms_gen   ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_gen, lepPlus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_gen, lepMinus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPlusCosTheta_gen, lepPlus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hlepMinusCosTheta_gen, lepMinus_costheta_cms_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hpseudorapiditydiff_gen, top_pseudorapiditydiff_cms_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiff_gen, top_rapiditydiff_cms_gen ,  weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiffMarco_gen, top_rapiditydiff_Marco_gen ,  weight, myType, jetBin, Nsolns);
                    }


                    if (m_top > 0 || applyNoCuts)
                    {
                        //2D unfolding requires ttbar solution even for the purely leptonic variables
                        fillHistos( hlepChargeAsym_gen2d, lep_charge_asymmetry_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym_gen2d, lep_azimuthal_asymmetry_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym2_gen2d, lep_azimuthal_asymmetry2_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( htopSpinCorr_gen2d, top_spin_correlation_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosOpeningAngle_gen2d, lep_cos_opening_angle_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( htopCosTheta_gen2d, top_costheta_cms_gen   ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_gen2d, lepPlus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_gen2d, lepMinus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPlusCosTheta_gen2d, lepPlus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepMinusCosTheta_gen2d, lepMinus_costheta_cms_gen  ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hpseudorapiditydiff_gen2d, top_pseudorapiditydiff_cms_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiff_gen2d, top_rapiditydiff_cms_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiffMarco_gen2d, top_rapiditydiff_Marco_gen ,  tt_mass_gen, weight, myType, jetBin, Nsolns);

                        fillHistos( hlepChargeAsym_ttpT_gen2d, lep_charge_asymmetry_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym_ttpT_gen2d, lep_azimuthal_asymmetry_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym2_ttpT_gen2d, lep_azimuthal_asymmetry2_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( htopSpinCorr_ttpT_gen2d, top_spin_correlation_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosOpeningAngle_ttpT_gen2d, lep_cos_opening_angle_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( htopCosTheta_ttpT_gen2d, top_costheta_cms_gen   ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_ttpT_gen2d, lepPlus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_ttpT_gen2d, lepMinus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPlusCosTheta_ttpT_gen2d, lepPlus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hlepMinusCosTheta_ttpT_gen2d, lepMinus_costheta_cms_gen  ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hpseudorapiditydiff_ttpT_gen2d, top_pseudorapiditydiff_cms_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiff_ttpT_gen2d, top_rapiditydiff_cms_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiffMarco_ttpT_gen2d, top_rapiditydiff_Marco_gen ,  tt_pT_gen, weight, myType, jetBin, Nsolns);

                        fillHistos( hlepChargeAsym_ttRapidity2_gen2d, lep_charge_asymmetry_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym_ttRapidity2_gen2d, lep_azimuthal_asymmetry_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hlepAzimAsym2_ttRapidity2_gen2d, lep_azimuthal_asymmetry2_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( htopSpinCorr_ttRapidity2_gen2d, top_spin_correlation_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosOpeningAngle_ttRapidity2_gen2d, lep_cos_opening_angle_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( htopCosTheta_ttRapidity2_gen2d, top_costheta_cms_gen   ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_ttRapidity2_gen2d, lepPlus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hlepCosTheta_ttRapidity2_gen2d, lepMinus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hlepPlusCosTheta_ttRapidity2_gen2d, lepPlus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hlepMinusCosTheta_ttRapidity2_gen2d, lepMinus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hpseudorapiditydiff_ttRapidity2_gen2d, top_pseudorapiditydiff_cms_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiff_ttRapidity2_gen2d, top_rapiditydiff_cms_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                        fillHistos( hrapiditydiffMarco_ttRapidity2_gen2d, top_rapiditydiff_Marco_gen ,  abs(ttRapidity2_gen), weight, myType, jetBin, Nsolns);
                    }

                    if (from_gluon == true)
                    {
                        fillHistos( httMassGluongenp, tt_mass_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( httRapidityGluongenp, ttRapidity_gen  ,  weight, myType, jetBin, Nsolns);
                    }
                    else
                    {
                        fillHistos( httMassQuarkgenp, tt_mass_gen  ,  weight, myType, jetBin, Nsolns);
                        fillHistos( httRapidityQuarkgenp, ttRapidity_gen  ,  weight, myType, jetBin, Nsolns);
                    }

                }//only for mc

                //} //jet smearing loop, no longer used

            }//good hypothesis loop

        } // closes loop over events

        if (applyNoCuts) cout << "number of events (no cuts) before and after vertex weighting              =  " << nEvents_noCuts_novtxweight << "   " << nEvents_noCuts << endl;
        if (applyNoCuts) cout << "number of dilepton events (no cuts) before and after vertex weighting              =  " << nEvents_noCuts_novtxweight_dil << "   " << nEvents_noCuts_dil << endl;

        //float nEvents_primary = cms2.evt_nEvts();
        //cout << "acceptance                       =  " << (1.0*nSelectedEvents)/(nEvents_primary*kFactor * evt_scale1fb() * lumi) <<endl;

    }  // closes loop over files

    return;

} // closes myLooper function



void topAFB_looper::fillUnderOverFlow(TH1D *h1, float value, double weight, int Nsolns)
{
    double min = h1->GetXaxis()->GetXmin();
    double max = h1->GetXaxis()->GetXmax();

    if (value >= max) value = h1->GetBinCenter(h1->GetNbinsX());
    if (value <= min) value = h1->GetBinCenter(1);

    int bin_number = h1->FindBin(value);
    double orig_content = h1->GetBinContent(bin_number);
    double orig_error = h1->GetBinError(bin_number);

    //h1->Fill(value, weight);
    h1->SetBinContent( bin_number, orig_content + weight );
    h1->SetBinError( bin_number, sqrt( orig_error * orig_error + weight * weight * double(Nsolns) ) );
}

//--------------------------------------------------------------------

void topAFB_looper::fillUnderOverFlow(TH2D *h2, float xvalue, float yvalue, double weight, int Nsolns)
{
    double maxx = h2->GetXaxis()->GetXmax();
    double minx = h2->GetXaxis()->GetXmin();
    double maxy = h2->GetYaxis()->GetXmax();
    double miny = h2->GetYaxis()->GetXmin();

    if (xvalue >= maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
    if (xvalue <= minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
    if (yvalue >= maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
    if (yvalue <= miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

    int bin_number = h2->FindBin(xvalue, yvalue);
    double orig_content = h2->GetBinContent(bin_number);
    double orig_error = h2->GetBinError(bin_number);

    //h2->Fill(xvalue, yvalue, weight);
    h2->SetBinContent( bin_number, orig_content + weight );
    h2->SetBinError( bin_number, sqrt( orig_error * orig_error + weight * weight * double(Nsolns) ) );
}


//--------------------------------------------------------------------

void topAFB_looper::fillOverFlow(TH1D *h1, float value, float weight)
{
    float max = h1->GetXaxis()->GetXmax();
    if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
    h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void topAFB_looper::fillOverFlow(TH2D *h2, float xvalue, float yvalue, float weight)
{
    float maxx = h2->GetXaxis()->GetXmax();
    float maxy = h2->GetYaxis()->GetXmax();

    if (xvalue > maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
    if (yvalue > maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());

    h2->Fill(xvalue, yvalue, weight);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TH1D *h1[4][4], float value, double weight, int myType, int nJetsIdx, int Nsolns)
{
    //fillUnderOverFlow(h1[myType][nJetsIdx], value, weight, Nsolns);
    fillUnderOverFlow(h1[myType][3],        value, weight, Nsolns);
    //fillUnderOverFlow(h1[3][nJetsIdx],      value, weight, Nsolns);
    fillUnderOverFlow(h1[3][3],             value, weight, Nsolns);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TH2D *h2[4][4], float xvalue, float yvalue, double weight, int myType, int nJetsIdx, int Nsolns)
{
    //fillUnderOverFlow(h2[myType][nJetsIdx], xvalue, yvalue, weight, Nsolns);
    fillUnderOverFlow(h2[myType][3],        xvalue, yvalue, weight, Nsolns);
    //fillUnderOverFlow(h2[3][nJetsIdx],      xvalue, yvalue, weight, Nsolns);
    fillUnderOverFlow(h2[3][3],             xvalue, yvalue, weight, Nsolns);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TProfile *h2[4][4], float xvalue, float yvalue, int myType, int nJetsIdx)
{
    //h2[myType][nJetsIdx] -> Fill(xvalue, yvalue);
    h2[myType][3]        -> Fill(xvalue, yvalue);
    //h2[3][nJetsIdx]      -> Fill(xvalue, yvalue);
    h2[3][3]             -> Fill(xvalue, yvalue);
}

