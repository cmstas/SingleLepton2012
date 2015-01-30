#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <exception>
#include <set>
#include <algorithm>
#include <map>

#include "LHAPDF/LHAPDF.h"
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
#include "./topAFB_looper.h"
#include "Histograms.cc"

using namespace std;
using namespace tas;

typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;


double topAFB_looper::TopPtWeight(double topPt)
{
    if ( topPt < 0 ) return 1;
    if (topPt > 400) topPt = 400;

    double result = exp(0.156 - 0.00137 * topPt); //8 TeV fit (l+j and 2l combined)
    //note this fit is for data/madgraph, and we are using MC@NLO

    return result;
}



//--------------------------------------------------------------------

struct DorkyStatus1Identifier {
  int id;
  float px, py, pz, E;
  bool operator < (const DorkyStatus1Identifier &) const;
  bool operator == (const DorkyStatus1Identifier &) const;
};

//--------------------------------------------------------------------
/*
bool DorkyStatus1Identifier::operator < (const DorkyStatus1Identifier &other) const
{
  if (id != other.id)
    return id < other.id;
  if (px != other.px)
    return px < other.px;
  if (py != other.py)
    return py < other.py;
  if (pz != other.pz)
    return pz < other.pz;
  if (E != other.E)
    return E < other.E;
  return false;
}
*/

// != works because the floats are identical, but this is safer
bool DorkyStatus1Identifier::operator < (const DorkyStatus1Identifier &other) const
{
  if (id != other.id)
    return id < other.id;
  if (fabs(1. - px/other.px)>1e-6)
    return px < other.px;
  if (fabs(1. - py/other.py)>1e-6)
    return py < other.py;
  if (fabs(1. - pz/other.pz)>1e-6)
    return pz < other.pz;
  if (fabs(1. - E/other.E)>1e-6)
    return E < other.E;
  return false;
}

//--------------------------------------------------------------------
/*
bool DorkyStatus1Identifier::operator == (const DorkyStatus1Identifier &other) const
{
  if (id != other.id)
    return false;
  if (px != other.px)
    return false;
  if (py != other.py)
    return false;
  if (pz != other.pz)
    return false;
  if (E != other.E)
    return false;
  return true;
}
*/

// != works because the floats are identical, but this is safer
bool DorkyStatus1Identifier::operator == (const DorkyStatus1Identifier &other) const
{
  if (id != other.id)
    return false;
  if (fabs(1. - px/other.px)>1e-6)
    return false;
  if (fabs(1. - py/other.py)>1e-6)
    return false;
  if (fabs(1. - pz/other.pz)>1e-6)
    return false;
  if (fabs(1. - E/other.E)>1e-6)
    return false;
  return true;
}


//--------------------------------------------------------------------

std::set<DorkyStatus1Identifier> already_seen_stat1;
bool is_duplicate_stat1 (const DorkyStatus1Identifier &id) {
  std::pair<std::set<DorkyStatus1Identifier>::const_iterator, bool> ret =
    already_seen_stat1.insert(id);
  return !ret.second;
}

std::set<DorkyStatus1Identifier> already_seen_stat1_t;
bool is_duplicate_stat1_t (const DorkyStatus1Identifier &id) {
  std::pair<std::set<DorkyStatus1Identifier>::const_iterator, bool> ret =
    already_seen_stat1_t.insert(id);
  return !ret.second;
}

std::set<DorkyStatus1Identifier> already_seen_stat1_b;
bool is_duplicate_stat1_b (const DorkyStatus1Identifier &id) {
  std::pair<std::set<DorkyStatus1Identifier>::const_iterator, bool> ret =
    already_seen_stat1_b.insert(id);
  return !ret.second;
}

std::set<DorkyStatus1Identifier> already_seen_stat2_bhadron;
bool is_duplicate_stat2_bhadron (const DorkyStatus1Identifier &id) {
  std::pair<std::set<DorkyStatus1Identifier>::const_iterator, bool> ret =
    already_seen_stat2_bhadron.insert(id);
  return !ret.second;
}

//--------------------------------------------------------------------






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
    TString prefixT = prefix;

    bookHistos(prefix.c_str(), 1, 1);
    bool isData = false;
    bool applyNoCuts = true;

    LHAPDF::initPDFSetM(1,"../pdfs/cteq6mE.LHgrid");
    LHAPDF::initPDFM(1, 0);
    LHAPDF::initPDFSetM(2,"../pdfs/cteq6mE.LHgrid");

    //--------------------------
    // File and Event Loop
    //---------------------------
    TObjArray *listOfFiles = chain->GetListOfFiles();
    unsigned int nEventsChain = 0;
    unsigned int nEvents = chain->GetEntries();
    nEventsChain = nEvents;
    unsigned int nEventsTotal = 0;
    ULong64_t nEvents_noCuts = 0;
    ULong64_t nEvents_noCuts_dil = 0;
    ULong64_t nEvents_noCuts_nonegweight = 0;
    ULong64_t nEvents_noCuts_nonegweight_dil = 0;
    double nSelectedEvents = 0;
    unsigned int npreSelectedEvents = 0;
    unsigned int npreSelectedEvents_genmatch1 = 0;
    unsigned int npreSelectedEvents_genmatch2 = 0;
    TIter fileIter(listOfFiles);
    map<int, int> m_events;
    while (TChainElement *currentFile = (TChainElement *)fileIter.Next() )
    {

        TString filename = currentFile->GetTitle();

        bool ismcatnlo = false;

        if(filename.Contains("mcatnlo") || prefixT.Contains("mcatnlo") ) {
          cout<<"Processing MC@NLO sample: using status=2 Ws (before FSR from b)"<<endl;
          ismcatnlo = true; 
        }

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
            int ndavtx = 0;
            double weight = 1.0;
            double weights[41];
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
                nleps = leptonGenpCount_lepTauDecays_status3only(nels, nmus, ntaus);

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
                    if (fillEvent) fillHistos( hlepPlusCosThetaTau_gen,  cosTheta_lepPlus_status1, weighttaudecay ? weight_lepPlus : 1., myType, 0);
                    if (fillEvent) fillHistos( hlepPlusxTau_gen,  EoverEmax_lepPlus, weighttaudecay ? weight_lepPlus : 1., myType, 0);
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
                    if (fillEvent) fillHistos( hlepMinusCosThetaTau_gen,  cosTheta_lepMinus_status1, weighttaudecay ? weight_lepMinus : 1., myType, 0);
                    if (fillEvent) fillHistos( hlepMinusxTau_gen,  EoverEmax_lepMinus, weighttaudecay ? weight_lepMinus : 1., myType, 0);
                    weight_taudecay *= weight_lepMinus;
                }

            } //ntaus>0
            //cout<<weight_taudecay<<endl;

            if (!isData) weight = evt_scale1fb() * lumi;
            //negative weights for MC@NLO
            if (prefix == "ttdil" || prefix == "ttotr") weight = weight * (fabs(genps_weight()) / genps_weight());
            //tau decay cosTheta* weighting
            if (weighttaudecay && (prefix == "ttdil" || prefix == "ttotr")  && ntaus > 0) weight *= weight_taudecay;

            nEvents_noCuts_nonegweight += 1;
            if (genps_weight()>0) nEvents_noCuts += 1;
            else if (genps_weight()<0) nEvents_noCuts -= 1;
            else cout<<"ambiguous genps_weight()"<<endl;

            if (prefix == "ttdil"    &&  nleps != 2) continue;
            if (prefix == "ttotr"    &&  nleps == 2) continue;
            if (prefix == "DYee"     &&  nels != 2) continue;
            if (prefix == "DYmm"     &&  nmus != 2) continue;
            if (prefix == "DYtautau" &&  ntaus != 2) continue;

            nEvents_noCuts_nonegweight_dil += 1;
            if (genps_weight()>0) nEvents_noCuts_dil += 1;
            else if (genps_weight()<0) nEvents_noCuts_dil -= 1;
            else cout<<"ambiguous genps_weight() in dilepton event"<<endl;


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

                float m_topminus_gen;
                float m_topplus_gen;
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



/* //replaced by fillgenlevel
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

                    m_topminus_gen = topminus_genp_p4.M();
                    m_topplus_gen = topplus_genp_p4.M();

                    tt_mass_gen = (topplus_genp_p4 + topminus_genp_p4).M();
                    ttRapidity_gen = topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity();
                    ttRapidity2_gen = (topplus_genp_p4 + topminus_genp_p4).Rapidity();

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

*/

                    int ntops = 0;

                    for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) if( abs( genps_id().at(igen) ) == 6 && genps_status().at(igen) == 3 ) ntops++;

                    if ( ntops != 2 ) { cout<<"*** skipping event with ntops = "<<ntops<<" ***"<<endl; continue; } //only events with 2 tops count as ttdl

                    fillgenlevel(ismcatnlo, nleps, ntaus, ntops);

                    if ( lepPlus_status1_id_ == -11 && lepMinus_status1_id_ == 11 ) myType = 0;
                    else if ( lepPlus_status1_id_ == -13 && lepMinus_status1_id_ == 11 ) myType = 2;
                    else if ( lepPlus_status1_id_ == -11 && lepMinus_status1_id_ == 13 ) myType = 2;
                    else if ( lepPlus_status1_id_ == -13 && lepMinus_status1_id_ == 13 ) myType = 1;
                    else cout<<"indeterminate dilepton type: "<<lepPlus_status1_id_<<" "<<lepMinus_status1_id_<<endl;

/* confirmed fillgenlevel gives same results
                    if(m_topminus_gen-m_topminus_gen_origleps_!=0) cout<<m_topminus_gen-m_topminus_gen_origleps_<<endl;
                    if(m_topplus_gen-m_topplus_gen_origleps_!=0) cout<<m_topplus_gen-m_topplus_gen_origleps_<<endl;
                    if(tt_mass_gen-tt_mass_gen_origleps_!=0) cout<<tt_mass_gen-tt_mass_gen_origleps_<<endl;
                    if(ttRapidity_gen-ttRapidity_gen_origleps_!=0) cout<<ttRapidity_gen-ttRapidity_gen_origleps_<<endl;
                    if(ttRapidity2_gen-ttRapidity2_gen_origleps_!=0) cout<<ttRapidity2_gen-ttRapidity2_gen_origleps_<<endl;
                    if(top_rapiditydiff_cms_gen-top_rapiditydiff_cms_gen_origleps_!=0) cout<<top_rapiditydiff_cms_gen-top_rapiditydiff_cms_gen_origleps_<<endl;
                    if(top_pseudorapiditydiff_cms_gen-top_pseudorapiditydiff_cms_gen_origleps_!=0) cout<<top_pseudorapiditydiff_cms_gen-top_pseudorapiditydiff_cms_gen_origleps_<<endl;
                    if(top_rapiditydiff_Marco_gen-top_rapiditydiff_Marco_gen_origleps_!=0) cout<<top_rapiditydiff_Marco_gen-top_rapiditydiff_Marco_gen_origleps_<<endl;
                    if(tt_pT_gen-tt_pT_gen_origleps_!=0) cout<<tt_pT_gen-tt_pT_gen_origleps_<<endl;
                    if(top_costheta_cms_gen-top_costheta_cms_gen_origleps_!=0) cout<<top_costheta_cms_gen-top_costheta_cms_gen_origleps_<<endl;
                    if(lep_charge_asymmetry_gen-lep_charge_asymmetry_gen_origleps_!=0) cout<<lep_charge_asymmetry_gen-lep_charge_asymmetry_gen_origleps_<<endl;
                    if(lep_azimuthal_asymmetry_gen-lep_azimuthal_asymmetry_gen_origleps_!=0) cout<<lep_azimuthal_asymmetry_gen-lep_azimuthal_asymmetry_gen_origleps_<<endl;
                    if(lep_azimuthal_asymmetry2_gen-lep_azimuthal_asymmetry2_gen_origleps_!=0) cout<<lep_azimuthal_asymmetry2_gen-lep_azimuthal_asymmetry2_gen_origleps_<<endl;
                    if(lepPlus_costheta_cms_gen-lepPlus_costheta_cms_gen_origleps_!=0) cout<<lepPlus_costheta_cms_gen-lepPlus_costheta_cms_gen_origleps_<<endl;
                    if(lepMinus_costheta_cms_gen-lepMinus_costheta_cms_gen_origleps_!=0) cout<<lepMinus_costheta_cms_gen-lepMinus_costheta_cms_gen_origleps_<<endl;
                    if(top_spin_correlation_gen-top_spin_correlation_gen_origleps_!=0) cout<<top_spin_correlation_gen-top_spin_correlation_gen_origleps_<<endl;
                    if(lep_cos_opening_angle_gen-lep_cos_opening_angle_gen_origleps_!=0) cout<<lep_cos_opening_angle_gen-lep_cos_opening_angle_gen_origleps_<<endl;
*/


                    m_topminus_gen = m_topminus_gen_;
                    m_topplus_gen = m_topplus_gen_;
                    tt_mass_gen = tt_mass_gen_;
                    ttRapidity_gen = ttRapidity_gen_;
                    ttRapidity2_gen = ttRapidity2_gen_;
                    top_rapiditydiff_cms_gen = top_rapiditydiff_cms_gen_;
                    top_pseudorapiditydiff_cms_gen = top_pseudorapiditydiff_cms_gen_;
                    top_rapiditydiff_Marco_gen = top_rapiditydiff_Marco_gen_;
                    tt_pT_gen = tt_pT_gen_;
                    top_costheta_cms_gen = top_costheta_cms_gen_;
                    lep_charge_asymmetry_gen = lep_charge_asymmetry_gen_;
                    lep_azimuthal_asymmetry_gen = lep_azimuthal_asymmetry_gen_;
                    lep_azimuthal_asymmetry2_gen = lep_azimuthal_asymmetry2_gen_;
                    lepPlus_costheta_cms_gen = lepPlus_costheta_cms_gen_;
                    lepMinus_costheta_cms_gen = lepMinus_costheta_cms_gen_;
                    top_spin_correlation_gen = top_spin_correlation_gen_;
                    lep_cos_opening_angle_gen = lep_cos_opening_angle_gen_;


                    // Calculate PDF weights
                    if ( prefix == "ttdil"  ){


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
                          weights[subset] = ((fx1Qi*fx2Qi)/(fx1Q0gen*fx2Q0gen));
                          //cout<<subset<<" "<<weights[subset]<<endl;
                          weights[subset] *= weight;
                          //cout << fx1Qi <<endl;
                          // cout << fx1Q0gen <<endl;
                          // cout << pdf_weight <<endl;


                        }// end of loop over subset of PDFs
                    }



                    for (unsigned int PDFset = 0; PDFset < 41; PDFset++)
                    {

                        fillHistos( htopMass_plus_gen, m_topplus_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                        fillHistos( htopMass_minus_gen, m_topminus_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                        fillHistos( httMass_gen, tt_mass_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                        fillHistos( httpT_gen, tt_pT_gen ,  weights[PDFset], myType, PDFset, Nsolns);

                        fillHistos( hlepChargeAsym_gen, lep_charge_asymmetry_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                        fillHistos( hlepAzimAsym_gen, lep_azimuthal_asymmetry_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                        fillHistos( hlepAzimAsym2_gen, lep_azimuthal_asymmetry2_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                        if (m_top > 0 || applyNoCuts)
                        {
                            fillHistos( htopSpinCorr_gen, top_spin_correlation_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosOpeningAngle_gen, lep_cos_opening_angle_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( htopCosTheta_gen, top_costheta_cms_gen   ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_gen, lepPlus_costheta_cms_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_gen, lepMinus_costheta_cms_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepPlusCosTheta_gen, lepPlus_costheta_cms_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepMinusCosTheta_gen, lepMinus_costheta_cms_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hpseudorapiditydiff_gen, top_pseudorapiditydiff_cms_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiff_gen, top_rapiditydiff_cms_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiffMarco_gen, top_rapiditydiff_Marco_gen ,  weights[PDFset], myType, PDFset, Nsolns);
                        }


                        if (m_top > 0 || applyNoCuts)
                        {
                            //2D unfolding requires ttbar solution even for the purely leptonic variables
                            fillHistos( hlepChargeAsym_gen2d, lep_charge_asymmetry_gen ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepAzimAsym_gen2d, lep_azimuthal_asymmetry_gen ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepAzimAsym2_gen2d, lep_azimuthal_asymmetry2_gen ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( htopSpinCorr_gen2d, top_spin_correlation_gen  ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosOpeningAngle_gen2d, lep_cos_opening_angle_gen  ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( htopCosTheta_gen2d, top_costheta_cms_gen   ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_gen2d, lepPlus_costheta_cms_gen  ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_gen2d, lepMinus_costheta_cms_gen  ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepPlusCosTheta_gen2d, lepPlus_costheta_cms_gen  ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepMinusCosTheta_gen2d, lepMinus_costheta_cms_gen  ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hpseudorapiditydiff_gen2d, top_pseudorapiditydiff_cms_gen ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiff_gen2d, top_rapiditydiff_cms_gen ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiffMarco_gen2d, top_rapiditydiff_Marco_gen ,  tt_mass_gen, weights[PDFset], myType, PDFset, Nsolns);

                            fillHistos( hlepChargeAsym_ttpT_gen2d, lep_charge_asymmetry_gen ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepAzimAsym_ttpT_gen2d, lep_azimuthal_asymmetry_gen ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepAzimAsym2_ttpT_gen2d, lep_azimuthal_asymmetry2_gen ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( htopSpinCorr_ttpT_gen2d, top_spin_correlation_gen  ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosOpeningAngle_ttpT_gen2d, lep_cos_opening_angle_gen  ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( htopCosTheta_ttpT_gen2d, top_costheta_cms_gen   ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_ttpT_gen2d, lepPlus_costheta_cms_gen  ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_ttpT_gen2d, lepMinus_costheta_cms_gen  ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepPlusCosTheta_ttpT_gen2d, lepPlus_costheta_cms_gen  ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepMinusCosTheta_ttpT_gen2d, lepMinus_costheta_cms_gen  ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hpseudorapiditydiff_ttpT_gen2d, top_pseudorapiditydiff_cms_gen ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiff_ttpT_gen2d, top_rapiditydiff_cms_gen ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiffMarco_ttpT_gen2d, top_rapiditydiff_Marco_gen ,  tt_pT_gen, weights[PDFset], myType, PDFset, Nsolns);

                            fillHistos( hlepChargeAsym_ttRapidity2_gen2d, lep_charge_asymmetry_gen ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepAzimAsym_ttRapidity2_gen2d, lep_azimuthal_asymmetry_gen ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepAzimAsym2_ttRapidity2_gen2d, lep_azimuthal_asymmetry2_gen ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( htopSpinCorr_ttRapidity2_gen2d, top_spin_correlation_gen  ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosOpeningAngle_ttRapidity2_gen2d, lep_cos_opening_angle_gen  ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( htopCosTheta_ttRapidity2_gen2d, top_costheta_cms_gen   ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_ttRapidity2_gen2d, lepPlus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepCosTheta_ttRapidity2_gen2d, lepMinus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepPlusCosTheta_ttRapidity2_gen2d, lepPlus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hlepMinusCosTheta_ttRapidity2_gen2d, lepMinus_costheta_cms_gen  ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hpseudorapiditydiff_ttRapidity2_gen2d, top_pseudorapiditydiff_cms_gen ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiff_ttRapidity2_gen2d, top_rapiditydiff_cms_gen ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( hrapiditydiffMarco_ttRapidity2_gen2d, top_rapiditydiff_Marco_gen ,  abs(ttRapidity2_gen), weights[PDFset], myType, PDFset, Nsolns);
                        }
                        /*
                        if (from_gluon == true)
                        {
                            fillHistos( httMassGluongenp, tt_mass_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( httRapidityGluongenp, ttRapidity_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                        }
                        else
                        {
                            fillHistos( httMassQuarkgenp, tt_mass_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                            fillHistos( httRapidityQuarkgenp, ttRapidity_gen  ,  weights[PDFset], myType, PDFset, Nsolns);
                        }
                        */

                    }

                }//only for mc

                //} //jet smearing loop, no longer used

            }//good hypothesis loop

        } // closes loop over events

        if (applyNoCuts) cout << "number of events (no cuts) with and without accounting for negative weights              =  " << nEvents_noCuts_nonegweight << "   " << nEvents_noCuts << endl;
        if (applyNoCuts) cout << "number of dilepton events (no cuts) with and without accounting for negative weights              =  " << nEvents_noCuts_nonegweight_dil << "   " << nEvents_noCuts_dil << endl;

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

void topAFB_looper::fillHistos(TH1D *h1[4][41], float value, double weight, int myType, int PDFset, int Nsolns)
{
    fillUnderOverFlow(h1[myType][PDFset], value, weight, Nsolns);
    fillUnderOverFlow(h1[3][PDFset],      value, weight, Nsolns);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TH2D *h2[4][41], float xvalue, float yvalue, double weight, int myType, int PDFset, int Nsolns)
{
    fillUnderOverFlow(h2[myType][PDFset], xvalue, yvalue, weight, Nsolns);
    fillUnderOverFlow(h2[3][PDFset],      xvalue, yvalue, weight, Nsolns);
}

//--------------------------------------------------------------------

void topAFB_looper::fillHistos(TProfile *h2[4][41], float xvalue, float yvalue, int myType, int PDFset)
{
    h2[myType][PDFset] -> Fill(xvalue, yvalue);
    h2[3][PDFset]      -> Fill(xvalue, yvalue);
}


int topAFB_looper::leptonGenpCount_lepTauDecays_status3only(int &nele, int &nmuon, int &ntau){
    nele = 0;
    nmuon = 0;
    ntau = 0;
    int size = cms2.genps_id().size();
    for (int jj = 0; jj < size; jj++)
    {
        if ( cms2.genps_status().at(jj) != 3 ) continue;
        if (abs(cms2.genps_id().at(jj)) == 11) nele++;
        if (abs(cms2.genps_id().at(jj)) == 13) nmuon++;
        if (abs(cms2.genps_id().at(jj)) == 15)
        {
            for (unsigned int kk = 0; kk < cms2.genps_lepdaughter_id()[jj].size(); kk++)
            {
                int daughter = abs(cms2.genps_lepdaughter_id()[jj][kk]);
                //we count neutrino's because that guarantees that
                //there is a corresponding lepton and that it comes from
                // a leptonic tau decay. You can get electrons from converted photons
                //which are radiated by charged pions from the tau decay but thats
                //hadronic and we don't care for those
                if ( daughter == 12 || daughter == 14)
                    ntau++;
            }//daughter loop
        }//if tau
    }//genps loop

    return nele + nmuon + ntau;
}









//below this line is a direct copy and paste from singleLeptonLooper.cc (with mcmtln_ and mcmln_ commented out)

void topAFB_looper::fillgenlevel(bool ismcatnlo, int nleps, int ntaus, int ntops) {

  TLorentzVector WPlus_status2_T(0, 0, 0, 0), WMinus_status2_T(0, 0, 0, 0);
  LorentzVector WPlus_status2(0, 0, 0, 0), WMinus_status2(0, 0, 0, 0);
  LorentzVector topPlus_status2(0, 0, 0, 0), topMinus_status2(0, 0, 0, 0);

  vt_stat1.SetXYZT(0,0,0,0);
  vtbar_stat1.SetXYZT(0,0,0,0);
  vb_stat1.SetXYZT(0,0,0,0);
  vbbar_stat1.SetXYZT(0,0,0,0);
  WPlus_status3_orig.SetXYZT(0,0,0,0);
  WMinus_status3_orig.SetXYZT(0,0,0,0);
  topPlus_status3_bW.SetXYZT(0,0,0,0);
  topMinus_status3_bW.SetXYZT(0,0,0,0);
  lepPlus_status3_orig.SetXYZT(0,0,0,0);
  lepMinus_status3_orig.SetXYZT(0,0,0,0);
  nuPlus_status3_orig.SetXYZT(0,0,0,0);
  nuMinus_status3_orig.SetXYZT(0,0,0,0);
  ttpair.SetXYZT(0,0,0,0);

  //int ntops = 0;
  int foundstat2WPlus = 0;
  int foundstat2WMinus = 0;
  int foundstat3WPlus = 0;
  int foundstat3WMinus = 0;
  int foundbPlus = 0;
  int foundbMinus = 0;
  int foundtopPlus = 0;
  int foundtopMinus = 0;
  int ntopPlusDaughters = 0;
  int ntopMinusDaughters = 0;

  bool foundwz = false;
  bool foundlep = false;
  bool foundnu  = false;

  //vector <LorentzVector> t_daughters;
  //vector <LorentzVector> tbar_daughters;
  int t_dup_stat1 = 0;
  int tbar_dup_stat1 = 0;

  //vector <LorentzVector> b_daughters;
  //vector <LorentzVector> bbar_daughters;
  int b_dup_stat1 = 0;
  int bbar_dup_stat1 = 0;

  int n_bhadrons = 0;

  already_seen_stat1_t.clear();
  already_seen_stat1_b.clear();
  already_seen_stat2_bhadron.clear();

  for ( int igen = 0 ; igen < (int)genps_id().size() ; igen++ ) { 

    int id = genps_id().at(igen);
    int pid = abs( genps_id().at(igen) );
    int mothid = abs(genps_id_mother().at(igen));

    //find final status=2 stable b-hadrons for use in genjet clustering
    if( genps_status().at(igen) == 2 && isBHadron(igen) ) {
      DorkyStatus1Identifier ident = { cms2.genps_id().at(igen), cms2.genps_p4().at(igen).Px(), cms2.genps_p4().at(igen).Py(), cms2.genps_p4().at(igen).Pz(), cms2.genps_p4().at(igen).E() };
      if(  !is_duplicate_stat2_bhadron(ident) ) {
        n_bhadrons++;
        //cout<<genps_status().at(igen)<<" "<<n_bhadrons<<" "<<genps_p4().at(igen).Pt()<<" "<<id<<" "<<genps_id_mother().at(igen)<<endl;
      }
    }


    //For MC@NLO. Find first status=2 Ws, which are before FSR. Note, unlike the W the status=3 b is before FSR.
    if( id == 24 ){
      if( genps_status().at(igen) == 2 && !foundstat2WPlus ) {
        foundstat2WPlus = 1;
        WPlus_status2.SetXYZT(genps_p4()[igen].x(),
                            genps_p4()[igen].y(),
                            genps_p4()[igen].z(),
                            genps_p4()[igen].t()
                           );
        WPlus_status2_T.SetXYZT(genps_p4()[igen].x(),
                            genps_p4()[igen].y(),
                            genps_p4()[igen].z(),
                            genps_p4()[igen].t()
                           );
      }
      if( genps_status().at(igen) == 3 && !foundstat3WPlus ) {
        foundstat3WPlus = 1;
        WPlus_status3_ = &(genps_p4().at(igen));
        WPlus_status3_orig = genps_p4().at(igen);
        WPlus_status3_orig_ = &(WPlus_status3_orig);
      }
    }
    if( id == -24 ){
      if( genps_status().at(igen) == 2 && !foundstat2WMinus ) {
        foundstat2WMinus = 1;
        WMinus_status2.SetXYZT(genps_p4()[igen].x(),
                            genps_p4()[igen].y(),
                            genps_p4()[igen].z(),
                            genps_p4()[igen].t()
                           );
        WMinus_status2_T.SetXYZT(genps_p4()[igen].x(),
                            genps_p4()[igen].y(),
                            genps_p4()[igen].z(),
                            genps_p4()[igen].t()
                           );
      }
      if( genps_status().at(igen) == 3 && !foundstat3WMinus ) {
        foundstat3WMinus = 1;
        WMinus_status3_ = &(genps_p4().at(igen));
        WMinus_status3_orig = genps_p4().at(igen);
        WMinus_status3_orig_ = &(WMinus_status3_orig);
      }
    }

    //take last status=2 top (for testing only)
    if( id == 6 && genps_status().at(igen) == 2 ){
        topPlus_status2.SetXYZT(genps_p4()[igen].x(),
                            genps_p4()[igen].y(),
                            genps_p4()[igen].z(),
                            genps_p4()[igen].t()
                           );
    }
    if( id == -6 && genps_status().at(igen) == 2 ){
        topMinus_status2.SetXYZT(genps_p4()[igen].x(),
                            genps_p4()[igen].y(),
                            genps_p4()[igen].z(),
                            genps_p4()[igen].t()
                           );
    }


    if( genps_status().at(igen) != 3 ) continue;

    // only count/store massive b quarks
    // - in massless b samples, final state b's will have the mass set to 4.8 GeV
    //    while initial state b's have mass 0
    if ( abs(id) == 5 ) {
      // mass (from dumpDocLines)
      float m2 = cms2.genps_p4().at(igen).M2();
      float m  = m2 >= 0 ? sqrt(m2) : 0.0;
      if (m > 0.) {
        ++nbs_;
        genbs_.push_back(genps_p4().at(igen));
      }
    }




    if ( abs(genps_id_mother()[igen]) == 6 ) {
      //cout<<id<<" ";
      if ( (genps_id_mother()[igen]) == 6 ) ntopPlusDaughters++;
      if ( (genps_id_mother()[igen]) == -6 ) ntopMinusDaughters++;
      //if(fabs(id)==21) cout<<genps_id_mother()[igen]<<" gluon daughter"<<endl;  //MC@NLO has tops with status=3 gluon daughters
      if( ( id == 5 || (ismcatnlo && (id == 3 || id == 1)) ) && !foundbPlus ){   //id = 1,3 is for MC@NLO where some tops decay to dW/sW instead of bW
        foundbPlus = 1;
        b_         = &(genps_p4().at(igen));

        bPlus_status3_ = &(genps_p4().at(igen));

        //Create status=1 b. This only works properly for mc@nlo+herwig.
        already_seen_stat1.clear();
        for (unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++)
        {
            DorkyStatus1Identifier ident = { cms2.genps_lepdaughter_id()[igen][kk], cms2.genps_lepdaughter_p4()[igen][kk].Px(), cms2.genps_lepdaughter_p4()[igen][kk].Py(), cms2.genps_lepdaughter_p4()[igen][kk].Pz(), cms2.genps_lepdaughter_p4()[igen][kk].E() };
            if ( is_duplicate_stat1(ident) ) {b_dup_stat1++; continue;}
            if ( is_duplicate_stat1_b(ident) ) {cout<<"***********this should be impossible************"<<endl;}
            //b_daughters.push_back(cms2.genps_lepdaughter_p4()[igen][kk]);
            vb_stat1 += cms2.genps_lepdaughter_p4()[igen][kk];
        }
        ////cout<<nbs_<<" b_: "<<b_daughters.size()<<" "<<b_dup_stat1<<endl;
        //cout<<" "<<b_->Px()<<" "<<b_->Py()<<" "<<b_->Pz()<<" "<<b_->E()<<endl;
        //cout<<" "<<vb_stat1.Px()<<" "<<vb_stat1.Py()<<" "<<vb_stat1.Pz()<<" "<<vb_stat1.E()<<endl;

      }
      if( (id == -5 || (ismcatnlo && (id == -3 || id == -1)) ) && !foundbMinus ){ //id = -1,-3 is for MC@NLO where some tops decay to dW/sW instead of bW
        foundbMinus = 1;
        bbar_      = &(genps_p4().at(igen));

        bMinus_status3_ = &(genps_p4().at(igen));

        //Create status=1 bbar. This only works properly for mc@nlo+herwig.
        already_seen_stat1.clear();
        for (unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++)
        {
            DorkyStatus1Identifier ident = { cms2.genps_lepdaughter_id()[igen][kk], cms2.genps_lepdaughter_p4()[igen][kk].Px(), cms2.genps_lepdaughter_p4()[igen][kk].Py(), cms2.genps_lepdaughter_p4()[igen][kk].Pz(), cms2.genps_lepdaughter_p4()[igen][kk].E() };
            if ( is_duplicate_stat1(ident) ) {bbar_dup_stat1++; continue;}
            if ( is_duplicate_stat1_b(ident) && ismcatnlo ) {cout<<"***********bbar shares daughter with b************"<<endl; cout<<" "<<cms2.genps_lepdaughter_p4()[igen][kk].Px()<<" "<<cms2.genps_lepdaughter_p4()[igen][kk].Py()<<" "<<cms2.genps_lepdaughter_p4()[igen][kk].Pz()<<" "<<cms2.genps_lepdaughter_p4()[igen][kk].E()<<endl;}
            //bbar_daughters.push_back(cms2.genps_lepdaughter_p4()[igen][kk]);
            vbbar_stat1 += cms2.genps_lepdaughter_p4()[igen][kk];
        }
        ////cout<<nbs_<<" bbar_: "<<bbar_daughters.size()<<" "<<bbar_dup_stat1<<endl;
        //cout<<" "<<bbar_->Px()<<" "<<bbar_->Py()<<" "<<bbar_->Pz()<<" "<<bbar_->E()<<endl;
        //cout<<" "<<vbbar_stat1.Px()<<" "<<vbbar_stat1.Py()<<" "<<vbbar_stat1.Pz()<<" "<<vbbar_stat1.E()<<endl;

      }
    }



    if( id == 6 && !foundtopPlus ){
      foundtopPlus = 1;
      t_         = &(genps_p4().at(igen));
      ptt_       = genps_p4().at(igen).pt();
      //ntops++;

      topPlus_status3_ = &(genps_p4().at(igen));
      topPlus_status3_bW = genps_p4().at(igen);
      topPlus_status3_bW_ = &(topPlus_status3_bW);

      //Create status=1 top. This only works properly for mc@nlo+herwig.
      already_seen_stat1.clear();
      for (unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++)
      {
          DorkyStatus1Identifier ident = { cms2.genps_lepdaughter_id()[igen][kk], cms2.genps_lepdaughter_p4()[igen][kk].Px(), cms2.genps_lepdaughter_p4()[igen][kk].Py(), cms2.genps_lepdaughter_p4()[igen][kk].Pz(), cms2.genps_lepdaughter_p4()[igen][kk].E() };
          if ( is_duplicate_stat1(ident) ) {t_dup_stat1++; continue;}
          if ( is_duplicate_stat1_t(ident) ) {cout<<"***********this should be impossible************"<<endl;}
          //t_daughters.push_back(cms2.genps_lepdaughter_p4()[igen][kk]);
          vt_stat1 += cms2.genps_lepdaughter_p4()[igen][kk];
      }
      ////cout<<"t_: "<<t_daughters.size()<<" "<<t_dup_stat1<<endl;
      //cout<<" "<<t_->Px()<<" "<<t_->Py()<<" "<<t_->Pz()<<" "<<t_->E()<<endl;
      //cout<<" "<<vt_stat1.Px()<<" "<<vt_stat1.Py()<<" "<<vt_stat1.Pz()<<" "<<vt_stat1.E()<<endl;

    }
    if( id == -6 && !foundtopMinus ){
      foundtopMinus = 1;
      tbar_      = &(genps_p4().at(igen));
      pttbar_    = genps_p4().at(igen).pt();
      //ntops++;

      topMinus_status3_ = &(genps_p4().at(igen));
      topMinus_status3_bW = genps_p4().at(igen);
      topMinus_status3_bW_ = &(topMinus_status3_bW);

      //Create status=1 tbar. This only works properly for mc@nlo+herwig.
      already_seen_stat1.clear();
      for (unsigned int kk = 0; kk < cms2.genps_lepdaughter_id().at(igen).size(); kk++)
      {
          DorkyStatus1Identifier ident = { cms2.genps_lepdaughter_id()[igen][kk], cms2.genps_lepdaughter_p4()[igen][kk].Px(), cms2.genps_lepdaughter_p4()[igen][kk].Py(), cms2.genps_lepdaughter_p4()[igen][kk].Pz(), cms2.genps_lepdaughter_p4()[igen][kk].E() };
          if ( is_duplicate_stat1(ident) ) {tbar_dup_stat1++; continue;}
          if ( is_duplicate_stat1_t(ident) && ismcatnlo ) {cout<<"***********tbar shares daughter with t************"<<endl;}
          //tbar_daughters.push_back(cms2.genps_lepdaughter_p4()[igen][kk]);
          vtbar_stat1 += cms2.genps_lepdaughter_p4()[igen][kk];
      }
      ////cout<<"tbar_: "<<tbar_daughters.size()<<" "<<tbar_dup_stat1<<endl;
      //cout<<" "<<tbar_->Px()<<" "<<tbar_->Py()<<" "<<tbar_->Pz()<<" "<<tbar_->E()<<endl;
      //cout<<" "<<vtbar_stat1.Px()<<" "<<vtbar_stat1.Py()<<" "<<vtbar_stat1.Pz()<<" "<<vtbar_stat1.E()<<endl;

    }

                                if ( genps_id_mother()[igen] == 24 )
                                {
                                    if ( (genps_id()[igen] == -11 || genps_id()[igen] == -13 ||  genps_id()[igen] == -15) )
                                    {
                                        lepPlus_status3_ = &(genps_p4().at(igen));
                                        lepPlus_status3_id_ = genps_id().at(igen);
                                        lepPlus_status3_orig = genps_p4().at(igen);
                                        lepPlus_status3_orig_ = &(lepPlus_status3_orig);

                                        //status = 1 lepton
                                        lepPlus_status3_nDaughters_ = genps_lepdaughter_id().at(igen).size();
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[igen].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[igen][kk];
                                            if ( daughterID == -11 || daughterID == -13 )
                                            {
                                                lepPlus_status1_ = &(genps_lepdaughter_p4()[igen][kk]);
                                                lepPlus_status1_id_ = daughterID;
                                                break;
                                            }
                                            //need to add all status=1 photons in a DR<0.1 cone around the lepton.
                                        }

                                    }
                                    else if ( (genps_id()[igen] == 12 || genps_id()[igen] == 14 ||  genps_id()[igen] == 16) )
                                    {
                                        nuPlus_status3_ = &(genps_p4().at(igen));
                                        nuPlus_status3_id_ = genps_id().at(igen);
                                        nuPlus_status3_orig = genps_p4().at(igen);
                                        nuPlus_status3_orig_ = &(nuPlus_status3_orig);

                                        //status = 1 neutrino
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[igen].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[igen][kk];
                                            if ( daughterID == genps_id()[igen] )
                                            {
                                                nuPlus_status1_ = &(genps_lepdaughter_p4()[igen][kk]);
                                                nuPlus_status1_id_ = daughterID;
                                                break;
                                            }
                                        }

                                    }

                                }
                                else if ( genps_id_mother()[igen] == -24 )
                                {
                                    if ( (genps_id()[igen] == 11 || genps_id()[igen] == 13 ||  genps_id()[igen] == 15) )
                                    {
                                        lepMinus_status3_ = &(genps_p4().at(igen));
                                        lepMinus_status3_id_ = genps_id().at(igen);
                                        lepMinus_status3_orig = genps_p4().at(igen);
                                        lepMinus_status3_orig_ = &(lepMinus_status3_orig);

                                        //status = 1 lepton
                                        lepMinus_status3_nDaughters_ = genps_lepdaughter_id().at(igen).size();
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[igen].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[igen][kk];
                                            if ( daughterID == 11 || daughterID == 13 )
                                            {
                                                lepMinus_status1_ = &(genps_lepdaughter_p4()[igen][kk]);
                                                lepMinus_status1_id_ = daughterID;
                                                break;
                                            }
                                            //need to add all status=1 photons in a DR<0.1 cone around the lepton.
                                        }

                                    }
                                    else if ( (genps_id()[igen] == -12 || genps_id()[igen] == -14 ||  genps_id()[igen] == -16) )
                                    {
                                        nuMinus_status3_ = &(genps_p4().at(igen));
                                        nuMinus_status3_id_ = genps_id().at(igen);
                                        nuMinus_status3_orig = genps_p4().at(igen);
                                        nuMinus_status3_orig_ = &(nuMinus_status3_orig);

                                        //status = 1 neutrino
                                        for (unsigned int kk = 0; kk < genps_lepdaughter_id()[igen].size(); kk++)
                                        {
                                            int daughterID = genps_lepdaughter_id()[igen][kk];
                                            if ( daughterID == genps_id()[igen] )
                                            {
                                                nuMinus_status1_ = &(genps_lepdaughter_p4()[igen][kk]);
                                                nuMinus_status1_id_ = daughterID;
                                                break;
                                            }
                                        }

                                    }

                                }






    //store stop
    if ( id == 1000006)
      stop_t_ = &(genps_p4().at(igen));   
    else if ( id == -1000006 )
      stop_tbar_ = &(genps_p4().at(igen));   

    //store neutralino
    if ( genps_id_mother().at(igen) == 1000006  && ( abs(id) == 1000022 ) ) {
      neutralino_t_ = &(genps_p4().at(igen));
    }
    if ( genps_id_mother().at(igen) == -1000006 && ( abs(id) == 1000022 ) ) {
      neutralino_tbar_ = &(genps_p4().at(igen));
    }
    
    //store c1/n2 for WH
    if ( abs(id) == 1000024)
      genc1_ = &(genps_p4().at(igen));   
    else if ( abs(id) == 1000023 )
      genn2_ = &(genps_p4().at(igen));   

    //store neutralinos for WH
    if ( abs(genps_id_mother().at(igen)) == 1000024  && ( abs(id) == 1000022 ) ) {
      neutralino_c1_ = &(genps_p4().at(igen));
    }
    if ( abs(genps_id_mother().at(igen)) == 1000023 && ( abs(id) == 1000022 ) ) {
      neutralino_n2_ = &(genps_p4().at(igen));
    }
    
    //store daughter lepton
    if ( abs(mothid) == 24 && (abs(id) == 11 || abs(id) == 13 || abs(id) ==15)) {

      if (genps_id_mother().at(igen)>0) {
        // lept 1 is the particle 
        lep_t_id_ = genps_id().at(igen);
        lep_t_ = &(genps_p4().at(igen));
      } else {
        // lept 2 is the anti-particle
        lep_tbar_id_ = genps_id().at(igen);
        lep_tbar_ = &(genps_p4().at(igen));
      }
    }

    // store lepton, neutrino and W for single lepton events                     
    if (pid==11 || pid==13) {
      foundlep = true;
      mclep_ = &genps_p4()[igen];
    }
    if (pid==12 || pid==14) {
      foundnu = true;
      mcnu_  = &genps_p4()[igen];
    }

    // store W or Z pT 
    // ignoring cases where have more than 1 boson for now
    if ( pid == 24 ) {
      ptwgen_ = genps_p4().at(igen).pt();
      foundwz = true;
      nwzpartons_  = 0;
    }
    if ( pid == 23 ) {
      ptzgen_ = genps_p4().at(igen).pt();
      foundwz = true;
      nwzpartons_  = 0;
    }

    if (foundwz && ( pid == 1 || pid == 2 || pid == 3 || pid == 4 || pid == 5 || pid == 6 || pid == 21 ) )   
      nwzpartons_++;

    // skip lines up to t and tbar
    if( igen < 8 ) continue;

    // require particle is a quark or a gluon
    if( !( pid==1 || pid==2 || pid==3 || pid==4 || pid==5 || pid==6 || pid == 21 ) ) continue;

    // save status 3 quarks/gluons
    //genqgs_.push_back(genps_p4().at(igen)); //for some reason this started crashing

    // require mother is not a top or W
    if( mothid == 6 || mothid == 24) continue;

    // found additional parton
    npartons_ ++;
    if( genps_p4().at(igen).pt() > maxpartonpt_ ) maxpartonpt_ = genps_p4().at(igen).pt();
    //    cout << "found parton, igen " << igen << " id " << pid << " motherid " << mothid << " pt " << genps_p4().at(igen).pt() << endl;
    
  }


  if (ismcatnlo && ( *topPlus_status3_!=topPlus_status2 || *topMinus_status3_!=topMinus_status2 ) ) cout<<" final top different from status=3 top, ntopPlusDaughters: "<<ntopPlusDaughters<<" ntopMinusDaughters: "<<ntopMinusDaughters<<endl; //no status=2 tops in pythia
  //if (ismcatnlo && ( ntops != 2 ) ) cout<<" ntops = "<<ntops<<" ntopPlusDaughters: "<<ntopPlusDaughters<<" ntopMinusDaughters: "<<ntopMinusDaughters<<endl; //no status=2 tops in pythia


  if(!foundbPlus || !foundbMinus) {
    cout<<"One of the bs is missing! "<< nbs_ <<endl;
    dumpDocLines();
  }

  // check explicitly for t and tbar, in case model has multiple tops etc
  if (topPlus_status3_ && topMinus_status3_ && ntops==2) {
    ttpair = *topPlus_status3_ + *topMinus_status3_;
    ttbar_    = &ttpair;
    ptttbar_  = ttbar_->pt();
    mttbar_   = ttbar_->mass();
    etattbar_ = ttbar_->eta();
    rapidityttbar_ = ttbar_->Rapidity();
  }

//  if (foundlep && foundnu) {
//    mcmln_ = (*mclep_+*mcnu_).mass();
//    mcmtln_ = getMT( mclep_->Pt() , mclep_->Phi() , mcnu_->Pt() , mcnu_->Phi() );
//  }

  //remainder of function is only for ttdl
  if( !(nleps==2 && ntops==2) ) return; //in MC@NLO there are a few events with 4 tops that cause a crash

  //For MC@NLO. Boost status=3 W back to status=2
  if(ismcatnlo) {
    //if(ntaus==0 && (*lepPlus_status1_!=*lepPlus_status3_ || *lepMinus_status1_!=*lepMinus_status3_) ) cout<<" status 1 and 3 leptons not identical "<<lepPlus_status3_->E()-lepPlus_status1_->E()<<" "<<lepMinus_status3_->E()-lepMinus_status1_->E()<<endl;

    topPlus_status1_ = &(vt_stat1);
    topMinus_status1_ = &(vtbar_stat1);
    bPlus_status1_ = &(vb_stat1);
    bMinus_status1_ = &(vbbar_stat1);

    //if(ntaus==0 && ntopPlusDaughters==2 && fabs( (*bPlus_status1_+*lepPlus_status1_+*nuPlus_status1_).E() - topPlus_status1_->E() )>1e-3 ) cout<<"Ndaughters_topPlus: "<<ntopPlusDaughters<<" "<< (*bPlus_status1_+*lepPlus_status1_+*nuPlus_status1_).E() - topPlus_status1_->E() <<endl;
    //if(ntaus==0 && ntopMinusDaughters==2 && fabs( (*bMinus_status1_+*lepMinus_status1_+*nuMinus_status1_).E() - topMinus_status1_->E() )>1e-3 ) cout<<"Ndaughters_topMinus: "<<ntopMinusDaughters<<" "<< (*bMinus_status1_+*lepMinus_status1_+*nuMinus_status1_).E() - topMinus_status1_->E() <<endl;

    //status=1 top should not include the ISR jet, so sum the b+l+nu (this method doesn't work when there are taus, but then we can't easily define a status=1 top in any case)
    //if(ntaus==0) topPlus_status1_ = bPlus_status1_+lepPlus_status1_+nuPlus_status1_;
    if(ntaus==0) topPlus_status1_->SetXYZT(vb_stat1.x()+lepPlus_status1_->x()+nuPlus_status1_->x(),
                    vb_stat1.y()+lepPlus_status1_->y()+nuPlus_status1_->y(),
                    vb_stat1.z()+lepPlus_status1_->z()+nuPlus_status1_->z(),
                    vb_stat1.t()+lepPlus_status1_->t()+nuPlus_status1_->t()
                   );
    //if(ntaus==0) topMinus_status1_ = bMinus_status1_+lepMinus_status1_+nuMinus_status1_;
    if(ntaus==0) topMinus_status1_->SetXYZT(vbbar_stat1.x()+lepMinus_status1_->x()+nuMinus_status1_->x(),
                    vbbar_stat1.y()+lepMinus_status1_->y()+nuMinus_status1_->y(),
                    vbbar_stat1.z()+lepMinus_status1_->z()+nuMinus_status1_->z(),
                    vbbar_stat1.t()+lepMinus_status1_->t()+nuMinus_status1_->t()
                   );


    //WPlus
    TLorentzVector WPlus_status3B;
    WPlus_status3B.SetXYZT(WPlus_status3_->x(),
                    WPlus_status3_->y(),
                    WPlus_status3_->z(),
                    WPlus_status3_->t()
                   );
    TLorentzVector lepPlus_status3B;
    lepPlus_status3B.SetXYZT(lepPlus_status3_->x(),
                    lepPlus_status3_->y(),
                    lepPlus_status3_->z(),
                    lepPlus_status3_->t()
                   );
    TLorentzVector nuPlus_status3B;
    nuPlus_status3B.SetXYZT(nuPlus_status3_->x(),
                    nuPlus_status3_->y(),
                    nuPlus_status3_->z(),
                    nuPlus_status3_->t()
                   );


    WPlus_status3B.Boost(-WPlus_status2_T.BoostVector());

    lepPlus_status3B.Boost(-WPlus_status2_T.BoostVector());
    lepPlus_status3B.Boost(-WPlus_status3B.BoostVector());
    lepPlus_status3B.Boost(WPlus_status2_T.BoostVector());
    nuPlus_status3B.Boost(-WPlus_status2_T.BoostVector());
    nuPlus_status3B.Boost(-WPlus_status3B.BoostVector());
    nuPlus_status3B.Boost(WPlus_status2_T.BoostVector());

    WPlus_status3B.Boost(-WPlus_status3B.BoostVector());
    WPlus_status3B.Boost(WPlus_status2_T.BoostVector());


    //WMinus
    TLorentzVector WMinus_status3B;
    WMinus_status3B.SetXYZT(WMinus_status3_->x(),
                    WMinus_status3_->y(),
                    WMinus_status3_->z(),
                    WMinus_status3_->t()
                   );
    TLorentzVector lepMinus_status3B;
    lepMinus_status3B.SetXYZT(lepMinus_status3_->x(),
                    lepMinus_status3_->y(),
                    lepMinus_status3_->z(),
                    lepMinus_status3_->t()
                   );
    TLorentzVector nuMinus_status3B;
    nuMinus_status3B.SetXYZT(nuMinus_status3_->x(),
                    nuMinus_status3_->y(),
                    nuMinus_status3_->z(),
                    nuMinus_status3_->t()
                   );

    WMinus_status3B.Boost(-WMinus_status2_T.BoostVector());

    lepMinus_status3B.Boost(-WMinus_status2_T.BoostVector());
    lepMinus_status3B.Boost(-WMinus_status3B.BoostVector());
    lepMinus_status3B.Boost(WMinus_status2_T.BoostVector());
    nuMinus_status3B.Boost(-WMinus_status2_T.BoostVector());
    nuMinus_status3B.Boost(-WMinus_status3B.BoostVector());
    nuMinus_status3B.Boost(WMinus_status2_T.BoostVector());

    WMinus_status3B.Boost(-WMinus_status3B.BoostVector());
    WMinus_status3B.Boost(WMinus_status2_T.BoostVector());


    //cout<<WPlus_status3B.E() - WPlus_status2_T.E()<<" "<<WMinus_status3B.E() - WMinus_status2_T.E()<<endl;
    //cout<<(lepPlus_status3B+nuPlus_status3B).E() - WPlus_status2_T.E()<<" "<<(lepMinus_status3B+nuMinus_status3B).E() - WMinus_status2_T.E()<<endl<<endl;

    //if(topPlus_status3_->E()-(WPlus_status2_T + bPlus_status3_).E() > 1e-4) cout<<" top has >2 daughters "<<topPlus_status3_->E()-(WPlus_status2_T + bPlus_status3_).E()<<endl;

    lepPlus_status3_->SetXYZT(lepPlus_status3B.X(),
                    lepPlus_status3B.Y(),
                    lepPlus_status3B.Z(),
                    lepPlus_status3B.T()
                   );
    nuPlus_status3_->SetXYZT(nuPlus_status3B.X(),
                    nuPlus_status3B.Y(),
                    nuPlus_status3B.Z(),
                    nuPlus_status3B.T()
                   );
    WPlus_status3_->SetXYZT(WPlus_status2.X(),
                    WPlus_status2.Y(),
                    WPlus_status2.Z(),
                    WPlus_status2.T()
                   );


    lepMinus_status3_->SetXYZT(lepMinus_status3B.X(),
                    lepMinus_status3B.Y(),
                    lepMinus_status3B.Z(),
                    lepMinus_status3B.T()
                   );
    nuMinus_status3_->SetXYZT(nuMinus_status3B.X(),
                    nuMinus_status3B.Y(),
                    nuMinus_status3B.Z(),
                    nuMinus_status3B.T()
                   );
    WMinus_status3_->SetXYZT(WMinus_status2.X(),
                    WMinus_status2.Y(),
                    WMinus_status2.Z(),
                    WMinus_status2.T()
                   );


    //also recompute status=3 tops due to presence of events in MC@NLO with gluon FSR in the top decay. Note this means the effective top has lower mass. This is probably not what we want, because the hard FSR gluon is probably from the b (according to Mrenna).
    //if(ntopPlusDaughters>2) cout<< " Ndaughters_topPlus: "<<ntopPlusDaughters<<" "<<topPlus_status3_->M()<<" "<<(WPlus_status2_T + bPlus_status3_).M()<<endl;
    //topPlus_status3_ = WPlus_status2 + bPlus_status3_; 
    topPlus_status3_bW_->SetXYZT(WPlus_status2.X()+bPlus_status3_->X(),
                    WPlus_status2.Y()+bPlus_status3_->Y(),
                    WPlus_status2.Z()+bPlus_status3_->Z(),
                    WPlus_status2.T()+bPlus_status3_->T()
                   );
    //topMinus_status3_ = WMinus_status2 + bMinus_status3_;
    topMinus_status3_bW_->SetXYZT(WMinus_status2.X()+bMinus_status3_->X(),
                    WMinus_status2.Y()+bMinus_status3_->Y(),
                    WMinus_status2.Z()+bMinus_status3_->Z(),
                    WMinus_status2.T()+bMinus_status3_->T()
                   );

    //check for error in calculation
    if (  fabs( (*lepPlus_status3_+*nuPlus_status3_).M() - WPlus_status2_T.M() ) > 1e-4  ) {   //the stat2 and 3 Ws can have slightly different masses - in that case just check their betas (Ws and tops will no longer exactly match)
      if (  fabs( (*lepPlus_status3_+*nuPlus_status3_).Beta() - WPlus_status2_T.Beta() ) > 1e-4  ) cout<<"********* Inconsistent WPlus ***********"<<endl;
    }
    else if( fabs( (bPlus_status3_->E()+lepPlus_status3_->E()+nuPlus_status3_->E()) - topPlus_status3_->E() ) > 1e-2 && ntopPlusDaughters == 2 ) {
      cout<<" Top daughters don't match top. Ndaughters_topPlus: "<<ntopPlusDaughters<<" lepPlusID: "<<lepPlus_status3_id_<<" evt: "<<evt_event()<<endl;
      cout<<(bPlus_status3_->E()+lepPlus_status3_->E()+nuPlus_status3_->E()) - topPlus_status3_->E()<<" W: "<<(lepPlus_status3_->E()+nuPlus_status3_->E()) - WPlus_status2_T.E()<<" mW: "<<(*lepPlus_status3_ + *nuPlus_status3_).M() - WPlus_status2_T.M()<<" gammaW: "<< WPlus_status2_T.Gamma() <<endl; //here we expect exact agreement
      cout<<(bPlus_status3_->E()+lepPlus_status3_->E()+nuPlus_status3_->E()) - topPlus_status2.E()<<endl; //here we expect a difference when Ndaughters!=2
      cout<<" WPlus2 "<<WPlus_status2_T.Px()<<" "<<WPlus_status2_T.Py()<<" "<<WPlus_status2_T.Pz()<<" "<<WPlus_status2_T.M()<<endl;
      cout<<" WPlus3 "<<WPlus_status3B.Px()<<" "<<WPlus_status3B.Py()<<" "<<WPlus_status3B.Pz()<<" "<<WPlus_status3B.M()<<endl;
      cout<<" lep+nuPlus "<<(nuPlus_status3B+lepPlus_status3B).Px()<<" "<<(nuPlus_status3B+lepPlus_status3B).Py()<<" "<<(nuPlus_status3B+lepPlus_status3B).Pz()<<" "<<(nuPlus_status3B+lepPlus_status3B).M()<<endl;
    }
    if (  fabs( (*lepMinus_status3_+*nuMinus_status3_).M() - WMinus_status2_T.M() ) > 1e-4  ) {   //the stat2 and 3 Ws can have slightly different masses - in that case just check their betas (Ws and tops will no longer exactly match)
      if (  fabs( (*lepMinus_status3_+*nuMinus_status3_).Beta() - WMinus_status2_T.Beta() ) > 1e-4  ) cout<<"********* Inconsistent WMinus ***********"<<endl;
    }
    else if( fabs( (bMinus_status3_->E()+lepMinus_status3_->E()+nuMinus_status3_->E()) - topMinus_status3_->E() ) > 1e-2 && ntopMinusDaughters == 2 ) {
      cout<<" Top daughters don't match top. Ndaughters_topMinus: "<<ntopMinusDaughters<<" lepMinusID: "<<lepMinus_status3_id_<<" evt: "<<evt_event()<<endl;
      cout<<(bMinus_status3_->E()+lepMinus_status3_->E()+nuMinus_status3_->E()) - topMinus_status3_->E()<<" W: "<<(lepMinus_status3_->E()+nuMinus_status3_->E()) - WMinus_status2_T.E()<<" mW: "<<(*lepMinus_status3_ + *nuMinus_status3_).M() - WMinus_status2_T.M()<<" gammaW: "<< WMinus_status2_T.Gamma() <<endl; //here we expect exact agreement
      cout<<(bMinus_status3_->E()+lepMinus_status3_->E()+nuMinus_status3_->E()) - topMinus_status2.E()<<endl; //here we expect a difference when Ndaughters!=2
      cout<<" WMinus2 "<<WMinus_status2_T.Px()<<" "<<WMinus_status2_T.Py()<<" "<<WMinus_status2_T.Pz()<<" "<<WMinus_status2_T.M()<<endl;
      cout<<" WMinus3 "<<WMinus_status3B.Px()<<" "<<WMinus_status3B.Py()<<" "<<WMinus_status3B.Pz()<<" "<<WMinus_status3B.M()<<endl;
      cout<<" lep+nuMinus "<<(nuMinus_status3B+lepMinus_status3B).Px()<<" "<<(nuMinus_status3B+lepMinus_status3B).Py()<<" "<<(nuMinus_status3B+lepMinus_status3B).Pz()<<" "<<(nuMinus_status3B+lepMinus_status3B).M()<<endl;
    }

    //if(ntaus==0) cout<< " Ndaughters_topPlus: "<<ntopPlusDaughters<<" Ndaughters_topMinus: "<<ntopMinusDaughters<<" "<< topPlus_status3_->E() - topPlus_status1_->E() << " " << topMinus_status3_->E() - topMinus_status1_->E() <<endl;

  }


  //calculate gen-level quantities
  TLorentzVector topplus_genp_p4(0, 0, 0, 0), topminus_genp_p4(0, 0, 0, 0), lepPlus_gen(0, 0, 0, 0), lepMinus_gen(0, 0, 0, 0);
  TLorentzVector cms_gen(0, 0, 0, 0);

  topplus_genp_p4.SetXYZT(topPlus_status3_->x(),
                  topPlus_status3_->y(),
                  topPlus_status3_->z(),
                  topPlus_status3_->t()
                 );

  topminus_genp_p4.SetXYZT(topMinus_status3_->x(),
                  topMinus_status3_->y(),
                  topMinus_status3_->z(),
                  topMinus_status3_->t()
                 );

  lepPlus_gen.SetXYZT(lepPlus_status3_->x(),
                  lepPlus_status3_->y(),
                  lepPlus_status3_->z(),
                  lepPlus_status3_->t()
                 );

  lepMinus_gen.SetXYZT(lepMinus_status3_->x(),
                  lepMinus_status3_->y(),
                  lepMinus_status3_->z(),
                  lepMinus_status3_->t()
                 );


  m_topminus_gen_ = topminus_genp_p4.M();
  m_topplus_gen_ = topplus_genp_p4.M();

  tt_mass_gen_ = (topplus_genp_p4 + topminus_genp_p4).M();
  ttRapidity_gen_ = topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity();
  ttRapidity2_gen_ = (topplus_genp_p4 + topminus_genp_p4).Rapidity();

  top_rapiditydiff_cms_gen_ = (topplus_genp_p4.Rapidity() - topminus_genp_p4.Rapidity()) * (topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity());
  top_pseudorapiditydiff_cms_gen_ = abs(topplus_genp_p4.Eta()) - abs(topminus_genp_p4.Eta());
  top_rapiditydiff_Marco_gen_ = abs(topplus_genp_p4.Rapidity()) - abs(topminus_genp_p4.Rapidity());


  cms_gen = topplus_genp_p4 + topminus_genp_p4;
  tt_pT_gen_ = cms_gen.Pt();
  topplus_genp_p4.Boost(-cms_gen.BoostVector());
  topminus_genp_p4.Boost(-cms_gen.BoostVector());
  top_costheta_cms_gen_ = topplus_genp_p4.Vect().Dot(cms_gen.Vect()) / (topplus_genp_p4.Vect().Mag() * cms_gen.Vect().Mag());


  lep_charge_asymmetry_gen_ = abs(lepPlus_gen.Eta()) - abs(lepMinus_gen.Eta());
  lep_azimuthal_asymmetry_gen_ = lepPlus_gen.DeltaPhi(lepMinus_gen);
  lep_azimuthal_asymmetry2_gen_ = acos(cos(lep_azimuthal_asymmetry_gen_));

  lepPlus_gen.Boost(-cms_gen.BoostVector());
  lepPlus_gen.Boost(-topplus_genp_p4.BoostVector());
  lepMinus_gen.Boost(-cms_gen.BoostVector());
  lepMinus_gen.Boost(-topminus_genp_p4.BoostVector());

  lepPlus_costheta_cms_gen_ = lepPlus_gen.Vect().Dot(topplus_genp_p4.Vect()) / (lepPlus_gen.Vect().Mag() * topplus_genp_p4.Vect().Mag());
  lepMinus_costheta_cms_gen_ = lepMinus_gen.Vect().Dot(topminus_genp_p4.Vect()) / (lepMinus_gen.Vect().Mag() * topminus_genp_p4.Vect().Mag());

  top_spin_correlation_gen_ = lepPlus_costheta_cms_gen_ * lepMinus_costheta_cms_gen_ ;
  lep_cos_opening_angle_gen_ = lepPlus_gen.Vect().Dot(lepMinus_gen.Vect()) / (lepPlus_gen.Vect().Mag() * lepMinus_gen.Vect().Mag());



  //recalculate using original leptons
  topplus_genp_p4.SetXYZT(topPlus_status3_->x(),
                  topPlus_status3_->y(),
                  topPlus_status3_->z(),
                  topPlus_status3_->t()
                 );

  topminus_genp_p4.SetXYZT(topMinus_status3_->x(),
                  topMinus_status3_->y(),
                  topMinus_status3_->z(),
                  topMinus_status3_->t()
                 );

  lepPlus_gen.SetXYZT(lepPlus_status3_orig_->x(),
                  lepPlus_status3_orig_->y(),
                  lepPlus_status3_orig_->z(),
                  lepPlus_status3_orig_->t()
                 );

  lepMinus_gen.SetXYZT(lepMinus_status3_orig_->x(),
                  lepMinus_status3_orig_->y(),
                  lepMinus_status3_orig_->z(),
                  lepMinus_status3_orig_->t()
                 );


  m_topminus_gen_origleps_ = topminus_genp_p4.M();
  m_topplus_gen_origleps_ = topplus_genp_p4.M();

  tt_mass_gen_origleps_ = (topplus_genp_p4 + topminus_genp_p4).M();
  ttRapidity_gen_origleps_ = topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity();
  ttRapidity2_gen_origleps_ = (topplus_genp_p4 + topminus_genp_p4).Rapidity();

  top_rapiditydiff_cms_gen_origleps_ = (topplus_genp_p4.Rapidity() - topminus_genp_p4.Rapidity()) * (topplus_genp_p4.Rapidity() + topminus_genp_p4.Rapidity());
  top_pseudorapiditydiff_cms_gen_origleps_ = abs(topplus_genp_p4.Eta()) - abs(topminus_genp_p4.Eta());
  top_rapiditydiff_Marco_gen_origleps_ = abs(topplus_genp_p4.Rapidity()) - abs(topminus_genp_p4.Rapidity());


  cms_gen = topplus_genp_p4 + topminus_genp_p4;
  tt_pT_gen_origleps_ = cms_gen.Pt();
  topplus_genp_p4.Boost(-cms_gen.BoostVector());
  topminus_genp_p4.Boost(-cms_gen.BoostVector());
  top_costheta_cms_gen_origleps_ = topplus_genp_p4.Vect().Dot(cms_gen.Vect()) / (topplus_genp_p4.Vect().Mag() * cms_gen.Vect().Mag());


  lep_charge_asymmetry_gen_origleps_ = abs(lepPlus_gen.Eta()) - abs(lepMinus_gen.Eta());
  lep_azimuthal_asymmetry_gen_origleps_ = lepPlus_gen.DeltaPhi(lepMinus_gen);
  lep_azimuthal_asymmetry2_gen_origleps_ = acos(cos(lep_azimuthal_asymmetry_gen_origleps_));

  lepPlus_gen.Boost(-cms_gen.BoostVector());
  lepPlus_gen.Boost(-topplus_genp_p4.BoostVector());
  lepMinus_gen.Boost(-cms_gen.BoostVector());
  lepMinus_gen.Boost(-topminus_genp_p4.BoostVector());

  lepPlus_costheta_cms_gen_origleps_ = lepPlus_gen.Vect().Dot(topplus_genp_p4.Vect()) / (lepPlus_gen.Vect().Mag() * topplus_genp_p4.Vect().Mag());
  lepMinus_costheta_cms_gen_origleps_ = lepMinus_gen.Vect().Dot(topminus_genp_p4.Vect()) / (lepMinus_gen.Vect().Mag() * topminus_genp_p4.Vect().Mag());

  top_spin_correlation_gen_origleps_ = lepPlus_costheta_cms_gen_origleps_ * lepMinus_costheta_cms_gen_origleps_ ;
  lep_cos_opening_angle_gen_origleps_ = lepPlus_gen.Vect().Dot(lepMinus_gen.Vect()) / (lepPlus_gen.Vect().Mag() * lepMinus_gen.Vect().Mag());


}

bool topAFB_looper::isBHadronPdgId(int PdgId)
{
  //borrowed from https://github.com/cms-kr/KrAFT/blob/master/GeneratorTools/plugins/GhostGenParticleProducer.cc#L166
  int absPdgId = abs(PdgId);
  if ( absPdgId <= 100 ) return false; // Fundamental particles and MC internals
  if ( absPdgId >= 1000000000 ) return false; // Nuclei, +-10LZZZAAAI

  // General form of PDG ID is 7 digit form
  // +- n nr nL nq1 nq2 nq3 nJ
  //const int nJ = absPdgId % 10; // Spin
  const int nq3 = (absPdgId / 10) % 10;
  const int nq2 = (absPdgId / 100) % 10;
  const int nq1 = (absPdgId / 1000) % 10;

  if ( nq3 == 0 ) return false; // Diquarks
  if ( nq1 == 0 and nq2 == 5 ) return true; // B mesons
  if ( nq1 == 5 ) return true; // B baryons

  return false;
}


bool topAFB_looper::isBHadron(int igen)
{
  if( genps_p4().at(igen).Pt()<5 ) return false; //Pt>5 required: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/ParticleLevelTopDefinitions
  if ( !isBHadronPdgId( genps_id().at(igen) ) ) return false;
  for (unsigned int kk = 0; kk < genps_lepdaughter_id().at(igen).size(); ++kk)
  {
    //cout<<genps_lepdaughter_id()[igen][kk]<<endl;
    if(isBHadronPdgId(genps_lepdaughter_id()[igen][kk])) return false; //this doesn't work with the current ntuples because the daughters are the final status=1 particles rather than the immediate status=2 daughters: need to re-ntuplise
  }

  return true;
}




