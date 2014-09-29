#ifndef STOPTREELOOPER_H
#define STOPTREELOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"

#include <iostream>
#include "Math/LorentzVector.h"

#include <cmath>
#include <map>

using namespace std;

class StopTreeLooper
{

public:
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

    StopTreeLooper();
    ~StopTreeLooper();

    void setOutFileName(string filename);
    void loop(TChain *chain, TString name);

    //plotting
    void makeSIGPlots(float evtweight, std::map<std::string, TH1D *> &h_1d,
                      string tag_selection, string flav_tag );
    void makettPlots(float evtweight, std::map<std::string, TH1D *> &h_1d, std::map<std::string, TH2D *> &h_2d,
                     string tag_selection, string flav_tag );
    void makeAccPlots(float evtweight, std::map<std::string, TH1D *> &h_1d, std::map<std::string, TH2D *> &h_2d,
                      string tag_selection, string flav_tag );
    void makeCR1Plots(float evtweight, std::map<std::string, TH1D *> &h_1d,
                      string tag_selection, string flav_tag );
    void makeCR2Plots(float evtweight, std::map<std::string, TH1D *> &h_1d,
                      string tag_selection, string flav_tag_dl );
    void makeCR3Plots(float evtweight, std::map<std::string, TH1D *> &h_1d,
                      string tag_selection, string flav_tag_dl );
    void makeNJPlots( float evtweight, std::map<std::string, TH1D *> &h_1d,
                      string tag_selection, string flav_tag );
    void makeZPlots(  float evtweight, std::map<std::string, TH1D *> &h_1d,
                      string tag_selection, string flav_tag );

    //ttbar solver
    void solvettbar();
    double get_pdf_weight( TLorentzVector &t1, TLorentzVector &t2 );
    double get_dalitz_prob( TLorentzVector &lep, TLorentzVector &top );

    //baby ntuples
    void MakeBabyNtuple(const char *babyFilename);
    void FillBabyNtuple();
    void CloseBabyNtuple();

    //selection
    bool passFullSelection(bool isData);

    //top pT reweighting
    double TopPtWeight(double topPt);


private:

    TFile *babyFile_;
    TTree *babyTree_;

    string m_outfilename_;
    // njets requirement
    int min_njets;
    //for phi corrected met
    float t1metphicorr;
    float t1metphicorrphi;
    float pfcalo_metratio;
    float pfcalo_metdphi;
    float pfcalo_deltamet;

    float mlb_1;
    float mlb_2;
    float mlb_3;
    float mlb_4;
    float mlb_min;

    //jets information
    int n_jets;
    int n_bjets;
    int n_ljets;
    vector<LorentzVector> jets;
    vector<LorentzVector> bjets;
    vector<LorentzVector> nonbjets;
    vector<LorentzVector> bcandidates;
    vector<float> btag;
    vector<float> sigma_jets;
    vector<int> mc;
    bool tobtecveto_;

    //variables for baby ntuples
    Int_t   run;
    Int_t   ls;
    Int_t   evt;
    Int_t   channel;
    Int_t   genchannel;
    double weight;

    float lep_charge_asymmetry;
    float lep_azimuthal_asymmetry;
    float lep_azimuthal_asymmetry2;
    float top_rapiditydiff_cms;
    float top_pseudorapiditydiff_cms;
    float top_rapiditydiff_Marco;
    float top_costheta_cms;
    float lepPlus_costheta_cms;
    float lepMinus_costheta_cms;
    float top_spin_correlation;
    float lep_cos_opening_angle;
    float dilmass;
    float tt_mass;
    float ttRapidity2;
    float tt_pT;
    float top1_pt;
    float top2_pt;
    float top1_p_CM;
    float top2_p_CM;
    float top_rapiditydiffsigned_cms;

    float lep_charge_asymmetry_gen;
    float lep_azimuthal_asymmetry_gen;
    float lep_azimuthal_asymmetry2_gen;
    float top_rapiditydiff_cms_gen;
    float top_pseudorapiditydiff_cms_gen;
    float top_rapiditydiff_Marco_gen;
    float top_costheta_cms_gen;
    float lepPlus_costheta_cms_gen;
    float lepMinus_costheta_cms_gen;
    float top_spin_correlation_gen;
    float lep_cos_opening_angle_gen;
    float tt_mass_gen;
    float ttRapidity2_gen;
    float tt_pT_gen;
    float top1_pt_gen;
    float top2_pt_gen;


    //variables used by solver
    TLorentzVector lepPlus;
    TLorentzVector lepMinus;
    TLorentzVector jet1;
    TLorentzVector jet2;
    TLorentzVector top1_p4;
    TLorentzVector top2_p4;
    TLorentzVector nusum;
    TLorentzVector cms;

    TLorentzVector lepPlus_gen;
    TLorentzVector lepMinus_gen;
    TLorentzVector topplus_genp_p4;
    TLorentzVector topminus_genp_p4;
    TLorentzVector cms_gen;
    TLorentzVector nuPlus_gen;
    TLorentzVector nuMinus_gen;

    vector <TLorentzVector> nu1_vecs;
    vector <TLorentzVector> nu2_vecs;
    vector <TLorentzVector> top1_vecs;
    vector <TLorentzVector> top2_vecs;
    vector <double> AMWT_weights;
    float m_top;
    float m_top_B;
    float closestDeltaMET_maxwcombo;
    float closestDeltaMET_othercombo;
    float closestDeltaMET_bestcombo;

    int imaxweight;
    bool closestApproach;

};

#endif
