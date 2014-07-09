#ifndef topAFB_looper_h
#define topAFB_looper_h

#include <stdint.h>

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;
typedef vector<ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > > VofP4;

class topAFB_looper
{
public:
    topAFB_looper();
    ~topAFB_looper();
    void ScanChain (TChain *, vector<TString> , string, float);
    void bookHistos(const char *sample, int nchannels, int nhistsets);
    double TopPtWeight(double topPt);
    void fillUnderOverFlow(TH1D *h1, float value, double weight = 1., int Nsolns = 1);
    void fillUnderOverFlow(TH2D *h2, float xvalue, float yvalue, double weight = 1., int Nsolns = 1);
    //void fillUnderOverFlow(TProfile *h2, float xvalue, float yvalue);
    void fillOverFlow(TH1D *h1, float value, float weight = 1.);
    void fillOverFlow(TH2D *h2, float xvalue, float yvalue, float weight = 1.);
    void fillHistos(TH1D *h1[4][4], float value, double weight, int myType, int nJetsIdx, int Nsolns = 1);
    void fillHistos(TH2D *h2[4][4], float xvalue, float yvalue, double weight, int myType, int nJetsIdx, int Nsolns = 1);
    void fillHistos(TProfile *h2[4][4], float xvalue, float yvalue,  int myType, int nJetsIdx);
    int leptonGenpCount_lepTauDecays_status3only(int &nele, int &nmuon, int &ntau);
    void fillgenlevel(bool ismcatnlo, int nleps, int ntaus);
    bool isBHadronPdgId(int PdgId);
    bool isBHadron(int igen);

private:

    bool applyTopPtWeighting;
    bool weighttaudecay;

    // event identification

    Int_t   run_;
    Int_t   ls_;
    Int_t   evt_;
    Int_t ndavtx_;
    float   t_mass_;
    double weight_;
    Int_t Nsolns_;


    float massltb_;
    float massllb_;
    float dr_ltjet_gen_;
    float dr_lljet_gen_;
    float tt_mass_ ;
    float ttRapidity_ ;
    float ttRapidity2_ ;
    float ttPt_ ;
    float lep_charge_asymmetry_ ;
    float lep_pseudorap_diff_ ;
    float lep_azimuthal_asymmetry_ ;
    float lep_azimuthal_asymmetry2_ ;
    float top_spin_correlation_ ;
    float lep_cos_opening_angle_ ;
    float top_costheta_cms_    ;
    float top_rapiditydiff_cms_    ;
    float top_rapiditydiff_Marco_    ;
    float top_pseudorapiditydiff_cms_    ;
    float lepPlus_costheta_cms_ ;
    float lepMinus_costheta_cms_ ;
    /*
    float tt_mass_gen_ ;
    float ttRapidity_gen_ ;
    float ttRapidity2_gen_ ;
    float ttPt_gen_ ;
    float lep_charge_asymmetry_gen_ ;
    float lep_azimuthal_asymmetry_gen_ ;
    float lep_azimuthal_asymmetry2_gen_ ;
    float top_spin_correlation_gen_ ;
    float lep_cos_opening_angle_gen_ ;
    float top_costheta_cms_gen_    ;
    float lepPlus_costheta_cms_gen_ ;
    float lepMinus_costheta_cms_gen_ ;
    float top_rapiditydiff_cms_gen_    ;
    float top_rapiditydiff_Marco_gen_    ;
    float top_pseudorapiditydiff_cms_gen_;
    */

    TH1D *hnJet[4][4];                   // Njet distributions
    TH1D *hnBtagJet[4][4];                   // NBTagjet distributions
    TH1D *hnVtx[4][4];
    TH1D *hAMWTweightnojetsmear[4][4];
    //Top Mass Plots
    TH1D *httRapidity[4][4];
    TH1D *httRapidity2[4][4];
    TH1D *httpT[4][4];
    TH1D *httMass[4][4];
    TH1D *httMass_pull[4][4];
    TH1D *httMass_diff[4][4];
    TH1D *htopMass_diff_plus[4][4];
    TH1D *htopMass_diff_minus[4][4];
    TH1D *htopPCM_diff_plus[4][4];
    TH1D *htopPCM_diff_minus[4][4];

    TH1D *htopMass[4][4];
    TH1D *httpT_gen[4][4];
    TH1D *httMass_gen[4][4];
    TH1D *htopMass_plus_gen[4][4];
    TH1D *htopMass_minus_gen[4][4];
    TH1D *hmassllb[4][4];
    TH1D *hmassltb[4][4];
    TH1D *hmassllb1Dmasscut[4][4];
    TH1D *hmassltb1Dmasscut[4][4];
    TH1D *htheSumJetPt[4][4];
    TH1D *htheSumBtagJetPt[4][4];
    TH1D *hthefirstJetPt[4][4];
    TH1D *hthesecondJetPt[4][4];
    TH1D *htheleadinglepPt[4][4];
    TH1D *hthesecondlepPt[4][4];
    TH1D *hthesumlepPt[4][4];
    TH1D *hlepEta[4][4];
    TH1D *hjetPt[4][4];
    TH1D *hjetEta[4][4];
    TH1D *hMET[4][4];
    TH1D *hdMET[4][4];
    TH1D *htheSumLepPt[4][4];
    TH1D *htopCosTheta[4][4];
    TH1D *hpseudorapiditydiff[4][4];
    TH1D *hrapiditydiff[4][4];
    TH1D *hlepCosTheta[4][4];
    TH1D *hlepChargeAsym[4][4];
    TH1D *hlepAzimAsym[4][4];
    TH1D *hlepAzimAsym2[4][4];
    TH1D *htopSpinCorr[4][4];
    TH1D *hlepCosOpeningAngle[4][4];
    TH1D *htopCosTheta_gen[4][4];
    TH1D *hlepCosTheta_gen[4][4];
    TH1D *hlepChargeAsym_gen[4][4];
    TH1D *hlepAzimAsym_gen[4][4];
    TH1D *hlepAzimAsym2_gen[4][4];
    TH1D *htopSpinCorr_gen[4][4];
    TH1D *hlepCosOpeningAngle_gen[4][4];

    TH1D *hlepRapDiff[4][4];
    TH1D *hlepAngleBetween[4][4];
    TH1D *hlepAngleBetweenCMS[4][4];
    TH1D *hpseudorapiditydiff2[4][4];
    TH1D *hrapiditydiff2[4][4];
    TH1D *hlepPlusCosTheta[4][4];
    TH1D *hlepMinusCosTheta[4][4];
    TH1D *hjetAzimAsym[4][4];
    TH1D *hjetRapDiff[4][4];
    TH1D *hjetAngleBetween[4][4];
    TH1D *hjetAngleBetweenCMS[4][4];
    TH1D *hlepPhi[4][4];
    TH1D *hlepPlusPhi[4][4];
    TH1D *hlepMinusPhi[4][4];
    TH1D *hjetPhi[4][4];
    TH1D *hlepPlusEta[4][4];
    TH1D *hlepMinusEta[4][4];
    TH1D *hlepPlusPt[4][4];
    TH1D *hlepMinusPt[4][4];

    TH1D *hlepPt[4][4];
    TH1D *hrapiditydiffMarco[4][4];

    TH1D *hlepPlusCosTheta_gen[4][4];
    TH1D *hlepMinusCosTheta_gen[4][4];
    TH1D *hpseudorapiditydiff_gen[4][4];
    TH1D *hrapiditydiff_gen[4][4];
    TH1D *hrapiditydiffMarco_gen[4][4];

    TH2D *hlepChargeAsym_gen2d[4][4];
    TH2D *hlepAzimAsym_gen2d[4][4];
    TH2D *hlepAzimAsym2_gen2d[4][4];
    TH2D *htopSpinCorr_gen2d[4][4];
    TH2D *hlepCosOpeningAngle_gen2d[4][4];
    TH2D *htopCosTheta_gen2d[4][4];
    TH2D *hlepCosTheta_gen2d[4][4];
    TH2D *hlepPlusCosTheta_gen2d[4][4];
    TH2D *hlepMinusCosTheta_gen2d[4][4];
    TH2D *hpseudorapiditydiff_gen2d[4][4];
    TH2D *hrapiditydiff_gen2d[4][4];
    TH2D *hrapiditydiffMarco_gen2d[4][4];

    TH2D *hlepChargeAsym_ttpT_gen2d[4][4];
    TH2D *hlepAzimAsym_ttpT_gen2d[4][4];
    TH2D *hlepAzimAsym2_ttpT_gen2d[4][4];
    TH2D *htopSpinCorr_ttpT_gen2d[4][4];
    TH2D *hlepCosOpeningAngle_ttpT_gen2d[4][4];
    TH2D *htopCosTheta_ttpT_gen2d[4][4];
    TH2D *hlepCosTheta_ttpT_gen2d[4][4];
    TH2D *hlepPlusCosTheta_ttpT_gen2d[4][4];
    TH2D *hlepMinusCosTheta_ttpT_gen2d[4][4];
    TH2D *hpseudorapiditydiff_ttpT_gen2d[4][4];
    TH2D *hrapiditydiff_ttpT_gen2d[4][4];
    TH2D *hrapiditydiffMarco_ttpT_gen2d[4][4];


    TH2D *hlepChargeAsym_ttRapidity2_gen2d[4][4];
    TH2D *hlepAzimAsym_ttRapidity2_gen2d[4][4];
    TH2D *hlepAzimAsym2_ttRapidity2_gen2d[4][4];
    TH2D *htopSpinCorr_ttRapidity2_gen2d[4][4];
    TH2D *hlepCosOpeningAngle_ttRapidity2_gen2d[4][4];
    TH2D *htopCosTheta_ttRapidity2_gen2d[4][4];
    TH2D *hlepCosTheta_ttRapidity2_gen2d[4][4];
    TH2D *hlepPlusCosTheta_ttRapidity2_gen2d[4][4];
    TH2D *hlepMinusCosTheta_ttRapidity2_gen2d[4][4];
    TH2D *hpseudorapiditydiff_ttRapidity2_gen2d[4][4];
    TH2D *hrapiditydiff_ttRapidity2_gen2d[4][4];
    TH2D *hrapiditydiffMarco_ttRapidity2_gen2d[4][4];


    TH1D *hlepChargeAsymGenDiff[4][4];
    TH1D *hlepAzimAsymGenDiff[4][4];
    TH1D *hlepAzimAsym2GenDiff[4][4];
    TH1D *htopSpinCorrGenDiff[4][4];
    TH1D *hlepCosOpeningAngleGenDiff[4][4];
    TH1D *htopCosThetaGenDiff[4][4];
    TH1D *hlepCosThetaGenDiff[4][4];
    TH1D *hlepPlusCosThetaGenDiff[4][4];
    TH1D *hlepMinusCosThetaGenDiff[4][4];
    TH1D *hpseudorapiditydiffGenDiff[4][4];
    TH1D *hrapiditydiffGenDiff[4][4];
    TH1D *hrapiditydiffMarcoGenDiff[4][4];


    TH2D *htopCosTheta_2d[4][4];
    TH2D *hlepCosTheta_2d[4][4];
    TH2D *hlepChargeAsym_2d[4][4];
    TH2D *hlepAzimAsym_2d[4][4];
    TH2D *htopSpinCorr_2d[4][4];
    TH2D *hlepCosOpeningAngle_2d[4][4];
    TH2D *httMass_2d[4][4];
    TH2D *httpT_2d[4][4];
    TH2D *htopP_2d[4][4];

    TH2D *hmasslb_2d[4][4];


    TH1D *httMassGluongenp[4][4];
    TH1D *httMassQuarkgenp[4][4];
    TH1D *httRapidityGluongenp[4][4];
    TH1D *httRapidityQuarkgenp[4][4];


    TH1D *hlepPlusCosThetaTau_gen[4][4];
    TH1D *hlepMinusCosThetaTau_gen[4][4];

    TH1D *hlepPlusxTau_gen[4][4];
    TH1D *hlepMinusxTau_gen[4][4];

    // For the DY Estimation
    TH1D *hmetInDYEst[4][4];         // MET in Z window for the DY Estimation
    TH1D *hmetOutDYEst[4][4];        // MET outside Z window for the DY Estimation
    TH1D *hdilMassWithMetDYEst[4][4];// Dilepton mass with MET requirement for DY estimation
    TH1D *hdilMassNoMetDYEst[4][4];  // Dilepton mass without MET requirement for DY estimation





    // copied from singleLeptonLooper.h
    Int_t lepPlus_status3_id_;
    Int_t lepMinus_status3_id_;
    Int_t lepPlus_status3_nDaughters_;
    Int_t lepMinus_status3_nDaughters_;
    Int_t nuPlus_status3_id_;
    Int_t nuMinus_status3_id_;

    Int_t   nbs_;
    VofP4 genbs_;
    VofP4 genqgs_;
    Int_t   lep_t_id_;
    Int_t   lep_tbar_id_;
    LorentzVector*  mcnu_;
    LorentzVector*  mclep_;
    Float_t ptzgen_;
    Float_t ptwgen_;

    Int_t   npartons_;
    Int_t   nwzpartons_;
    Float_t maxpartonpt_;
    Float_t ptt_;
    Float_t pttbar_;
    Float_t ptttbar_;
    Float_t mttbar_;
    Float_t etattbar_;
    Float_t rapidityttbar_;

    Float_t m_topminus_gen_;
    Float_t m_topplus_gen_;
    Float_t tt_mass_gen_;
    Float_t ttRapidity_gen_;
    Float_t ttRapidity2_gen_;
    Float_t top_rapiditydiff_cms_gen_;
    Float_t top_pseudorapiditydiff_cms_gen_;
    Float_t top_rapiditydiff_Marco_gen_;
    Float_t tt_pT_gen_;
    Float_t top_costheta_cms_gen_;
    Float_t lep_charge_asymmetry_gen_;
    Float_t lep_azimuthal_asymmetry_gen_;
    Float_t lep_azimuthal_asymmetry2_gen_;
    Float_t lepPlus_costheta_cms_gen_;
    Float_t lepMinus_costheta_cms_gen_;
    Float_t top_spin_correlation_gen_;
    Float_t lep_cos_opening_angle_gen_;

    Float_t m_topminus_gen_origtops_;
    Float_t m_topplus_gen_origtops_;
    Float_t tt_mass_gen_origtops_;
    Float_t ttRapidity_gen_origtops_;
    Float_t ttRapidity2_gen_origtops_;
    Float_t top_rapiditydiff_cms_gen_origtops_;
    Float_t top_pseudorapiditydiff_cms_gen_origtops_;
    Float_t top_rapiditydiff_Marco_gen_origtops_;
    Float_t tt_pT_gen_origtops_;
    Float_t top_costheta_cms_gen_origtops_;
    Float_t lep_charge_asymmetry_gen_origtops_;
    Float_t lep_azimuthal_asymmetry_gen_origtops_;
    Float_t lep_azimuthal_asymmetry2_gen_origtops_;
    Float_t lepPlus_costheta_cms_gen_origtops_;
    Float_t lepMinus_costheta_cms_gen_origtops_;
    Float_t top_spin_correlation_gen_origtops_;
    Float_t lep_cos_opening_angle_gen_origtops_;

    // assorted p4's
    LorentzVector*  t_;   
    LorentzVector*  tbar_;   
    LorentzVector*  b_;   
    LorentzVector*  bbar_;   
    LorentzVector*  lep_t_;   
    LorentzVector*  lep_tbar_;   
    LorentzVector*  stop_t_;   
    LorentzVector*  stop_tbar_;
    LorentzVector*  neutralino_t_;   
    LorentzVector*  neutralino_tbar_;    
    LorentzVector*  genc1_;   
    LorentzVector*  genn2_;
    LorentzVector*  neutralino_c1_;   
    LorentzVector*  neutralino_n2_;    
    LorentzVector*  ttbar_;   
    LorentzVector*  mlep_;   
    LorentzVector*  mclep1_;   
    LorentzVector*  mclep2_;   
    LorentzVector*  mctaud1_;   
    LorentzVector*  mctaud2_;  
    LorentzVector*  mctaudvis1_;   
    LorentzVector*  mctaudvis2_;  

    LorentzVector* lepPlus_status3_;
    LorentzVector* lepMinus_status3_;
    LorentzVector* bPlus_status3_;
    LorentzVector* bMinus_status3_;
    LorentzVector* nuPlus_status3_;
    LorentzVector* nuMinus_status3_;
    LorentzVector* topPlus_status3_;
    LorentzVector* topMinus_status3_;
    LorentzVector* WPlus_status3_;
    LorentzVector* WMinus_status3_;
    LorentzVector* lepPlus_status1_;
    LorentzVector* lepMinus_status1_;
    LorentzVector* bPlus_status1_;
    LorentzVector* bMinus_status1_;
    LorentzVector* nuPlus_status1_;
    LorentzVector* nuMinus_status1_;
    LorentzVector* topPlus_status1_;
    LorentzVector* topMinus_status1_;
    //LorentzVector* WPlus_status1_;
    //LorentzVector* WMinus_status1_;
    LorentzVector* WPlus_status3_orig_;
    LorentzVector* WMinus_status3_orig_;
    LorentzVector* topPlus_status3_orig_;
    LorentzVector* topMinus_status3_orig_;
    LorentzVector* lepPlus_status3_orig_;
    LorentzVector* lepMinus_status3_orig_;
    LorentzVector* nuPlus_status3_orig_;
    LorentzVector* nuMinus_status3_orig_;


    LorentzVector vt_stat1;
    LorentzVector vtbar_stat1;
    LorentzVector vb_stat1;
    LorentzVector vbbar_stat1;
    LorentzVector WPlus_status3_orig;
    LorentzVector WMinus_status3_orig;
    LorentzVector topPlus_status3_orig;
    LorentzVector topMinus_status3_orig;
    LorentzVector lepPlus_status3_orig;
    LorentzVector lepMinus_status3_orig;
    LorentzVector nuPlus_status3_orig;
    LorentzVector nuMinus_status3_orig;
    LorentzVector ttpair;



};



#endif
