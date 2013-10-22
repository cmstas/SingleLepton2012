#ifndef WHLOOPER_H
#define WHLOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include "Math/LorentzVector.h"
 
#include <cmath>
#include <map>

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace std;

class WHLooper {

 public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

  enum csvpoint { CSVL, CSVM, CSVT };

  WHLooper();
  ~WHLooper();

  void setOutFileName(string filename); 
  void loop(TChain *chain, TString name);

 private:

  void fillHists1DWrapper(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "");
  void fillHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");
  void fillHists2D(std::map<std::string, TH2F*>& h_2d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");
  void fillFlavorHists1D(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");
  void fillJetAccHists(std::map<std::string, TH1F*>& h_1d, const float evtweight = 1., const std::string& dir = "", const std::string& suffix = "");

  float getCSVCut(const csvpoint csv = WHLooper::CSVM);

  void dumpEventInfo(const std::string& comment);

  string m_outfilename_;
  TFile* outfile_;
  //for phi corrected met
  float t1metphicorr;
  float t1metphicorrphi;
  float t1metphicorrmt;

  // store main cut values
  float CUT_BBMASS_LOW_;
  float CUT_BBMASS_HIGH_;
  float CUT_BBMASS_CR1_LOW_;
  float CUT_BBMASS_CR1_HIGH_;
  float CUT_BBMASS_CR8_LOW_;
  float CUT_BBMASS_CR8_HIGH_;
  float CUT_MET_PRESEL_;
  float CUT_MET_;
  float CUT_MT_PRESEL_;
  float CUT_MT_;
  float CUT_MT_CR13_LOW_;
  float CUT_MT_CR13_HIGH_;
  float CUT_MT2BL_;

  // internal vars to store
  std::vector<LorentzVector> jets_;
  std::vector<LorentzVector> bjets_;
  std::vector<LorentzVector> jets_fwd_;
  std::vector<float> jets_csv_;
  std::vector<int> jets_idx_;
  std::vector<int> jets_fwd_idx_;
  std::vector<int> bjets_idx_;
  std::vector<float> jets_smearcorrs_;
  Int_t njets_;
  Int_t njetsfwd_;
  Int_t njets20_;
  Int_t njetsall_;
  Int_t nbjets_;
  Int_t nbjetst_;
  Int_t nbjetsl_;
  Float_t met_;
  Float_t metphi_;
  Float_t mt_;
  Float_t mt2b_;
  Float_t mt2bl_;
  Float_t mt2w_;
  Float_t mct_;
  LorentzVector bb_;
  Float_t met_soft_;
  Float_t sumet_;
  Float_t sumet_soft_;
  Float_t ht_;
  Float_t wpt_;
  Float_t lepmetdphi_;
  Float_t bbwdphi_;
  Float_t bbpt_;
  Float_t bbdR_;
  Float_t pt_J1_;
  Float_t pt_J2_;
  Float_t lep1pt_;
  // for ttbar samples
  Float_t genmt2bl_;
  bool tobtecveto_;
  // to hold BDT vals
  std::vector<float> bdtvals_;

  // for CR3 (with two leptons)
  LorentzVector lep_;
  Float_t pseudomet_lep_;
  Float_t pseudometphi_lep_;
  Float_t pseudomt_lep_;
  Float_t dphi_pseudomet_lep_;
  Float_t pseudomt2b_;
  Float_t pseudomt2bl_;
  Float_t pseudomt2w_;

  // internal flags
  bool isWjets_;
  bool isWNjets_;
  bool isWNjets_nobb_;
  bool isWNjets_onlybb_;
  bool isWbbMG_;
  bool isWbbMCatNLO_;
  bool isWZlight_;
  bool isWZbb_;
  bool isWHbb_;
  bool isWino_;
  bool isTChiwh_;
  bool isTChiWHMG_;
  bool isTChihhwwbb_;
  bool isScan_;
  bool isttmg_;
  bool isttsl_;
  bool isttdl_;
  bool istsl_;
  bool istdl_;
  bool isttvsl_;
  bool isttvdl_;
  bool isttvother_;

};

#endif
