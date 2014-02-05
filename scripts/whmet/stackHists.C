#include <vector>
#include <iostream>
#include <string>

#include "TFile.h"
#include "plotUtilities.C"
#include "stackHists.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TString.h"
#include "TMath.h"

bool alreadyInitialized_ = false;
std::string currentPath_;

void initialize( const char* path, bool doData, int doSig ) {

  if (currentPath_.size()) {
    if (currentPath_ == std::string(path)) {
      std::cout << "Using current dir" << std::endl;
      return;
    }

    std::cout << "Re-intializing inputs" << std::endl;

    if (doData && data_) data_->Close();
    wjets_->Close();
    wjetsdd_->Close();
    wbb_->Close();
    zjets_->Close();
    //    tt_->Close();
    //    tt0l_->Close();
    tt1l_->Close();
    tt2l_->Close();
    vv_->Close();
    //    t_->Close();
    ttV_->Close();
    vvv_->Close();
    rare_->Close();
    wzbb_->Close();
    whbb_->Close();
    diltop_->Close();
    lepb_->Close();

    if (doSig > 0) {
      TChiwh_250_1_->Close();
    }
    if (doSig > 1) {
      TChiwh_200_1_->Close();
    }
    if (doSig > 2) {
      TChiwh_300_1_->Close();
    }
    if (doSig > 3) {
      TChiwh_130_1_->Close();
      TChiwh_150_1_->Close();
      //      TChiwh_175_1_->Close();
      // TChiwh_250_25_->Close();
      // TChiwh_250_50_->Close();
      TChiwh_350_1_->Close();
      TChiwh_400_1_->Close();
    }

    mcfiles_.clear();
    mclabels_.clear();

    // should delete old pointers here..
  }

  std::cout << "Reading files from path: " << path << std::endl;

  if (doData) data_ = new TFile(Form("%s/data_histos.root",path));
  // wjets_ = new TFile(Form("%s/wjets_histos.root",path));
  wjets_ = new TFile(Form("%s/wjets_nobb_histos.root",path));
  wjetsdd_ = new TFile(Form("%s/wlight_dd_histos.root",path));
  // wjets_ = new TFile(Form("%s/wjets_comb_histos.root",path));
  //  wbb_ = new TFile(Form("%s/wjets_onlybb_histos.root",path));
  wbb_ = new TFile(Form("%s/wbb_histos.root",path));
  //  wbb_ = new TFile(Form("%s/wbb_nlo_histos.root",path));
  //  wbb_ = new TFile(Form("%s/wbb_all_histos.root",path));
  zjets_ = new TFile(Form("%s/zjets_histos.root",path));
  //  tt_ = new TFile(Form("%s/ttbar_histos.root",path));
  //  tt0l_ = new TFile(Form("%s/ttbar_0l_histos.root",path));
  // tt1l_ = new TFile(Form("%s/ttbar_1l_histos.root",path));
  // tt2l_ = new TFile(Form("%s/ttbar_2l_histos.root",path));
  tt1l_ = new TFile(Form("%s/ttbar_mg_1l_histos.root",path));
  tt2l_ = new TFile(Form("%s/ttbar_mg_2l_histos.root",path));
  vv_ = new TFile(Form("%s/VV_histos.root",path));
  //  t_ = new TFile(Form("%s/single_top_histos.root",path));
  t1l_ = new TFile(Form("%s/single_top_1l_histos.root",path));
  t2l_ = new TFile(Form("%s/single_top_2l_histos.root",path));
  ttV_ = new TFile(Form("%s/ttV_histos.root",path));
  vvv_ = new TFile(Form("%s/VVV_histos.root",path));
  rare_ = new TFile(Form("%s/rare_histos.root",path));
  rareall_ = new TFile(Form("%s/rare_all_histos.root",path));
  wzbb_ = new TFile(Form("%s/wzbb_histos.root",path));
  whbb_ = new TFile(Form("%s/whbb_histos.root",path));
  diltop_ = new TFile(Form("%s/diltop_histos.root",path));
  lepb_ = new TFile(Form("%s/bg_1lb_histos.root",path));
  top1l_ = new TFile(Form("%s/top1l_histos.root",path));

  // if (doSig > 0) {
  //   TChiwh_250_1_ = new TFile(Form("%s/TChiwh_250_1_histos.root",path));
  // }
  // if (doSig > 1) {
  //   TChiwh_200_1_ = new TFile(Form("%s/TChiwh_200_1_histos.root",path));
  // }
  // if (doSig > 2) {
  //   TChiwh_300_1_ = new TFile(Form("%s/TChiwh_300_1_histos.root",path));
  // }
  // if (doSig > 3) {
  //   TChiwh_150_1_ = new TFile(Form("%s/TChiwh_150_1_histos.root",path));
  //   TChiwh_250_25_ = new TFile(Form("%s/TChiwh_250_25_histos.root",path));
  //   TChiwh_250_50_ = new TFile(Form("%s/TChiwh_250_50_histos.root",path));
  //   TChiwh_350_1_ = new TFile(Form("%s/TChiwh_350_1_histos.root",path));
  // }

  // if (doSig > 0) {
  //   TChiwh_200_1_ = new TFile(Form("%s/Wino_200_1_histos.root",path));
  //   TChiwh_250_1_ = new TFile(Form("%s/Wino_250_1_histos.root",path));
  //   TChiwh_300_1_ = new TFile(Form("%s/Wino_300_1_histos.root",path));
  // }
  // if (doSig > 3) {
  //   TChiwh_130_1_ = new TFile(Form("%s/Wino_130_1_histos.root",path));
  //   TChiwh_150_1_ = new TFile(Form("%s/Wino_150_1_histos.root",path));
  //   TChiwh_175_1_ = new TFile(Form("%s/Wino_175_1_histos.root",path));
  //   TChiwh_350_1_ = new TFile(Form("%s/Wino_350_1_histos.root",path));
  //   TChiwh_400_1_ = new TFile(Form("%s/Wino_400_1_histos.root",path));
  // }

  if (doSig > 0) {
    TChiwh_130_1_ = new TFile(Form("%s/TChiWH_130_1_histos.root",path));
    TChiwh_150_1_ = new TFile(Form("%s/TChiWH_150_1_histos.root",path));
    TChiwh_200_1_ = new TFile(Form("%s/TChiWH_200_1_histos.root",path));
    TChiwh_250_1_ = new TFile(Form("%s/TChiWH_250_1_histos.root",path));
    TChiwh_300_1_ = new TFile(Form("%s/TChiWH_300_1_histos.root",path));
    TChiwh_350_1_ = new TFile(Form("%s/TChiWH_350_1_histos.root",path));
    TChiwh_400_1_ = new TFile(Form("%s/TChiWH_400_1_histos.root",path));
  }

  // --- control regions, wlight from MC
  // mcfiles_.push_back(tt2l_);    mclabels_.push_back("ttbar 2l");
  // mcfiles_.push_back(tt1l_);    mclabels_.push_back("ttbar 1l");
  //  mcfiles_.push_back(wzbb_);    mclabels_.push_back("wzbb");
  // mcfiles_.push_back(wjets_); mclabels_.push_back("wlight");
  // mcfiles_.push_back(wbb_); mclabels_.push_back("wbb");
  // mcfiles_.push_back(t1l_);     mclabels_.push_back("single top 1l");
  // mcfiles_.push_back(t2l_);     mclabels_.push_back("single top 2l");
  // mcfiles_.push_back(rare_);     mclabels_.push_back("rare");
  // mcfiles_.push_back(whbb_);     mclabels_.push_back("whbb");

  // --- control regions, wlight from dd
  // mcfiles_.push_back(tt2l_);    mclabels_.push_back("ttbar 2l");
  // mcfiles_.push_back(tt1l_);    mclabels_.push_back("ttbar 1l");
  //  mcfiles_.push_back(wzbb_);    mclabels_.push_back("wzbb");
  // mcfiles_.push_back(wjetsdd_); mclabels_.push_back("wlight");
  // mcfiles_.push_back(wbb_); mclabels_.push_back("wbb");
  // mcfiles_.push_back(t1l_);     mclabels_.push_back("single top 1l");
  // mcfiles_.push_back(t2l_);     mclabels_.push_back("single top 2l");
  // mcfiles_.push_back(rare_);     mclabels_.push_back("rare");
  // mcfiles_.push_back(whbb_);     mclabels_.push_back("whbb");

  // --- signal region, collapsed samples, wlight from dd
  mcfiles_.push_back(diltop_);    mclabels_.push_back("dilep top");
  mcfiles_.push_back(top1l_);    mclabels_.push_back("top 1l");
  mcfiles_.push_back(wzbb_);    mclabels_.push_back("wzbb");
  mcfiles_.push_back(wbb_);    mclabels_.push_back("wbb");
  mcfiles_.push_back(wjetsdd_); mclabels_.push_back("wlight");
  //mcfiles_.push_back(wjets_); mclabels_.push_back("wlight");
  mcfiles_.push_back(rareall_);     mclabels_.push_back("rare");

  //  mcfiles_.push_back(tt0l_);    mclabels_.push_back("ttbar 0l");
  //  mcfiles_.push_back(zjets_); mclabels_.push_back("zjets");
  //  mcfiles_.push_back(vv_);    mclabels_.push_back("VV");
  //  mcfiles_.push_back(ttV_);   mclabels_.push_back("ttV");
  //  mcfiles_.push_back(vvv_);   mclabels_.push_back("VVV");

  // if (doSig == 1) {
  //   mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
  // }
  // else if (doSig == 2) {
  //   mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("TChiwh_200_1");
  //   mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
  // }
  // else if (doSig == 3) {
  //   mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("TChiwh_200_1");
  //   mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
  //   mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("TChiwh_300_1");
  // }
  // else if (doSig > 3) {
  //   mcfiles_.push_back(TChiwh_150_1_);     mclabels_.push_back("TChiwh_150_1");
  //   mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("TChiwh_200_1");
  //   mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("TChiwh_250_1");
  //   mcfiles_.push_back(TChiwh_250_25_);     mclabels_.push_back("TChiwh_250_25");
  //   mcfiles_.push_back(TChiwh_250_50_);     mclabels_.push_back("TChiwh_250_50");
  //   mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("TChiwh_300_1");
  //   mcfiles_.push_back(TChiwh_350_1_);     mclabels_.push_back("TChiwh_350_1");
  // }

  // if (doSig == 1) {
  //   mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
  // }
  // else if (doSig == 2) {
  //   mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
  //   mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("Wino_250_1");
  // }
  // else if (doSig == 3) {
  //   mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
  //   mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("Wino_250_1");
  //   mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("Wino_300_1");
  // }
  // else if (doSig > 3) {
  //   mcfiles_.push_back(TChiwh_130_1_);     mclabels_.push_back("Wino_130_1");
  //   mcfiles_.push_back(TChiwh_150_1_);     mclabels_.push_back("Wino_150_1");
  //   mcfiles_.push_back(TChiwh_175_1_);     mclabels_.push_back("Wino_175_1");
  //   mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
  //   mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("Wino_250_1");
  //   mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("Wino_300_1");
  //   mcfiles_.push_back(TChiwh_350_1_);     mclabels_.push_back("Wino_350_1");
  //   mcfiles_.push_back(TChiwh_400_1_);     mclabels_.push_back("Wino_400_1");
  // }

  if (doSig == 1) {
    mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
  }
  else if (doSig == 2) {
    mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
    mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("Wino_250_1");
  }
  else if (doSig == 3) {
    mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
    //mcfiles_.push_back(TChiwh_150_1_);     mclabels_.push_back("Wino_150_1");
    mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("Wino_250_1");
    //mcfiles_.push_back(TChiwh_350_1_);     mclabels_.push_back("Wino_350_1");
    mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("Wino_300_1");
  }
  else if (doSig > 3) {
    mcfiles_.push_back(TChiwh_130_1_);     mclabels_.push_back("Wino_130_1");
    mcfiles_.push_back(TChiwh_150_1_);     mclabels_.push_back("Wino_150_1");
    mcfiles_.push_back(TChiwh_200_1_);     mclabels_.push_back("Wino_200_1");
    mcfiles_.push_back(TChiwh_250_1_);     mclabels_.push_back("Wino_250_1");
    mcfiles_.push_back(TChiwh_300_1_);     mclabels_.push_back("Wino_300_1");
    mcfiles_.push_back(TChiwh_350_1_);     mclabels_.push_back("Wino_350_1");
    mcfiles_.push_back(TChiwh_400_1_);     mclabels_.push_back("Wino_400_1");
  }

  currentPath_ = std::string(path);
  return;
}

void stackHist( const char* path, const char* hist, const char* flavor, const char* dir, int nbins, float xmin, float xmax, const char* xtitle ) {
  initialize(path);
  TCanvas* c = new TCanvas();
  compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir , 
		 nbins ,  xmin , xmax , xtitle ,
		 true , true , true , false );

  return;
}

//TH1F* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData, int rebin, bool normalize ) {
// TGraphErrors* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData, int rebin, bool normalize, float mcnorm ) {
TCanvas* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData, int rebin, bool normalize, float mcnorm, int doSig, bool doRatio, const char* scalesample, bool stacksig, float signorm, bool errband, const char* cmslabel ) {
  initialize(path, doData, doSig);
  histStyle style = getHistStyle(hist);
  TCanvas* c = new TCanvas();
  compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir ,
  		 style.nbins/rebin ,  style.xmin , style.xmax , (style.xtitle).c_str() ,
		 doData , doData && doRatio , true , style.log, normalize , false , mcnorm, scalesample,
		 stacksig, signorm, errband, cmslabel);
  // compareDataMC( mcfiles_ , mclabels_ , data_ , hist , flavor , dir ,
  // 		 style.nbins/rebin ,  style.xmin , style.xmax , (style.xtitle).c_str() ,
  // 			doData , doData , true , style.log, normalize , false , mcnorm );

  return c;
}

void printYieldsDir( const char* path, const char* dir, bool doData, int latex, int doSig, bool doWeighted ) {
  initialize(path, doData, doSig);
  printYields(mcfiles_, mclabels_, data_, dir, doData, latex, doWeighted);
}

void printCutflowRegion( const char* path, const char* dir, bool doData, int latex, int doSig, bool doWeighted ) {
  initialize(path, doData, doSig);
  std::vector<char*> dirs;
  TString tdir = TString(dir);
  if (tdir.EqualTo("sig")) {
    dirs.push_back("sig_presel");
    dirs.push_back("sig_met_nm1");
    dirs.push_back("sig_mt_nm1");
    dirs.push_back("sig_mt2bl_nm1");
    dirs.push_back("sig_final");
  } else if (tdir.EqualTo("inc")) {
    dirs.push_back("inc_presel");
    dirs.push_back("inc_1b");
    dirs.push_back("inc_2b");
  } else if (tdir.EqualTo("sig_bbmasslast")) {
    dirs.push_back("sig_bbmasslast_presel");
    dirs.push_back("sig_bbmasslast_met_nm1");
    dirs.push_back("sig_bbmasslast_mt_nm1");
    dirs.push_back("sig_bbmasslast_mt2bl_nm1");
    dirs.push_back("sig_bbmasslast_bbmass_nm1");
    dirs.push_back("sig_bbmasslast_final");
  } else if (tdir.EqualTo("sig_metlast")) {
    dirs.push_back("sig_metlast_presel");
    dirs.push_back("sig_metlast_mt2bl_nm1");
    dirs.push_back("sig_metlast_mt_nm1");
    dirs.push_back("sig_metlast_met_nm1");
    dirs.push_back("sig_metlast_final");
  } else if (tdir.EqualTo("cr1")) {
    dirs.push_back("cr1_presel");
    dirs.push_back("cr1_met_nm1");
    dirs.push_back("cr1_mt_nm1");
    dirs.push_back("cr1_mt2bl_nm1");
    dirs.push_back("cr1_final");
  } else if (tdir.EqualTo("cr1_metlast")) {
    dirs.push_back("cr1_metlast_presel");
    dirs.push_back("cr1_metlast_mt2bl_nm1");
    dirs.push_back("cr1_metlast_mt_nm1");
    dirs.push_back("cr1_metlast_met_nm1");
    dirs.push_back("cr1_metlast_met100");
    dirs.push_back("cr1_metlast_met125");
    dirs.push_back("cr1_metlast_met150");
    dirs.push_back("cr1_metlast_final");
  } else if (tdir.EqualTo("cr2")) {
    dirs.push_back("cr2_presel");
    dirs.push_back("cr2_bbmass_nm1");
    dirs.push_back("cr2_mt_nm1"); 
    dirs.push_back("cr2_met_nm1");
    dirs.push_back("cr2_mt2bl_nm1");
    dirs.push_back("cr2_mt2blcut");
    dirs.push_back("cr2_final");
  } else if (tdir.EqualTo("cr3")) {
    dirs.push_back("cr3_presel");
    dirs.push_back("cr3_bbmass_nm1");
    dirs.push_back("cr3_mt_nm1"); 
    dirs.push_back("cr3_met_nm1");
    dirs.push_back("cr3_mt2bl_nm1");
    dirs.push_back("cr3_mt2blcut");
    dirs.push_back("cr3_final");
  } else if (tdir.EqualTo("cr23")) {
    dirs.push_back("cr23_presel");
    dirs.push_back("cr23_bbmass_nm1");
    dirs.push_back("cr23_mt_nm1"); 
    dirs.push_back("cr23_mt2bl_nm1");
    dirs.push_back("cr23_met_nm1");
    dirs.push_back("cr23_met100");
    dirs.push_back("cr23_met125");
    dirs.push_back("cr23_met150");
    dirs.push_back("cr23_final");
  // } else if (tdir.EqualTo("cr2")) {
  //   dirs.push_back("cr2_presel");
  //   dirs.push_back("cr2_bbmass_nm1");
  //   dirs.push_back("cr2_mt2bl_nm1");
  //   dirs.push_back("cr2_met_nm1");
  //   dirs.push_back("cr2_mt_nm1"); 
  //   dirs.push_back("cr2_final");
  // } else if (tdir.EqualTo("cr3")) {
  //   dirs.push_back("cr3_presel");
  //   dirs.push_back("cr3_bbmass_nm1");
  //   dirs.push_back("cr3_mt2bl_nm1");
  //   dirs.push_back("cr3_met_nm1");
  //   dirs.push_back("cr3_mt_nm1"); 
  //   dirs.push_back("cr3_final");
  } else if (tdir.EqualTo("cr5")) {
    dirs.push_back("cr5_presel");
    dirs.push_back("cr5_bbmass_nm1");
    dirs.push_back("cr5_met_nm1");
    dirs.push_back("cr5_mt_nm1"); 
    dirs.push_back("cr5_mt2bl_nm1");
    dirs.push_back("cr5_final");
  } else if (tdir.EqualTo("cr5_metlast")) {
    dirs.push_back("cr5_metlast_presel");
    dirs.push_back("cr5_metlast_bbmass_nm1");
    dirs.push_back("cr5_metlast_mt2bl_nm1");
    dirs.push_back("cr5_metlast_mt_nm1");
    dirs.push_back("cr5_metlast_met_nm1");
    dirs.push_back("cr5_metlast_met100");
    dirs.push_back("cr5_metlast_met125");
    dirs.push_back("cr5_metlast_met150");
    dirs.push_back("cr5_metlast_final");
  } else if (tdir.EqualTo("cr6_metlast")) {
    dirs.push_back("cr6_metlast_presel");
    dirs.push_back("cr6_metlast_bbmass_nm1");
    dirs.push_back("cr6_metlast_mt2bl_nm1");
    dirs.push_back("cr6_metlast_mt_nm1");
    dirs.push_back("cr6_metlast_met_nm1");
    dirs.push_back("cr6_metlast_met100");
    dirs.push_back("cr6_metlast_met150");
    dirs.push_back("cr6_metlast_final");
  } else if (tdir.EqualTo("cr7")) {
    dirs.push_back("cr7_presel");
    dirs.push_back("cr7_met_nm1");
    dirs.push_back("cr7_mt_nm1"); 
    dirs.push_back("cr7_mt2bl_nm1");
    dirs.push_back("cr7_final");
  } else if (tdir.EqualTo("cr8")) {
    dirs.push_back("cr8_presel");
    dirs.push_back("cr8_met_nm1");
    dirs.push_back("cr8_mt_nm1"); 
    dirs.push_back("cr8_mt2bl_nm1");
    dirs.push_back("cr8_final");
  } else if (tdir.EqualTo("cr8_metlast")) {
    dirs.push_back("cr8_metlast_presel");
    dirs.push_back("cr8_metlast_mt2bl_nm1");
    dirs.push_back("cr8_metlast_mt_nm1");
    dirs.push_back("cr8_metlast_met_nm1");
    dirs.push_back("cr8_metlast_met100");
    dirs.push_back("cr8_metlast_met125");
    dirs.push_back("cr8_metlast_met150");
    dirs.push_back("cr8_metlast_final");
  } else if (tdir.EqualTo("cr13")) {
    dirs.push_back("cr13_presel");
    dirs.push_back("cr13_bbmass_nm1");
    dirs.push_back("cr13_mt2bl_nm1");
    dirs.push_back("cr13_met_nm1");
    dirs.push_back("cr13_met100");
    dirs.push_back("cr13_met150");
    dirs.push_back("cr13_final");
  } else if (tdir.EqualTo("cr14")) {
    dirs.push_back("cr14_presel");
    dirs.push_back("cr14_mt2bl_nm1");
    dirs.push_back("cr14_mt_nm1");
    dirs.push_back("cr14_met_nm1");
    dirs.push_back("cr14_met100");
    dirs.push_back("cr14_met125");
    dirs.push_back("cr14_met150");
    dirs.push_back("cr14_final");
  } else {
    std::cout << "WARNING: didn't recognize region: " << dir << std::endl;
    return;
  }

  printCutflow(mcfiles_, mclabels_, data_, dirs, doData, false, latex, doWeighted);
}

void printSignalTable( const char* path, bool doData, int latex, int doSig, bool doWeighted ) {
  initialize(path, doData, doSig);
  std::vector<char*> dirs;
  dirs.push_back("sig_metlast_met100");
  dirs.push_back("sig_metlast_met125");
  dirs.push_back("sig_metlast_met150");
  dirs.push_back("sig_metlast_final");
  // dirs.push_back("sig_2loose_met100");
  // dirs.push_back("sig_2loose_met150");
  // dirs.push_back("sig_2loose_final");

  printSigRegions(mcfiles_, mclabels_, data_, dirs, doData, latex, doWeighted);
}

void printSignalTableExcl( const char* path, bool doData, int latex, int doSig, bool doSystErr ) {
  initialize(path, doData, doSig);

  char* dir = "sig_metlast_met_nm1";

  printSigRegionsExcl(mcfiles_, mclabels_, data_, dir, doData, latex, doSystErr);
}

void saveAllHists( const char* path, const char* dir, const char* flavor, const char* outpath, bool doData, int rebin, bool normalize, float mcnorm ) {

  std::vector<std::string> histnames;

  //  histnames.push_back("h_bbmass");
  histnames.push_back("h_met");
  histnames.push_back("h_lep1mt");
  histnames.push_back("h_mt2bl");
  // histnames.push_back("h_mt2b");
  // histnames.push_back("h_mt2w");
  // histnames.push_back("h_bbpt");
  // histnames.push_back("h_wpt");
  histnames.push_back("h_lep1pt");
  // histnames.push_back("h_lep1metdphi");
  // histnames.push_back("h_njets");
  //  histnames.push_back("h_njetsalleta");
  // histnames.push_back("h_nbjets");
  histnames.push_back("h_jet1pt");
  histnames.push_back("h_jet2pt");
  // if (!TString(dir).Contains("cr5")) {
  //   histnames.push_back("h_bjet1pt");
  //   histnames.push_back("h_bjet2pt");
  // }
  // // histnames.push_back("h_jet1csv");
  // // histnames.push_back("h_jet2csv");
  // if (TString(dir).Contains("cr2")) {
  //   histnames.push_back("h_isotrkpt");
  //   histnames.push_back("h_lep1isotrkdphi");
  // }
  // if (TString(dir).Contains("cr3") || TString(dir).Contains("cr4")) {
  //   //  if (TString(dir).Contains("cr3")) {
  //   //    histnames.clear();
  //   histnames.push_back("h_pseudomet_lep");
  //   histnames.push_back("h_pseudomt_lep");
  //   histnames.push_back("h_pseudomt2bl");
  //   histnames.push_back("h_leppt");
  //   histnames.push_back("h_dphi_pseudomet_lep");
  //   histnames.push_back("h_dildphi");
  // }

  // for investigating met components
  // histnames.push_back("h_met");
  // histnames.push_back("h_trkmet");
  // histnames.push_back("h_met_soft");
  // histnames.push_back("h_sumet");
  // histnames.push_back("h_sumet_soft");
  // histnames.push_back("h_ht");
  // histnames.push_back("h_htlep");
  // histnames.push_back("h_ljetspt");
  // histnames.push_back("h_bbpt");
  // histnames.push_back("h_lep1pt");
  // histnames.push_back("h_jet1pt");
  // histnames.push_back("h_jet2pt");
  // histnames.push_back("h_rhovor");
  // if (!TString(dir).Contains("cr5") && !TString(dir).Contains("cr11")) {
  //   histnames.push_back("h_bjet1pt");
  //   histnames.push_back("h_bjet2pt");
  //   histnames.push_back("h_lbbpt");
  // }

  for (unsigned int i = 0; i < histnames.size(); ++i) {
    std::string label_flavor = std::string(flavor);
    if (std::string(flavor).length() > 0) histnames[i] += std::string("_") + std::string(flavor);
    else if (TString(dir).Contains("cr3") || TString(dir).Contains("cr4")) label_flavor = "all";
    TCanvas* c = stackHistAuto(path, histnames[i].c_str(), label_flavor.c_str(), dir, doData, rebin, normalize, mcnorm);
    c->SaveAs(Form("%s/%s_%s.eps",outpath,dir,TString(histnames[i]).ReplaceAll("h_","").Data()));
  }

}

void savePublicPlotsCRs(const char* path, const char* outpath) {

  // !! these should use MC W+light jets

  // CR0b hists
  TCanvas* c = stackHistAuto(path,"h_bbmass","","cr5_metlast_bbmass_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crzb_presel_mbb.eps",outpath));
  c->SaveAs(Form("%s/crzb_presel_mbb.pdf",outpath));
  c->SaveAs(Form("%s/crzb_presel_mbb.png",outpath));

  c = stackHistAuto(path,"h_mt2bl","","cr5_metlast_mt2bl_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crzb_presel_mt2bl.eps",outpath));
  c->SaveAs(Form("%s/crzb_presel_mt2bl.pdf",outpath));
  c->SaveAs(Form("%s/crzb_presel_mt2bl.png",outpath));

  c = stackHistAuto(path,"h_lep1mt","","cr5_metlast_mt2bl_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crzb_presel_mt.eps",outpath));
  c->SaveAs(Form("%s/crzb_presel_mt.pdf",outpath));
  c->SaveAs(Form("%s/crzb_presel_mt.png",outpath));

  c = stackHistAuto(path,"h_met","","cr5_metlast_mt2bl_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crzb_presel_met.eps",outpath));
  c->SaveAs(Form("%s/crzb_presel_met.pdf",outpath));
  c->SaveAs(Form("%s/crzb_presel_met.png",outpath));

  // CR2l hists
  c = stackHistAuto(path,"h_bbmass","","cr23_bbmass_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crdl_presel_mbb.eps",outpath));
  c->SaveAs(Form("%s/crdl_presel_mbb.pdf",outpath));
  c->SaveAs(Form("%s/crdl_presel_mbb.png",outpath));

  c = stackHistAuto(path,"h_mt2bl","","cr23_mt_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crdl_presel_mt2bl.eps",outpath));
  c->SaveAs(Form("%s/crdl_presel_mt2bl.pdf",outpath));
  c->SaveAs(Form("%s/crdl_presel_mt2bl.png",outpath));

  c = stackHistAuto(path,"h_lep1mt","","cr23_mt_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crdl_presel_mt.eps",outpath));
  c->SaveAs(Form("%s/crdl_presel_mt.pdf",outpath));
  c->SaveAs(Form("%s/crdl_presel_mt.png",outpath));

  c = stackHistAuto(path,"h_met","","cr23_mt_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/crdl_presel_met.eps",outpath));
  c->SaveAs(Form("%s/crdl_presel_met.pdf",outpath));
  c->SaveAs(Form("%s/crdl_presel_met.png",outpath));


}

void savePublicPlotsCRMbb(const char* path, const char* outpath) {

  // CRMbb mass hists (should make with scale factors applied)
  TCanvas* c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_mt_nm1",1,1,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/sig_bbmasslast_lepbsfs_mt_nm1_bbmass.eps",outpath));

  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met_nm1",1,2,0,-1,0,1,"",0,-1,0,"Unpublished");
  c->SaveAs(Form("%s/sig_bbmasslast_lepbsfs_met_nm1_bbmass.eps",outpath));

  // c = stackHistAuto(path,"h_mt2bl","","cr14_mt2bl_nm1",1,1,0,-1,0,1);
  // c->SaveAs(Form("%s/cr14_mt2bl_nm1_mt2bl.eps",outpath));

  // c = stackHistAuto(path,"h_lep1mt","","cr14_mt2bl_nm1",1,1,0,-1,0,1);
  // c->SaveAs(Form("%s/cr14_mt2bl_nm1_lep1mt.eps",outpath));

  // c = stackHistAuto(path,"h_met","","cr14_mt2bl_nm1",1,1,0,-1,0,1);
  // c->SaveAs(Form("%s/cr14_mt2bl_nm1_met.eps",outpath));

}

void savePublicPlotsSig1(const char* path, const char* outpath) {

  // sig presel plots, bkg norm to data, with 1 sig overlaid, scaled up
  TCanvas* c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_mt2bl_nm1",1,1,1,-1,1,1,"",0,20.,0,"Unpublished");
  c->SaveAs(Form("%s/sig_presel_mbb.eps",outpath));
  c->SaveAs(Form("%s/sig_presel_mbb.pdf",outpath));
  c->SaveAs(Form("%s/sig_presel_mbb.png",outpath));

  c = stackHistAuto(path,"h_mt2bl","","sig_metlast_mt2bl_nm1",1,1,1,-1,1,1,"",0,10.,0,"Unpublished");
  c->SaveAs(Form("%s/sig_presel_mt2bl.eps",outpath));
  c->SaveAs(Form("%s/sig_presel_mt2bl.pdf",outpath));
  c->SaveAs(Form("%s/sig_presel_mt2bl.png",outpath));

  c = stackHistAuto(path,"h_lep1mt","","sig_metlast_mt2bl_nm1",1,1,1,-1,1,1,"",0,10.,0,"Unpublished");
  c->SaveAs(Form("%s/sig_presel_mt.eps",outpath));
  c->SaveAs(Form("%s/sig_presel_mt.pdf",outpath));
  c->SaveAs(Form("%s/sig_presel_mt.png",outpath));

  c = stackHistAuto(path,"h_met","","sig_metlast_mt2bl_nm1",1,1,1,-1,1,1,"",0,10.,0,"Unpublished");
  c->SaveAs(Form("%s/sig_presel_met.eps",outpath));
  c->SaveAs(Form("%s/sig_presel_met.pdf",outpath));
  c->SaveAs(Form("%s/sig_presel_met.png",outpath));

}

void savePublicPlotsSig3(const char* path, const char* outpath) {

  // sig mass hists with signals overlaid
  TCanvas* c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met100",1,2,0,-1,3,0,"",1);
  c->SaveAs(Form("%s/sig_met100_mbb.eps",outpath));

  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met125",1,2,0,-1,3,0,"",1);
  c->SaveAs(Form("%s/sig_met125_mbb.eps",outpath));

  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met150",1,2,0,-1,3,0,"",1);
  c->SaveAs(Form("%s/sig_met150_mbb.eps",outpath));

  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_bbmass_nm1",1,2,0,-1,3,0,"",1);
  c->SaveAs(Form("%s/sig_met175_mbb.eps",outpath));

  // .. and again with error bands.  Currently stacked version
  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met100",1,2,0,-1,3,0,"",1,-1,1,"");
  //c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met100",1,2,0,-1,3,0,"",0,-1,1,"");
  c->SaveAs(Form("%s/sig_met100_mbb_errband.eps",outpath));
  c->SaveAs(Form("%s/sig_met100_mbb_errband.pdf",outpath));
  c->SaveAs(Form("%s/sig_met100_mbb_errband.png",outpath));

  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met125",1,2,0,-1,3,0,"",1,-1,1,"");
  //c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met125",1,2,0,-1,3,0,"",0,-1,1,"");
  c->SaveAs(Form("%s/sig_met125_mbb_errband.eps",outpath));
  c->SaveAs(Form("%s/sig_met125_mbb_errband.pdf",outpath));
  c->SaveAs(Form("%s/sig_met125_mbb_errband.png",outpath));

  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met150",1,2,0,-1,3,0,"",1,-1,1,"");
  //c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_met150",1,2,0,-1,3,0,"",0,-1,1,"");
  c->SaveAs(Form("%s/sig_met150_mbb_errband.eps",outpath));
  c->SaveAs(Form("%s/sig_met150_mbb_errband.pdf",outpath));
  c->SaveAs(Form("%s/sig_met150_mbb_errband.png",outpath));

  c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_bbmass_nm1",1,2,0,-1,3,0,"",1,-1,1,"");
  //c = stackHistAuto(path,"h_bbmass","","sig_bbmasslast_bbmass_nm1",1,2,0,-1,3,0,"",0,-1,1,"");
  c->SaveAs(Form("%s/sig_met175_mbb_errband.eps",outpath));
  c->SaveAs(Form("%s/sig_met175_mbb_errband.pdf",outpath));
  c->SaveAs(Form("%s/sig_met175_mbb_errband.png",outpath));

}



histStyle getHistStyle( const char* hist ) {

  histStyle style;
  std::string histString = std::string(hist);

  if (matchHistName(histString,"h_ht30")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep1pt")) {
    //    style.nbins = 15;
    style.nbins = 12;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep1mt")) {
    style.nbins = 14;
    style.xmin = 0.;
    style.xmax = 350.;
    // style.nbins = 20;
    // style.xmin = 0.;
    //    style.xmax = 500.;
    style.xtitle = "M_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_mt2bl")) {
    style.nbins = 14;
    style.xmin = 0.;
    style.xmax = 350.;
    // style.nbins = 20;
    // style.xmin = 0.;
    //    style.xmax = 500.;
    style.xtitle = "M_{T2}^{bl} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet1pt")) {
    // style.nbins = 25;
    // style.xmin = 0.;
    //    style.xmax = 500.;
    style.nbins = 14;
    style.xmin = 0.;
    style.xmax = 350.;
    style.xtitle = "p_{T}(jet1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet2pt")) {
    // style.nbins = 15;
    // style.xmin = 0.;
    // style.xmax = 300.;
    style.nbins = 12;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(jet2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_met")) {
    style.nbins = 14;
    style.xmin = 0.;
    style.xmax = 350.;
    // style.nbins = 20;
    // style.xmin = 0.;
    //    style.xmax = 500.;
    // style.nbins = 40;
    // style.xmin = 0.;
    // style.xmax = 1000.;
    style.xtitle = "E_{T}^{miss} [GeV]";
    //    style.xtitle = "type 1 phi cor pf MET [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_metsignif")) {
    style.nbins = 10;
    style.xmin = 0.;
    style.xmax = 250.;
    style.xtitle = "E_{T}^{miss} significance";
    style.log = true;
  } else if (matchHistName(histString,"h_bdt0")) {
    style.nbins = 24;
    style.xmin = -0.2;
    style.xmax = 1.;
    style.xtitle = "BDT 0 Weight";
    style.log = true;
  } else if (matchHistName(histString,"h_bdt1")) {
    style.nbins = 24;
    style.xmin = -0.2;
    style.xmax = 1.;
    style.xtitle = "BDT 1 Weight";
    style.log = true;
  } else if (matchHistName(histString,"h_bdt2")) {
    style.nbins = 24;
    style.xmin = -0.2;
    style.xmax = 1.;
    style.xtitle = "BDT 2 Weight";
    style.log = true;
  } else if (matchHistName(histString,"h_bdt3")) {
    style.nbins = 24;
    style.xmin = -0.2;
    style.xmax = 1.;
    style.xtitle = "BDT 3 Weight";
    style.log = true;
  } else if (matchHistName(histString,"h_bdt4")) {
    style.nbins = 24;
    style.xmin = -0.2;
    style.xmax = 1.;
    style.xtitle = "BDT 4 Weight";
    style.log = true;
  } else if (matchHistName(histString,"h_mct")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    // style.nbins = 14;
    // style.xmin = 0.;
    //    style.xmax = 350.;
    style.xtitle = "M_{CT} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ht")) {
    style.nbins = 32;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "H_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_htlep")) {
    style.nbins = 32;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "H_{T}+p_{T}(lep) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ht40")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ht40up")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T}^{up} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ht40dn")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "H_{T}^{dn} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilmass")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 200.;
    style.xtitle = "m(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    // style.nbins = 160;
    // style.xmin = 0.;
    // style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt_ee")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt_mm")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(ll) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dileta")) {
    style.nbins = 30;
    style.xmin = -6.0;
    style.xmax = 6.0;
    style.xtitle = "#eta(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dilrapidity")) {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dilrapidity_lowdilpt")) {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dilrapidity_highdilpt")) {
    // style.nbins = 25;
    // style.xmin = -3.0;
    // style.xmax = 3.0;
    style.nbins = 21;
    style.xmin = -2.52;
    style.xmax = 2.52;
    style.xtitle = "y(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dildeltar")) {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_dildeltaphi")) {
    style.nbins = 50;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(ll)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep2pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep3pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep3) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep2mt")) {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "m_{T}(lep2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep3mt")) {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "m_{T}(lep3) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lep1phi")) {
    style.nbins = 25;
    style.xmin = -TMath::Pi();
    style.xmax = TMath::Pi();
    style.xtitle = "#phi(lep1)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep2phi")) {
    style.nbins = 25;
    style.xmin = -TMath::Pi();
    style.xmax = TMath::Pi();
    style.xtitle = "#phi(lep2)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep1isopf")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 0.5;
    style.xtitle = "isopf(lep1) [GeV]";
    style.log = false;
  } else if (matchHistName(histString,"h_lep1metdphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(lep1,MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1metdphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(jet1,MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet2metdphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(jet2,MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet1metdphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(bjet1,MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet2metdphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(bjet2,MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_isotrkpt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(iso track) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_mt2b")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "M_{T2}^{b} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_mt2w")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "M_{T2}^{W} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet1eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1rapidity")) {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1eta_lowdilpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1rapidity_lowdilpt")) {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1eta_highdilpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1rapidity_highdilpt")) {
    style.nbins = 25;
    style.xmin = -3.0;
    style.xmax = 3.0;
    style.xtitle = "y(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet3pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(jet3) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet4pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(jet4) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet1pumva")) {
    style.nbins = 25;
    style.xmin = -1.0;
    style.xmax = 1.0;
    style.xtitle = "PU MVA(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet2pumva")) {
    style.nbins = 25;
    style.xmin = -1.0;
    style.xmax = 1.0;
    style.xtitle = "PU MVA(jet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_fwdjet1pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(fwd jet1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_fwdjet2pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(fwd jet2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bjet1pt")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(bjet1) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bjet2pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(bjet2) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bjet3pt")) {
    style.nbins = 15;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(bjet3) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_jet1eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet2eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet3eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(jet3)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet1eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet2eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet3eta")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet3)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet1eta_lowpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet2eta_lowpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet1eta_highpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_bjet2eta_highpt")) {
    style.nbins = 25;
    style.xmin = -2.5;
    style.xmax = 2.5;
    style.xtitle = "#eta(bjet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_fwdjet1eta")) {
    style.nbins = 25;
    style.xmin = -5.0;
    style.xmax = 5.0;
    style.xtitle = "#eta(fwd jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_fwdjet2eta")) {
    style.nbins = 25;
    style.xmin = -5.0;
    style.xmax = 5.0;
    style.xtitle = "#eta(fwd jet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_fwdjet1pumva")) {
    style.nbins = 25;
    style.xmin = -1.0;
    style.xmax = 1.0;
    style.xtitle = "PU MVA(fwd jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_fwdjet2pumva")) {
    style.nbins = 25;
    style.xmin = -1.0;
    style.xmax = 1.0;
    style.xtitle = "PU MVA(fwdjet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_bbpt")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(b#bar{b}) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_bbmass")) {
    // style.nbins = 50;
    // style.xmin = 0.;
    // style.xmax = 1000.;
    // style.xtitle = "M(b#bar{b}) [GeV]";
    // style.log = true;
    //    style.nbins = 25;
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "M_{b#bar{b}} [GeV]";
    style.log = false;
  } else if (matchHistName(histString,"h_bbmdrpt")) {
    style.nbins = 25;
    style.xmin = 0.0;
    //    style.xmax = 2*TMath::Pi();
    style.xmax = 5.;
    style.xtitle = "M(b#bar{b})*#DeltaR(b#bar{b})/p_{T}(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdr")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdphi")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdeta")) {
    style.nbins = 20;
    style.xmin = -6.;
    style.xmax = 6.;
    style.xtitle = "#Delta#eta(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdr_lowpt")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdeta_lowpt")) {
    style.nbins = 20;
    style.xmin = -6.;
    style.xmax = 6.;
    style.xtitle = "#Delta#eta(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdr_highpt")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_bbdeta_highpt")) {
    style.nbins = 20;
    style.xmin = -6.;
    style.xmax = 6.;
    style.xtitle = "#Delta#eta(b#bar{b})";
    style.log = false;
  } else if (matchHistName(histString,"h_jjdr")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(jet1,jet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_jjdphi")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(jet1,jet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep1bjet1dr")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(lep1,bjet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep1bjet2dr")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#DeltaR(lep1,bjet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_wpt")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(W) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ljetspt")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(l+jets) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_lbbpt")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(lbb) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_simplemetlbb")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(lbb) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_pfmet")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pf MET, uncorrected [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_trkmet")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "track MET [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_met_soft")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pf MET, soft terms [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_t1met10")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "type 1 pf MET [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumet")) {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 3000.;
    style.xtitle = "type1 pf sum E_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumet_soft")) {
    style.nbins = 60;
    style.xmin = 0.;
    style.xmax = 3000.;
    style.xtitle = "pf sum E_{T}, soft terms [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_metsig")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 15.;
    style.xtitle = "pf MET / #sqrt{sum E_{T}}";
    style.log = true;
  } else if (matchHistName(histString,"h_ngoodlep")) {
    style.nbins = 5;
    style.xmin = 0;
    style.xmax = 5;
    style.xtitle = "N(good leptons)";
    style.log = false;
  } else if (matchHistName(histString,"h_njets")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets, |#eta| < 2.4)";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsfwd")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets, 2.4 < |#eta| < 4.7)";
    style.log = true;
  } else if (matchHistName(histString,"h_njets20")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets, 20 < p_{T} < 30 GeV && |#eta| < 2.4)";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsall")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets)";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsalleta")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets, |#eta| < 4.7)";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsup")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets), JES up";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsdn")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets), JES down";
    style.log = true;
  } else if (matchHistName(histString,"h_njetsupdiff")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets) up - N(jets)";
    style.log = true;
  } else if (matchHistName(histString,"h_nbjets")) {
    style.nbins = 5;
    style.xmin = 0;
    style.xmax = 5;
    style.xtitle = "N(b jets)";
    style.log = true;
  } else if (matchHistName(histString,"h_nnonbjets")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets)";
    style.log = true;
  } else if (matchHistName(histString,"h_nnonbjetsup")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets), JES up";
    style.log = true;
  } else if (matchHistName(histString,"h_nnonbjetsdn")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(nonbtagged jets), JES down";
    style.log = true;
  } else if (matchHistName(histString,"h_njets40")) {
    style.nbins = 8;
    style.xmin = 0;
    style.xmax = 8;
    style.xtitle = "N(jets)";
    style.log = true; 
  } else if (matchHistName(histString,"h_nbcsvm")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(btags, csv Medium)";
    style.log = true;
  } else if (matchHistName(histString,"h_nbcsvl")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(btags, csv Loose)";
    style.log = true;
  } else if (matchHistName(histString,"h_nmistags")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(mistags, csv Medium)";
    style.log = true;
  } else if (matchHistName(histString,"h_ngenb")) {
    style.nbins = 4;
    style.xmin = 0;
    style.xmax = 4;
    style.xtitle = "N(reco jets matched to gen b quarks)";
    style.log = true;
  } else if (matchHistName(histString,"h_nlep")) {
    style.nbins = 5;
    style.xmin = 0;
    style.xmax = 5;
    style.xtitle = "N(lep)";
    style.log = true;
  } else if (matchHistName(histString,"h_sumpt")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumpt_ee")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumpt_mm")) {
    // style.nbins = 100;
    // style.xmin = 0.;
    // style.xmax = 500.;
    style.nbins = 160;
    style.xmin = 0.;
    style.xmax = 800.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptup")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptdn")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptup10")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up x 10 [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_sumptdn10")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down x 10 [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ptpar")) {
    style.nbins = 100;
    style.xmin = -100.;
    style.xmax = 400.;
    style.xtitle = "jet vec sum p_{T}^{#parallel} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_ptperp")) {
    style.nbins = 100;
    style.xmin = -250.;
    style.xmax = 250.;
    style.xtitle = "jet vec sum p_{T}^{#perp} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_nonbsumpt")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_nonbsumptup")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES up [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_nonbsumptdn")) {
    style.nbins = 100;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "p_{T}(jet system), JES down [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dilpt_minus_ptpar_b2b")) {
    style.nbins = 100;
    style.xmin = -250.;
    style.xmax = 250.;
    style.xtitle = "dilep p_{T} - jet vec sum p_{T}^{#parallel} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_deltaphi_dil_sum")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(dilep, jet vec sum)";
    style.log = false;
  } else if (matchHistName(histString,"h_deltaphi_z_sum")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(gen Z, jet system)";
    style.log = false;
  } else if (matchHistName(histString,"h_deltaphi_ttbar_nonbsum")) {
    style.nbins = 25;
    style.xmin = 0.;
    style.xmax = 2*TMath::Pi();
    style.xtitle = "#Delta#phi(true ttbar, non-btagged jet vec sum)";
    style.log = false;
  } else if (matchHistName(histString,"h_rhovor")) {
    style.nbins = 80;
    style.xmin = 0.;
    style.xmax = 40;
    style.xtitle = "#rho [GeV/unit area]";
    style.log = true;
  } else if (matchHistName(histString,"h_nvtx")) {
    style.nbins = 40;
    style.xmin = 0.;
    style.xmax = 40;
    style.xtitle = "N(vtx)";
    style.log = false;
  } else if (matchHistName(histString,"h_vtx")) {
    style.nbins = 40;
    style.xmin = 0.;
    style.xmax = 40;
    style.xtitle = "N(vtx)";
    style.log = false;
  } else if (matchHistName(histString,"h_pseudomet_lep")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pseudo MET [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_pseudomt_lep")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pseudo m_{T} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_pseudomt2bl")) {
    style.nbins = 20;
    style.xmin = 0.;
    style.xmax = 500.;
    style.xtitle = "pseudo M_{T2}^{bl} [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_leppt")) {
    style.nbins = 30;
    style.xmin = 0.;
    style.xmax = 300.;
    style.xtitle = "p_{T}(lep) [GeV]";
    style.log = true;
  } else if (matchHistName(histString,"h_dildphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(lep1,lep2)";
    style.log = false;
  } else if (matchHistName(histString,"h_dphi_pseudomet_lep")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(lep,pseudo MET)";
    style.log = false;
  } else if (matchHistName(histString,"h_lep1isotrkdphi")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(lep1,iso track)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet1csv")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 1.0;
    style.xtitle = "CSV(jet1)";
    style.log = false;
  } else if (matchHistName(histString,"h_jet2csv")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 1.0;
    style.xtitle = "CSV(jet2)";
    style.log = false;
  } else if (matchHistName(histString,"h_jetcsv")) {
    style.nbins = 20;
    style.xmin = 0.0;
    style.xmax = 1.0;
    style.xtitle = "CSV";
    style.log = false;
  } else if (matchHistName(histString,"h_pfmetcalometdphi_met100")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(pf MET, calo MET), MET > 100";
    style.log = true;
  } else if (matchHistName(histString,"h_pfmetcalometdphi_met150")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(pf MET, calo MET), MET > 150";
    style.log = true;
  } else if (matchHistName(histString,"h_pfmetcalometdphi_met175")) {
    style.nbins = 10;
    style.xmin = 0.0;
    style.xmax = TMath::Pi();
    style.xtitle = "#Delta#phi(pf MET, calo MET), MET > 175";
    style.log = true;
  } else if (matchHistName(histString,"h_passtauveto")) {
    style.nbins = 2;
    style.xmin = 0.;
    style.xmax = 2.;
    style.xtitle = "pass tau veto";
    style.log = false;

  } else {
    std::cout << "ERROR: hist name " << histString << " not recognized" << std::endl;
  }

  return style;
}

bool matchHistName(const TString& histname, const TString& matchname) {

  // first check for a simple match
  if (histname.EqualTo(matchname)) return true;

  // otherwise check for a name match with a channel suffix
  if (histname.Contains(matchname)) {
    TString suffix = histname;
    suffix.ReplaceAll(matchname,"");
    if (suffix.EqualTo("_e")) return true;
    if (suffix.EqualTo("_m")) return true;
    if (suffix.EqualTo("_ee")) return true;
    if (suffix.EqualTo("_mm")) return true;
    if (suffix.EqualTo("_em")) return true;
    if (suffix.EqualTo("_lowpu")) return true;
    if (suffix.EqualTo("_highpu")) return true;
    if (matchname.Contains("jetcsv") && suffix.Contains("_pt")) return true;
  }

  return false;
}

bool matchHistName(const TString& histname, const char* matchname) {
  return matchHistName(histname,TString(matchname));
}
