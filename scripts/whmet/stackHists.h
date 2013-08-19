#ifndef STACKHISTS
#define STACKHISTS

#include "TFile.h"
#include <vector>
#include <string>
#include "plotUtilities.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"

struct histStyle {
  int nbins;
  float xmin;
  float xmax;
  std::string xtitle;
  bool log;
};

TFile *data_;
TFile *tt_;
TFile *tt0l_;
TFile *tt1l_;
TFile *tt2l_;
TFile *wjets_;
TFile *wjetsdd_;
TFile *wbb_;
TFile *zjets_;
TFile *vv_;
TFile *t_;
TFile *t1l_;
TFile *t2l_;
TFile *ttV_;
TFile *vvv_;
TFile *rare_;
TFile *rareall_;
TFile *wzbb_;
TFile *whbb_;
TFile *diltop_;
TFile *lepb_;
TFile *top1l_;

TFile *TChiwh_130_1_;
TFile *TChiwh_150_1_;
TFile *TChiwh_175_1_;
TFile *TChiwh_200_1_;
TFile *TChiwh_250_1_;
TFile *TChiwh_250_25_;
TFile *TChiwh_250_50_;
TFile *TChiwh_300_1_;
TFile *TChiwh_350_1_;
TFile *TChiwh_400_1_;

vector<TFile*>  mcfiles_;
vector<char*>   mclabels_;

histStyle getHistStyle(const char* hist);
void initialize( const char* path, bool doData = true, int doSig = 0 );
//TH1F* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData = true, int rebin = 1, bool normalize = false );
// TGraphErrors* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData = true, 
//                              int rebin = 1, bool normalize = false, float mcnorm = -1. );
TCanvas* stackHistAuto( const char* path, const char* hist, const char* flavor, const char* dir, bool doData = true, 
			int rebin = 1, bool normalize = false, float mcnorm = -1., int doSig = 0 , bool doRatio = true, const char* scalesample = "",
			bool stacksig = false, float signorm = -1., bool errband = false);
void printYieldsDir( const char* path, const char* dir, bool doData, int latex = 0, int doSig = 9, bool doWeighted = true );
void saveAllHists( const char* path, const char* dir, const char* flavor, const char* outpath, bool doData, int rebin = 1, bool normalize = false, float mcnorm = -1. );
bool matchHistName(const TString& histname, const TString& matchname);
bool matchHistName(const TString& histname, const char* matchname);
void printCutflowRegion( const char* path, const char* dir, bool doData, int latex, int doSig, bool doWeighted = true );
void printSignalTable( const char* path, bool doData, int latex, int doSig, bool doWeighted = true );

#endif
