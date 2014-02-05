#ifndef PLOTUTILITIES
#define PLOTUTILITIES

#include "TFile.h"
#include "TChain.h"
#include "TLegend.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include <vector>
#include <string>

struct dataMCHists{
  TPad* comppad;
  TH1F* datahist;
  TH1F* mctothist;
  TH1F* ratio;
};

//TH1F* compareDataMC( std::vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , const char* flavor , const char* dir ,
TGraphErrors* compareDataMC( std::vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , const char* flavor , const char* dir ,
			     int nbins ,  float xmin , float xmax , const char* xtitle , 
			     bool overlayData = true , bool residual = false, bool drawLegend = true , 
			     bool log = false , bool normalize = false , bool fit = false, float mcnorm = -1., const char* scalesample = "",
			     bool stacksig = false, float signorm = -1., bool errband = false, const char* cmslabel = "Preliminary");

void initSymbols(int latex);
void  printLine(int latex);
void printHeader(bool dilep = false);
void print( TH1F* h , string label , bool dilep = false, bool correlatedError = false );
void printYields( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* dir , bool doData, int latex = 0, bool doWeighted = true );
void fillYieldHist( TH1F* hin, TH1F* hyield, int bin, bool doWeighted = true );
void printCutflow( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , vector<char*> dirs , bool doData, bool splitMC = false, int latex = 0, bool doWeighted = true );
void printSigRegions( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , vector<char*> dirs , bool doData, int latex = 0, bool doWeighted = true );
void printSigRegionsExcl( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , char* dir , bool doData, int latex = 0, bool doSystErr = false );
float getSystError( const TString sample );
std::string getTableName( const TString sample, int latex = 0 );

//void deleteHistos();
int getColor(const TString sample);
TLegend *getLegend( vector<char*> labels , bool overlayData, 
		    float x1 = 0.7, float y1 = 0.45 , float x2 = 0.87 , float y2 = 0.94, float signorm = -1., bool errband = false );

TGraphAsymmErrors* makeBand(TH1F* centhist, TH1F* uphist, TH1F* dnhist);
TCanvas* compareNormalized(std::string histname, TFile* file1, std::string label1, TFile* file2, std::string label2, int rebin = 1, bool norm = true, TFile* file3 = 0, std::string label3 = "", TFile* file4 = 0, std::string label4 = "");
TCanvas* compareNormalized(TH1F* h1, std::string label1, TH1F* h2, std::string label2, int rebin = 1, bool norm = true, TH1F* h3 = 0, std::string label3 = "", TH1F* h4 = 0, std::string label4 = "");

TH1F* cumulate (TH1F* in, bool increasing);
TGraphErrors* eff_rej (TH1F* signal, TH1F* background, bool normalize, bool increasing, bool print = false, float sf = 1.);
TGraph* s_over_rootb (TH1F* signal, TH1F* background, bool increasing, bool do_s_over_b = false, bool print = false, float sf = 1.);

TLegend* init_legend(float x1=0.5, float y1=0.5, float x2=0.88, float y2=0.8);
TLegend* legendize(TCanvas* c, const TString& opt = 'l', const TString& label1 = "", const TString& label2 = "", const TString& label3 = "");
float err_mult(float A, float B, float errA, float errB, float C);
float err_binomial(float A, float B, float errA, float errB);
TGraphErrors* makeMETGraph(TH1F* hdata, TH1F* hmc, float xoffset = 0.0, bool errOnData = false);

char* pm;         
char* delim;      
char* delimstart; 
char* delimend;   
char* ee;         
char* mm;         
char* em;         
char* e;         
char* m;         

int   startwidth;
int   width1;
int   width2;
int   linelength;

#endif
