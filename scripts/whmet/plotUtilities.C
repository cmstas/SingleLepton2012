#include "plotUtilities.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include "TCanvas.h"
#include "TObject.h"
#include "TLegend.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TPad.h"
#include "TCut.h"
#include "TProfile.h"
#include "THStack.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMath.h"
#include <sstream>
#include <iomanip>
#include "TF1.h"

//______________________________________________________________________________
float histError( TH1F* h , int minbin , int maxbin ){

  float err2 = 0;
  for( int i = minbin ; i <= maxbin ; ++i ){
    err2 += pow( h->GetBinError(i) , 2 );
  }

  return sqrt(err2);
}

//______________________________________________________________________________
void printLine( int latex ){

  if( latex == 1 ){
    cout << "\\hline" << endl;
  }
  else if( latex == 2 ){
    cout << "*\\hline" << endl;
  }
  else{
    for( int i = 0 ; i < linelength ; ++i ) cout << "-";
    cout << endl;
  }
}

//______________________________________________________________________________
void printHeader(bool dilep){

  if (dilep) {
    cout << delimstart << setw(startwidth) << "Sample"    << setw(width2)
	 << delim      << setw(width1) << ee          << setw(width2)
	 << delim      << setw(width1) << mm          << setw(width2)
	 << delim      << setw(width1) << em          << setw(width2)
	 << delim      << setw(width1) << "total"     << setw(width2) 
	 << delimend   << endl;
  } else {
    cout << delimstart << setw(startwidth) << "Sample"    << setw(width2)
	 << delim      << setw(width1) << e          << setw(width2)
	 << delim      << setw(width1) << m          << setw(width2)
	 << delim      << setw(width1) << "total"     << setw(width2) 
	 << delimend   << endl;
  }
}


//______________________________________________________________________________
void print( TH1F* h , string label , bool dilep, bool correlatedError ){

  stringstream se;
  stringstream sm;
  stringstream see;
  stringstream smm;
  stringstream sem;
  stringstream stot;

  if (dilep) {
    if( label == "Data" ){
      see  << Form( "%.0f" , h->GetBinContent(1) );
      smm  << Form( "%.0f" , h->GetBinContent(2) );
      sem  << Form( "%.0f" , h->GetBinContent(3) );
      stot << Form( "%.0f" , h->GetBinContent(4)       );
    } else if( label == "Data/Pred" ){
      see  << Form( "%.2f" , h->GetBinContent(1) ) << pm << Form( "%.2f" , h->GetBinError(1) );
      smm  << Form( "%.2f" , h->GetBinContent(2) ) << pm << Form( "%.2f" , h->GetBinError(2) );
      sem  << Form( "%.2f" , h->GetBinContent(3) ) << pm << Form( "%.2f" , h->GetBinError(3) );
      stot << Form( "%.2f" , h->GetBinContent(4) ) << pm << Form( "%.2f" , h->GetBinError(4) );
    }else{
      see  << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
      smm  << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
      sem  << Form( "%.1f" , h->GetBinContent(3) ) << pm << Form( "%.1f" , h->GetBinError(3) );
    
      float error = 0;
      if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2) + h->GetBinError(3);
      else                  error = h->GetBinError(4);
    
      stot << Form( "%.1f" , h->GetBinContent(4)       ) << pm << Form( "%.1f" , error  );
    }

    cout << delimstart << setw(startwidth) << label      << setw(width2)
	 << delim      << setw(width1) << see.str()  << setw(width2)
	 << delim      << setw(width1) << smm.str()  << setw(width2)
	 << delim      << setw(width1) << sem.str()  << setw(width2)
	 << delim      << setw(width1) << stot.str() << setw(width2)
	 << delimend   << endl;
  }

  else {
    if( label == "Data" ){
      se  << Form( "%.0f" , h->GetBinContent(1) );
      sm  << Form( "%.0f" , h->GetBinContent(2) );
      stot << Form( "%.0f" , h->GetBinContent(3)       );
    } else if( label == "Data/Pred" ){
      se  << Form( "%.2f" , h->GetBinContent(1) ) << pm << Form( "%.2f" , h->GetBinError(1) );
      sm  << Form( "%.2f" , h->GetBinContent(2) ) << pm << Form( "%.2f" , h->GetBinError(2) );
      stot << Form( "%.2f" , h->GetBinContent(3) ) << pm << Form( "%.2f" , h->GetBinError(3) );
    }else{
      se  << Form( "%.2f" , h->GetBinContent(1) ) << pm << Form( "%.2f" , h->GetBinError(1) );
      sm  << Form( "%.2f" , h->GetBinContent(2) ) << pm << Form( "%.2f" , h->GetBinError(2) );
    
      float error = 0;
      if( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2);
      else                  error = h->GetBinError(3);
    
      stot << Form( "%.2f" , h->GetBinContent(3)       ) << pm << Form( "%.2f" , error  );
    }

    cout << delimstart << setw(startwidth) << label      << setw(width2)
	 << delim      << setw(width1) << se.str()  << setw(width2)
	 << delim      << setw(width1) << sm.str()  << setw(width2)
	 << delim      << setw(width1) << stot.str() << setw(width2)
	 << delimend   << endl;
  }
  
}

//______________________________________________________________________________
void printYields( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* dir , bool doData, int latex, bool doWeighted ){

  initSymbols( latex );
  bool dilep = false;
  if (TString(dir).Contains("cr3") || TString(dir).Contains("cr4")) dilep = true;

  // h_events:
  //  bin 1: unweighted event yield
  //  bin 2: weighted event yield
  TString histname = "h_events";
  if (!TString(dir).Length() == 0) histname = Form("%s/",dir) + histname;
  TString histname_e = histname + "_e";
  TString histname_m = histname + "_m";
  TString histname_ee = histname + "_ee";
  TString histname_mm = histname + "_mm";
  TString histname_em = histname + "_em";

  printLine(latex);
  printHeader(dilep);
  printLine(latex);

  TH1F* hyield = 0;
  if (dilep) hyield = new TH1F("hyield","yield",4,0,4);
  else hyield = new TH1F("hyield","yield",3,0,3);
  hyield->Sumw2();
  TH1F* hmctot = 0;

  //----------------------
  // print SM MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < mcfiles.size() ; imc++){

    bool correlatedError = false;

    if( TString(labels[imc]).Contains("TChi") || TString(labels[imc]).Contains("Wino") ) continue;

    if (dilep) {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_ee = (TH1F*)mcfiles[imc]->Get(histname_ee.Data());
      TH1F* mchist_mm = (TH1F*)mcfiles[imc]->Get(histname_mm.Data());
      TH1F* mchist_em = (TH1F*)mcfiles[imc]->Get(histname_em.Data());

      fillYieldHist(mchist,hyield,4,doWeighted);
      fillYieldHist(mchist_ee,hyield,1,doWeighted);
      fillYieldHist(mchist_mm,hyield,2,doWeighted);
      fillYieldHist(mchist_em,hyield,3,doWeighted);
    }
    else {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_e = (TH1F*)mcfiles[imc]->Get(histname_e.Data());
      TH1F* mchist_m = (TH1F*)mcfiles[imc]->Get(histname_m.Data());

      fillYieldHist(mchist,hyield,3,doWeighted);
      fillYieldHist(mchist_e,hyield,1,doWeighted);
      fillYieldHist(mchist_m,hyield,2,doWeighted);
    }

    if( imc == 0 ) hmctot = (TH1F*) hyield->Clone();
    else           hmctot->Add(hyield);
    
    print( hyield , getTableName(labels[imc],latex) , dilep, correlatedError );

    //hyield->Reset();
  }

  printLine(latex);

  //-------------------------------
  // print sum of SM MC samples
  //-------------------------------

  print( hmctot , "Total SM Pred", dilep );

  printLine(latex);
 
  if (doData) {

    if (dilep) {
      TH1F* datahist = (TH1F*)datafile->Get(histname.Data());
      TH1F* datahist_ee = (TH1F*)datafile->Get(histname_ee.Data());
      TH1F* datahist_mm = (TH1F*)datafile->Get(histname_mm.Data());
      TH1F* datahist_em = (TH1F*)datafile->Get(histname_em.Data());

      fillYieldHist(datahist,hyield,4,doWeighted);
      fillYieldHist(datahist_ee,hyield,1,doWeighted);
      fillYieldHist(datahist_mm,hyield,2,doWeighted);
      fillYieldHist(datahist_em,hyield,3,doWeighted);
    }
    else {
      TH1F* datahist = (TH1F*)datafile->Get(histname.Data());
      TH1F* datahist_e = (TH1F*)datafile->Get(histname_e.Data());
      TH1F* datahist_m = (TH1F*)datafile->Get(histname_m.Data());

      fillYieldHist(datahist,hyield,3,doWeighted);
      fillYieldHist(datahist_e,hyield,1,doWeighted);
      fillYieldHist(datahist_m,hyield,2,doWeighted);
    }

    print( hyield , "Data", dilep );
    
    printLine(latex);

    TH1F* hratio = (TH1F*)hyield->Clone("hratio");
    hratio->Divide(hmctot);

    print( hratio , "Data/Pred", dilep );
    
    printLine(latex);
  }

  //----------------------
  // print SUSY MC samples
  //----------------------

  for(unsigned int imc = 0 ; imc < mcfiles.size() ; imc++){

    if( !TString(labels[imc]).Contains("TChi") && !TString(labels[imc]).Contains("Wino") )   continue;

    if (dilep) {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_ee = (TH1F*)mcfiles[imc]->Get(histname_ee.Data());
      TH1F* mchist_mm = (TH1F*)mcfiles[imc]->Get(histname_mm.Data());
      TH1F* mchist_em = (TH1F*)mcfiles[imc]->Get(histname_em.Data());

      fillYieldHist(mchist,hyield,4,doWeighted);
      fillYieldHist(mchist_ee,hyield,1,doWeighted);
      fillYieldHist(mchist_mm,hyield,2,doWeighted);
      fillYieldHist(mchist_em,hyield,3,doWeighted);

    }
    else {
      TH1F* mchist = (TH1F*)mcfiles[imc]->Get(histname.Data());
      TH1F* mchist_e = (TH1F*)mcfiles[imc]->Get(histname_e.Data());
      TH1F* mchist_m = (TH1F*)mcfiles[imc]->Get(histname_m.Data());

      fillYieldHist(mchist,hyield,3,doWeighted);
      fillYieldHist(mchist_e,hyield,1,doWeighted);
      fillYieldHist(mchist_m,hyield,2,doWeighted);
    }

    print( hyield , getTableName(labels[imc],1), dilep );

  }

  printLine(latex);

  delete hyield;
  delete hmctot;
}

//______________________________________________________________________________
void fillYieldHist( TH1F* hin, TH1F* hyield, int bin, bool doWeighted ) {

  int events_bin = 2;
  if (!doWeighted) events_bin = 1;
  Double_t err = 0.;
  if (hin) {
    // float nmc = hin->IntegralAndError(0,-1,err);
    // hyield->SetBinContent(bin,nmc);
    // hyield->SetBinError(bin,err);
    // these assume h_events with bin 2 having the weighted event yields
    hyield->SetBinContent(bin,hin->GetBinContent(events_bin));
    hyield->SetBinError(bin,hin->GetBinError(events_bin));
  } else {
    hyield->SetBinContent(bin,0.);
    hyield->SetBinError(bin,0.);
  }

}

//______________________________________________________________________________
void printCutflow( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , vector<char*> dirs , bool doData, bool splitMC, int latex, bool doWeighted ){

  initSymbols( latex );

  std::vector<TString> bg_labels;
  std::vector<TFile*> bg_files;
  std::vector<TString> sig_labels;
  std::vector<TFile*> sig_files;

  assert(mcfiles.size() == labels.size());

  // loop through mc files and separate into signals and bgs
  for(unsigned int imc = 0 ; imc < mcfiles.size() ; ++imc) {
    TString label = TString(labels[imc]);
    if( label.Contains("TChi") || label.Contains("Wino") ) {
      sig_labels.push_back(label);
      sig_files.push_back(mcfiles[imc]);
    } else {
      bg_labels.push_back(label);
      bg_files.push_back(mcfiles[imc]);
    }
  }

  // h_events:
  //  bin 1: unweighted event yield
  //  bin 2: weighted event yield
  TString histname = "h_events";

  int events_bin = 2;
  if (!doWeighted) events_bin = 1;

  printLine(latex);
  // header info
  cout << delimstart << setw(startwidth) << "Selection"    << setw(width2);
  if (splitMC) {
    for(unsigned int ibg = 0 ; ibg < bg_labels.size() ; ++ibg){
      cout << delim      << setw(width1) << getTableName(bg_labels[ibg],latex)    << setw(width2);
    }
  }
  cout << delim      << setw(width1) << "Total Pred"     << setw(width2);
  if (doData) cout << delim      << setw(width1) << "Data"          << setw(width2)
		   << delim      << setw(width1) << "Data/Pred"          << setw(width2);
  // print signals here
  for(unsigned int isig = 0 ; isig < sig_labels.size() ; ++isig){
    cout << delim      << setw(width1) << getTableName(sig_labels[isig],1)    << setw(width2);
  }
  cout << delimend   << endl;
  printLine(latex);

  // loop over dirs, print counts for each
  for (unsigned int idir = 0; idir < dirs.size(); ++idir) {
    TString histname_dir = histname;
    if (!TString(dirs[idir]).Length() == 0) histname_dir = Form("%s/",dirs[idir]) + histname;
    float nmc_tot = 0.;
    float errmc_tot = 0.;

    cout << delimstart << setw(startwidth) << dirs[idir]  << setw(width2);
    // loop over backgrounds
    for(unsigned int ibg = 0 ; ibg < bg_labels.size() ; ++ibg){
      TH1F* mchist = (TH1F*)bg_files[ibg]->Get(histname_dir.Data());
      float nmc = 0.;
      float errmc = 0.;
      if (mchist) {
	float nmc = mchist->GetBinContent(events_bin);
	float errmc = mchist->GetBinError(events_bin);
	nmc_tot += nmc;
	errmc_tot += pow(errmc,2);
      }

      if (splitMC) {
	cout << delim      << setw(width1) << Form( "%.1f" , nmc ) << pm << Form( "%.1f" , errmc )  << setw(width2);
      }
    } // loop over backgrounds
    errmc_tot = sqrt(errmc_tot);

    cout << delim      << setw(width1) << Form( "%.1f" , nmc_tot ) << pm << Form( "%.1f" , errmc_tot )      << setw(width2);

    if (doData) {
      TH1F* datahist = (TH1F*)datafile->Get(histname_dir.Data());
      float ndata = 0.;
      float errdata = 0.;
      float dataovermc = 0.;
      float errdataovermc = 0.;

      if (datahist) {
	ndata = datahist->GetBinContent(events_bin);
	errdata = datahist->GetBinError(events_bin);

	if (nmc_tot > 0.) {
	  dataovermc = ndata/nmc_tot;
	  if (ndata > 0.) {
	    errdataovermc = err_mult(ndata,nmc_tot,errdata,errmc_tot,dataovermc);
	  }
	}
      }

      cout << delim      << setw(width1) << Form("%d",int(ndata))  << setw(width2)
	   << delim      << setw(width1) << Form("%.2f",dataovermc) << pm << Form("%.2f", errdataovermc )  << setw(width2);
    }
    // print signals here
    for(unsigned int isig = 0 ; isig < sig_labels.size() ; ++isig){
      TH1F* mchist = (TH1F*)sig_files[isig]->Get(histname_dir.Data());
      float nmc = mchist->GetBinContent(events_bin);
      float errmc = mchist->GetBinError(events_bin);
      cout << delim      << setw(width1) << Form( "%.1f" , nmc ) << pm << Form( "%.1f" , errmc )  << setw(width2);
    }
    cout << delimend   << endl;

  } // loop over dirs

  printLine(latex);

}


//______________________________________________________________________________
void printSigRegions( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , vector<char*> dirs , bool doData, int latex, bool doWeighted ){

  initSymbols( latex );

  std::vector<TString> bg_labels;
  std::vector<TFile*> bg_files;
  std::vector<TString> sig_labels;
  std::vector<TFile*> sig_files;

  assert(mcfiles.size() == labels.size());

  // loop through mc files and separate into signals and bgs
  for(unsigned int imc = 0 ; imc < mcfiles.size() ; ++imc) {
    TString label = TString(labels[imc]);
    if( label.Contains("TChi") || label.Contains("Wino") ) {
      sig_labels.push_back(label);
      sig_files.push_back(mcfiles[imc]);
    } else {
      bg_labels.push_back(label);
      bg_files.push_back(mcfiles[imc]);
    }
  }

  // h_events:
  //  bin 1: unweighted event yield
  //  bin 2: weighted event yield
  TString histname = "h_events";

  int events_bin = 2;
  if (!doWeighted) events_bin = 1;

  printLine(latex);
  // header info
  cout << delimstart << setw(startwidth) << "Sample"    << setw(width2);
  for (unsigned int idir = 0; idir < dirs.size(); ++idir) {
    cout << delim      << setw(width1) << dirs[idir]    << setw(width2);
  }
  cout << delimend   << endl;
  printLine(latex);

  std::vector<float> nmc_tot(dirs.size(),0.);
  std::vector<float> errmc_tot(dirs.size(),0.);

  // loop over backgrounds
  for(unsigned int ibg = 0 ; ibg < bg_labels.size() ; ++ibg){
    cout << delimstart << setw(startwidth) << getTableName(bg_labels[ibg],latex)  << setw(width2);
    // loop over dirs, print counts for each
    for (unsigned int idir = 0; idir < dirs.size(); ++idir) {
      TString histname_dir = histname;
      if (!TString(dirs[idir]).Length() == 0) histname_dir = Form("%s/",dirs[idir]) + histname;

      TH1F* mchist = (TH1F*)bg_files[ibg]->Get(histname_dir.Data());
      float nmc = 0.;
      float errmc = 0.;
      if (mchist) {
	nmc = mchist->GetBinContent(events_bin);
	float errmc_stat = mchist->GetBinError(events_bin);
	errmc = sqrt(pow(getSystError(bg_labels[ibg]) * nmc, 2) + pow(errmc_stat,2));
	nmc_tot[idir] += nmc;
	errmc_tot[idir] += pow(errmc,2);
      }

      stringstream out;
      out << Form( "%.1f" , nmc ) << pm << Form( "%.1f" , errmc );
      cout << delim      << setw(width1) << out.str()  << setw(width2);

    } // loop over dirs
    cout << delimend   << endl;

  } // loop over bgs
  printLine(latex);

  // total predictions for each dir
  cout << delimstart << setw(startwidth) << "Total SM Pred"  << setw(width2);
  for (unsigned int idir = 0; idir < dirs.size(); ++idir) {
    errmc_tot[idir] = sqrt(errmc_tot[idir]);
    stringstream out;
    out << Form( "%.1f" , nmc_tot[idir] ) << pm << Form( "%.1f" , errmc_tot[idir] );
    cout << delim      << setw(width1) << out.str()  << setw(width2);
  }
  cout << delimend   << endl;
  printLine(latex);

  // print data (or placeholder if blinding)
  cout << delimstart << setw(startwidth) << "Data"  << setw(width2);
  for (unsigned int idir = 0; idir < dirs.size(); ++idir) {
    stringstream out;
    if (!doData) {
      out << "-";
      cout << delim      << setw(width1) << out.str()  << setw(width2);
      continue;
    }
    TString histname_dir = histname;
    if (!TString(dirs[idir]).Length() == 0) histname_dir = Form("%s/",dirs[idir]) + histname;
    TH1F* datahist = (TH1F*)datafile->Get(histname_dir.Data());
    float ndata = 0.;
    float errdata = 0.;
    if (datahist) {
      ndata = datahist->GetBinContent(events_bin);
      errdata = datahist->GetBinError(events_bin);
    } 
    out << Form("%d",int(ndata));
    cout << delim      << setw(width1) << out.str()  << setw(width2);
  }
  cout << delimend   << endl;
  printLine(latex);

  // print signals here
  for(unsigned int isig = 0 ; isig < sig_labels.size() ; ++isig){
    cout << delimstart << setw(startwidth) << getTableName(sig_labels[isig],latex)  << setw(width2);
    for (unsigned int idir = 0; idir < dirs.size(); ++idir) {
      TString histname_dir = histname;
      if (!TString(dirs[idir]).Length() == 0) histname_dir = Form("%s/",dirs[idir]) + histname;
      TH1F* mchist = (TH1F*)sig_files[isig]->Get(histname_dir.Data());
      float nmc = 0.;
      float errmc = 0.;
      if (mchist) {
	nmc = mchist->GetBinContent(events_bin);
	float errmc_stat = mchist->GetBinError(events_bin);
	errmc = sqrt(pow(getSystError(sig_labels[isig]) * nmc, 2) + pow(errmc_stat,2));
      }
      stringstream out;
      out << Form( "%.1f" , nmc ) << pm << Form( "%.1f" , errmc );
      cout << delim      << setw(width1) << out.str()  << setw(width2);

    } // loop over dirs
    cout << delimend   << endl;

  } // loop over sigs

  printLine(latex);

}

//______________________________________________________________________________
void printSigRegionsExcl( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , char* dir , bool doData, int latex, bool doSystErr ){

  initSymbols( latex );

  std::vector<int> metbins;
  //  metbins.push_back(50);
  metbins.push_back(100);
  metbins.push_back(125);
  metbins.push_back(150);
  metbins.push_back(175);

  std::vector<TString> bg_labels;
  std::vector<TFile*> bg_files;
  std::vector<TString> sig_labels;
  std::vector<TFile*> sig_files;

  assert(mcfiles.size() == labels.size());

  // loop through mc files and separate into signals and bgs
  for(unsigned int imc = 0 ; imc < mcfiles.size() ; ++imc) {
    TString label = TString(labels[imc]);
    if( label.Contains("TChi") || label.Contains("Wino") ) {
      sig_labels.push_back(label);
      sig_files.push_back(mcfiles[imc]);
    } else {
      bg_labels.push_back(label);
      bg_files.push_back(mcfiles[imc]);
    }
  }

  TString histname_dir = TString(dir) + "/h_met";

  printLine(latex);
  // header info
  cout << delimstart << setw(startwidth) << "Sample"    << setw(width2);
  for (unsigned int ibin = 0; ibin < metbins.size(); ++ibin) {
    cout << delim      << setw(width1) << metbins.at(ibin)    << setw(width2);
  }
  cout << delimend   << endl;
  printLine(latex);

  std::vector<float> nmc_tot(metbins.size(),0.);
  std::vector<float> errmc_tot(metbins.size(),0.);

  // loop over backgrounds
  for(unsigned int ibg = 0 ; ibg < bg_labels.size() ; ++ibg){
    cout << delimstart << setw(startwidth) << getTableName(bg_labels[ibg],latex)  << setw(width2);
    // loop over dirs, print counts for each
    for (unsigned int ibin = 0; ibin < metbins.size(); ++ibin) {
      TH1F* mchist = (TH1F*)bg_files[ibg]->Get(histname_dir.Data());
      float nmc = 0.;
      float errmc = 0.;
      if (mchist) {
	Double_t errmc_stat;
	int maxbin = -1;
	if (ibin < metbins.size()-1) maxbin = metbins.at(ibin+1);
	nmc = mchist->IntegralAndError(metbins.at(ibin)+1,maxbin,errmc_stat);
	if (doSystErr) errmc = sqrt(pow(getSystError(bg_labels[ibg]) * nmc, 2) + pow(errmc_stat,2));
	else errmc = errmc_stat;
	nmc_tot[ibin] += nmc;
	errmc_tot[ibin] += pow(errmc,2);
      }

      stringstream out;
      out << Form( "%.1f" , nmc ) << pm << Form( "%.1f" , errmc );
      cout << delim      << setw(width1) << out.str()  << setw(width2);

    } // loop over dirs
    cout << delimend   << endl;

  } // loop over bgs
  printLine(latex);

  // total predictions for each dir
  cout << delimstart << setw(startwidth) << "Total SM Pred"  << setw(width2);
  for (unsigned int ibin = 0; ibin < metbins.size(); ++ibin) {
    errmc_tot[ibin] = sqrt(errmc_tot[ibin]);
    stringstream out;
    out << Form( "%.1f" , nmc_tot[ibin] ) << pm << Form( "%.1f" , errmc_tot[ibin] );
    cout << delim      << setw(width1) << out.str()  << setw(width2);
  }
  cout << delimend   << endl;
  printLine(latex);

  // print data (or placeholder if blinding)
  cout << delimstart << setw(startwidth) << "Data"  << setw(width2);
  for (unsigned int ibin = 0; ibin < metbins.size(); ++ibin) {
    stringstream out;
    if (!doData) {
      out << "-";
      cout << delim      << setw(width1) << out.str()  << setw(width2);
      continue;
    }
    float ndata = 0.;
    float errdata = 0.;
    TH1F* datahist = (TH1F*)datafile->Get(histname_dir.Data());
    if (datahist) {
      Double_t errdata;
      int maxbin = -1;
      if (ibin < metbins.size()-1) maxbin = metbins.at(ibin+1);
      ndata = datahist->IntegralAndError(metbins.at(ibin)+1,maxbin,errdata);
    } 
    out << Form("%d",int(ndata));
    cout << delim      << setw(width1) << out.str()  << setw(width2);
  }
  cout << delimend   << endl;
  printLine(latex);

  // print signals here
  for(unsigned int isig = 0 ; isig < sig_labels.size() ; ++isig){
    cout << delimstart << setw(startwidth) << getTableName(sig_labels[isig],latex)  << setw(width2);
    for (unsigned int ibin = 0; ibin < metbins.size(); ++ibin) {
      TH1F* mchist = (TH1F*)sig_files[isig]->Get(histname_dir.Data());
      float nmc = 0.;
      float errmc = 0.;
      if (mchist) {
	Double_t errmc_stat;
	int maxbin = -1;
	if (ibin < metbins.size()-1) maxbin = metbins.at(ibin+1);
	nmc = mchist->IntegralAndError(metbins.at(ibin)+1,maxbin,errmc_stat);
	if (doSystErr) errmc = sqrt(pow(getSystError(sig_labels[isig]) * nmc, 2) + pow(errmc_stat,2));
	else errmc = errmc_stat;
      }
      stringstream out;
      out << Form( "%.1f" , nmc ) << pm << Form( "%.1f" , errmc );
      cout << delim      << setw(width1) << out.str()  << setw(width2);

    } // loop over dirs
    cout << delimend   << endl;

  } // loop over sigs

  printLine(latex);

}

//______________________________________________________________________________
float getSystError( const TString sample ){

  if (sample.Contains("ttbar 2l")) return 0.4;
  else if (sample.Contains("ttbar 1l")) return 0.5;
  else if (sample.Contains("wlight")) return 0.4;
  //  else if (sample.Contains("wbb")) return 0.5;
  else if (sample.Contains("wbb")) return 0.43;
  else if (sample.Contains("wzbb")) return 0.36;
  else if (sample.Contains("single top 1l")) return 0.5;
  else if (sample.Contains("single top 2l")) return 0.4;
  else if (sample.Contains("rare")) return 0.5;
  else if (sample.Contains("whbb")) return 0.5;
  else if (sample.Contains("dilep top")) return 0.4;
  else if (sample.Contains("lep plus b")) return 0.5;
  //  else if (sample.Contains("top 1l")) return 0.5;
  else if (sample.Contains("top 1l")) return 0.57;
  else if (sample.Contains("TChi")) return 0.12;
  else if (sample.Contains("Wino")) return 0.12;
  //  else if (sample.Contains("TChi")) return 0.1;

  std::cout << "ERROR: getSystError: didn't recognize sample: " << sample << std::endl;
  return 0.;

}

//______________________________________________________________________________
std::string getTableName( const TString sample, int latex ){

  if (latex == 0) {
    return std::string(sample.Data());
  }

  if (sample.Contains("ttbar 2l")) return "$t\\bar{t} \\rightarrow \\ell\\ell$";
  else if (sample.Contains("ttbar 1l")) return "$t\\bar{t} \\rightarrow \\ell+$ jets";
  else if (sample.Contains("wlight")) return "W + light";
  else if (sample.Contains("wzbb")) return "WZ $\\rightarrow \\ell\\nu b\\bar{b}$";
  else if (sample.Contains("wbb")) return "W + $b\\bar{b}$";
  else if (sample.Contains("single top 1l")) return "single top $1\\ell$";
  else if (sample.Contains("single top 2l")) return "$tW \\rightarrow \\ell\\ell$";
  else if (sample.Contains("rare")) return "Rare";
  else if (sample.Contains("whbb")) return "WH $\\rightarrow \\ell\\nu b\\bar{b}$";
  else if (sample.Contains("dilep top")) return "Dilepton top";
  else if (sample.Contains("lep plus b")) return "Single lepton + b";
  else if (sample.Contains("top 1l")) return "Single lepton top";
  else if (sample.Contains("_")) {
    TString copy(sample);
    return std::string(copy.ReplaceAll("_","\\_").Data());
  } 

  // no match: return input
  return std::string(sample.Data());

}

#include <TList.h>
#include <TIterator.h>

//______________________________________________________________________________
void deleteHistos() {
   // Delete all existing histograms in memory
   TObject* obj;
   TList* list = gDirectory->GetList() ;
   TIterator* iter = list->MakeIterator();
   while ((obj=iter->Next())) {
     if (obj->IsA()->InheritsFrom(TH1::Class()) ||
         obj->IsA()->InheritsFrom(TH2::Class()) ) {delete obj;}
   }
}

//______________________________________________________________________________
void initSymbols( int latex ){

  //-------------------------------------------------------
  // table format
  //-------------------------------------------------------

  startwidth  = 35;
  width1      = 22;
  width2      = 5;
  linelength  = (width1+width2)*5+1;

  //-------------------------------------------------------
  // symbols
  //-------------------------------------------------------
  
  if( latex == 1 ){
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "";
    delimend   = "\\\\";
    ee         = "ee";
    mm         = "$\\mu\\mu$";
    em         = "e$\\mu$";
    e         = "e";
    m         = "$\\mu$";
  } else if (latex == 2) {
    pm         = " $\\pm$ ";
    delim      = "&";
    delimstart = "*";
    delimend   = "\\\\";
    ee         = "ee";
    mm         = "$\\mu\\mu$";
    em         = "e$\\mu$";
    e         = "e";
    m         = "$\\mu$";
  }else{
    pm         = " +/- ";
    delim      = "|";
    delimstart = "|";
    delimend   = "|";
    ee         = "ee";
    mm         = "mm";
    em         = "em";
    e         = "e";
    m         = "m";
  }

}


//______________________________________________________________________________
int getColor(const TString sample) {

  // // original color scheme
  // if (sample.Contains("ttbar 2l")) return 7;
  // else if (sample.Contains("ttbar 1l")) return 4;
  // else if (sample.Contains("wlight")) return 2;
  // else if (sample.Contains("wbb")) return 5;
  // else if (sample.Contains("single top 1l")) return 8;
  // else if (sample.Contains("single top 2l")) return 9;
  // else if (sample.Contains("rare")) return 15;
  // else if (sample.Contains("whbb")) return 12;
  // else if (sample.Contains("dilep top")) return 7;
  // else if (sample.Contains("lep plus b")) return 4;
  // else if (sample.Contains("top 1l")) return 4;
  // else if (sample.Contains("_200_1")) return 2;
  // else if (sample.Contains("_250_1")) return 4;
  // else if (sample.Contains("_300_1")) return 3;

  // "standard" color scheme from top group
  if (sample.Contains("ttbar 2l")) return kRed-4; 
  else if (sample.Contains("ttbar 1l")) return kRed-3;
  else if (sample.Contains("wlight")) return kGreen-2;
  else if (sample.Contains("wbb")) return kGreen-3;
  else if (sample.Contains("single top 1l")) return kMagenta+2;
  else if (sample.Contains("single top 2l")) return kMagenta;
  else if (sample.Contains("rare")) return 15;
  else if (sample.Contains("wzbb")) return kBlue+1;
  else if (sample.Contains("whbb")) return 12;
  else if (sample.Contains("dilep top")) return kRed-4;
  else if (sample.Contains("lep plus b")) return kRed-3;
  else if (sample.Contains("top 1l")) return kRed-3;
  else if (sample.Contains("_150_1")) return kBlue-4;
  else if (sample.Contains("_200_1")) return kBlue-6;
  else if (sample.Contains("_250_1")) return kMagenta;
  else if (sample.Contains("_300_1")) return kCyan+1;
  else if (sample.Contains("_350_1")) return kCyan-6;

  std::cout << "ERROR: getColor: didn't recognize sample: " << sample << std::endl;
  return 1;
}

//______________________________________________________________________________
TLegend *getLegend( vector<char*> labels , bool overlayData, float x1, float y1, float x2, float y2, float signorm, bool errband){

  TLegend *leg = new TLegend(x1,y1,x2,y2);

  TH1F*    datahist = new TH1F("datahist","datahist",1,0,1);
  datahist->Sumw2();

  if( overlayData ) leg->AddEntry(datahist,"Data");

  const unsigned int nmc = labels.size();
  TH1F*    mchist[nmc];

  //-----------------
  // SM samples
  //-----------------

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    TString mclabel(labels.at(imc));

    if( mclabel.Contains("TChi") || mclabel.Contains("Wino") )   continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    // formatting
    mchist[imc]->SetFillColor( getColor(labels.at(imc)) );

    if( mclabel == "tt" ) mclabel = "t#bar{t}";
    else if( mclabel == "ttbar" ) mclabel = "t#bar{t}";
    else if( mclabel == "ttbar 1l" ) mclabel = "t#bar{t} #rightarrow 1l";
    else if( mclabel == "ttbar 2l" ) mclabel = "t#bar{t} #rightarrow 2l";
    else if( mclabel == "ttll" ) mclabel = "t#bar{t} #rightarrow ll";
    else if( mclabel == "tttau" ) mclabel = "t#bar{t} #rightarrow l#tau/#tau#tau";
    else if( mclabel == "ttfake" ) mclabel = "t#bar{t} #rightarrow fake";
    else if( mclabel == "t" ) mclabel = "single top";
    else if( mclabel == "single_top" ) mclabel = "single top";
    else if( mclabel == "wjets" ) mclabel = "W+jets";
    else if( mclabel == "wlight" ) mclabel = "W+light jets";
    else if( mclabel == "wbb" ) mclabel = "W+b#bar{b}";
    else if( mclabel == "zjets" ) mclabel = "Z+jets";
    else if( mclabel == "ww" ) mclabel = "WW";
    else if( mclabel == "wz" ) mclabel = "WZ";
    else if( mclabel == "zz" ) mclabel = "ZZ";
    else if( mclabel == "wzbb" ) mclabel = "WZ #rightarrow l#nub#bar{b}";
    else if( mclabel == "whbb" ) mclabel = "WH #rightarrow l#nub#bar{b}";
    else if( mclabel == "dilep top" ) mclabel = "2l top";
    else if( mclabel == "top 1l" ) mclabel = "1l top";
    else if( mclabel == "rare" ) mclabel = "Rare";

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],mclabel,"f");
    
  }

  if (errband) {
    TGraphErrors* bgerrgraph = new TGraphErrors();
    bgerrgraph->SetLineColor(0);
    bgerrgraph->SetFillColor(kGray+2);
    bgerrgraph->SetFillStyle(3005);
    leg->AddEntry(bgerrgraph,"Total Uncertainty","f");
  }

  //-----------------
  // sig samples
  //-----------------

  int nsigmc = 0;

  for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  //  for( int imc = nmc - 1 ; imc >= 0 ; imc-- ){

    TString mclabel(labels.at(imc));

    if( !mclabel.Contains("TChi") && !mclabel.Contains("Wino") )   continue;

    mchist[imc] = new TH1F(Form("mc_%i",imc),Form("mc_%i",imc),1,0,1);

    // formatting
    mchist[imc]->SetFillColor( 0 );
    mchist[imc]->SetLineColor( getColor(labels.at(imc)) );
    mchist[imc]->SetLineWidth(2);
    //    mchist[imc]->SetLineStyle(2);
    ++nsigmc;

    if (mclabel.Contains("Wino_")) {
      mclabel.ReplaceAll("_","/");
      mclabel.ReplaceAll("Wino/","#tilde{#chi}_{1}^{#pm}#tilde{#chi}_{2}^{0} #rightarrow (W#tilde{#chi}_{1}^{0})(H#tilde{#chi}_{1}^{0}) (");
      mclabel += ")";
      if (signorm > 0.) mclabel += Form(" x%d",int(signorm));
    }

    //leg->AddEntry(mchist[imc],labels.at(imc),"f");
    leg->AddEntry(mchist[imc],mclabel,"f");
    
  }

  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.026);
  
  return leg;

}

//______________________________________________________________________________
//TH1F* compareDataMC( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , const char* flavor , const char* dir ,
TGraphErrors* compareDataMC( vector<TFile*> mcfiles , vector<char*> labels , TFile* datafile , const char* histname , 
			     const char* flavor , const char* dir ,
			     int nbins ,  float xmin , float xmax ,  
			     const char* xtitle , bool overlayData , bool residual , bool drawLegend , bool log , 
			     bool normalize , bool fit, float mcnorm, const char* scalesample,
			     bool stacksig, float signorm, bool errband, const char* cmslabel ){

  dataMCHists datamchists;
  TH1F* ratio = 0;
  TGraphErrors* gratio = 0;

  TString fullhistname = histname;
  TString tsdir = TString(dir);
  if (!tsdir.Length() == 0) fullhistname = Form("%s/",dir) + fullhistname;

  TString ytitle = "Events";

  TPad* fullpad = new TPad();
  TPad* plotpad = new TPad();
  TPad* respad  = new TPad();

  if( residual ){
    fullpad = new TPad("fullpad","fullpad",0,0,1,1);
    //    fullpad->SetRightMargin(0.05);
    plotpad->SetBottomMargin(0.05);
    fullpad->Draw();
    fullpad->cd();

    plotpad = new TPad("plotpad","plotpad",0,0,1,0.8);
    plotpad->SetRightMargin(0.05);
    plotpad->Draw();
    plotpad->cd();
    if( log ) plotpad->SetLogy();
  }
  else{
    gPad->SetRightMargin(0.05);
    //    gPad->SetBottomMargin(0.05);
    if( log ) gPad->SetLogy();
  }


  cout << "Plotting " << fullhistname << endl;

  assert( mcfiles.size() == labels.size() );
  const unsigned int nmc = mcfiles.size();

  THStack* mcstack = new THStack("mcstack","mcstack");
  TH1F*    mctothist = 0;
  vector<TH1F*> mchist;
  mchist.resize(nmc);
  TH1F* datahist = 0;
  if (overlayData) {
    TH1F* datahist_tmp = (TH1F*)datafile->Get(fullhistname.Data());
    if (datahist_tmp) {
      datahist = (TH1F*)datahist_tmp->Clone(fullhistname+TString("_data"));
    } else {
      std::cout << "WARNING: no hist with name " << fullhistname << " found for data" << std::endl;
      datahist = 0;
    }
  }
  float nmctot = 0.0;
  float errmctot = 0.0;

  bool drawingSignal = false;
  for( int imc = nmc-1 ; imc > -1 ; imc-- ){
    TH1F* mchist_tmp = (TH1F*)mcfiles[imc]->Get(fullhistname.Data());
    if (mchist_tmp) {
      mchist[imc] = (TH1F*)mchist_tmp->Clone(fullhistname+Form("_%s",labels[imc]));
      Double_t err;
      float nmc = mchist[imc]->IntegralAndError(0,-1,err);
      if (TString(labels[imc]).Contains("TChiwh") || TString(labels[imc]).Contains("Wino")) {
	// signal: don't add to bg total
	cout << "MC signal yield " << labels[imc] << " " << Form("%.2f",nmc) << " +/- " << Form("%.2f",err) << endl;
	drawingSignal = true;
      } else {
	nmctot += mchist[imc]->Integral();
	Double_t err;
	float nmc = mchist[imc]->IntegralAndError(0,-1,err);
	errmctot = sqrt(pow(errmctot,2) + pow(err,2));
	cout << "MC yield " << labels[imc] << " " << Form("%.2f",nmc) << " +/- " << Form("%.2f",err) << endl;
      }
    } else {
      std::cout << "WARNING: no hist with name " << fullhistname << " found for sample " << labels[imc] << std::endl;
      mchist[imc] = 0;
    }
  }

  cout << "MC bg total yield " << Form("%.2f",nmctot) << " +/- " << Form("%.2f",errmctot) << endl;

  float ndata = 0.;

  if (overlayData) {
    ndata = datahist->Integral();
    cout << "data yield " << ndata << endl;
  }

  //  TString scalesample = "all";
  TString tscalesample = TString(scalesample);
  if (tscalesample.Length() == 0) tscalesample = "all";

  float SF = ndata/nmctot;
  if( normalize ){
    cout << "Data, MC, SF " << ndata << ", " << nmctot << ", " << SF << endl;
    if (mcnorm < 0.) cout << "Scaling MC: " << tscalesample << " by " << SF << "!!!!!!" << endl;
    else cout << "Scaling MC: " << tscalesample << " by input factor of " << mcnorm << "!!!!!!" << endl;
  }

  int nmcsig = 0;
  std::vector<TH1F*> mcsighist;
  int rebinFactor = 1;
  float xbinsize = 1.;

  //for( unsigned int imc = 0 ; imc < nmc ; imc++ ){
  for( int imc = nmc-1 ; imc > -1 ; imc-- ){

    // skip samples without hists
    if (!mchist[imc]) continue;

    mchist[imc]->Sumw2();
    float xmin_original = mchist[imc]->GetXaxis()->GetXmin();
    float xmax_original = mchist[imc]->GetXaxis()->GetXmax();
    int nbins_original = mchist[imc]->GetNbinsX();
    float xbinsize_original = (xmax_original - xmin_original)/float(nbins_original);
    xbinsize = (xmax - xmin)/float(nbins);
    rebinFactor = int(xbinsize/xbinsize_original);
    if (rebinFactor != 1) mchist[imc]->Rebin(rebinFactor);
    mchist[imc]->SetAxisRange(xmin,xmax-0.1*xbinsize,"X");
    Double_t err;
    float lastbin_and_overflow =  mchist[imc]->IntegralAndError(nbins,-1,err);
    if (!TString(histname).Contains("h_bdt")) {
      mchist[imc]->SetBinContent(nbins,lastbin_and_overflow);
      mchist[imc]->SetBinError(nbins,err);
    }

    //    if( normalize && ( tscalesample.EqualTo("all") || tscalesample.EqualTo(labels[imc]) ) ) {
    if( TString( labels.at(imc) ).Contains("TChiwh") || TString( labels.at(imc) ).Contains("Wino") ){
      // signal
      mchist[imc]->SetFillColor( 0 );
      mchist[imc]->SetLineWidth(2);
      mchist[imc]->SetLineColor( getColor(labels.at(imc)) );
      //      mchist[imc]->SetLineStyle(2);
      ++nmcsig;
      mcsighist.push_back(mchist[imc]);
    }else{
      // bg
      if( normalize && ( tscalesample.EqualTo("all") || TString(labels[imc]).Contains(tscalesample) ) ) {
	if (mcnorm < 0.) mchist[imc]->Scale(SF);
	else mchist[imc]->Scale(mcnorm);
      }
      mchist[imc]->SetLineColor(kBlack);
      mchist[imc]->SetFillColor( getColor(labels.at(imc)) );
      mcstack->Add( mchist[imc] );

      //if( imc == 0 ) mctothist = (TH1F*) mchist[imc]->Clone();
      //    if( imc == nmc-1 ) mctothist = (TH1F*) mchist[imc]->Clone();
      if( !mctothist ) mctothist = (TH1F*) mchist[imc]->Clone();
      else               mctothist->Add(mchist[imc]);

    }

  }

  mctothist->SetLineColor(kBlack);
  TGraphErrors* bgerrgraph = new TGraphErrors(mctothist);
  // loop over bins and set all errors to 25% (for now)
  if (errband) {
    bgerrgraph->SetFillColor(kGray+2);
    bgerrgraph->SetFillStyle(3005);
    const float errval = 0.25;
    cout << "Will draw error band with value: " << errval << endl;
    for (int i = 0; i < bgerrgraph->GetN(); ++i) {
      Double_t x,y;
      bgerrgraph->GetPoint(i,x,y);
      bgerrgraph->SetPointError(i,float(rebinFactor)/2.,errval*y);
    }
  }

  ytitle += Form(" / %d GeV",rebinFactor);
  float max = 1.;

  if( overlayData ){

    datahist->Sumw2();
    // float xmin_original = datahist->GetXaxis()->GetXmin();
    // float xmax_original = datahist->GetXaxis()->GetXmax();
    // int nbins_original = datahist->GetNbinsX();
    // float xbinsize_original = (xmax_original - xmin_original)/float(nbins_original);
    // float xbinsize = (xmax - xmin)/float(nbins);
    // int rebinFactor = int(xbinsize/xbinsize_original);
    if (rebinFactor != 1) datahist->Rebin(rebinFactor);
    datahist->SetAxisRange(xmin,xmax-0.1*xbinsize,"X");
    Double_t err;
    float lastbin_and_overflow =  datahist->IntegralAndError(nbins,-1,err);
    if (!TString(histname).Contains("h_bdt")) {
      datahist->SetBinContent(nbins,lastbin_and_overflow);
      datahist->SetBinError(nbins,err);
    }

    max = datahist->GetMaximum() + datahist->GetBinError(datahist->GetMaximumBin());
    if( mctothist->GetMaximum() > max ) max = mctothist->GetMaximum();

    datahist->GetXaxis()->SetTitle(xtitle);
    datahist->GetXaxis()->SetTitleSize(0.05);
    datahist->GetXaxis()->SetLabelSize(0.04);
    datahist->GetXaxis()->SetTitleOffset(1.0);
    datahist->GetYaxis()->SetTitle(ytitle);
    datahist->GetYaxis()->SetTitleSize(0.05);
    datahist->GetYaxis()->SetTitleOffset(1.20);
    datahist->SetLineColor(kBlack);
    datahist->SetMarkerColor(kBlack);
    datahist->SetMarkerStyle(20);
    datahist->Draw("E1");
    mcstack->Draw("samehist");
    if (errband) bgerrgraph->Draw("same 2");
    if (nmcsig > 0) {
      for (int i = (int) mcsighist.size() - 1; i > -1 ; i--) {
	if (signorm > 0.) mcsighist[i]->Scale(signorm);
	if (stacksig) mcsighist[i]->Add(mctothist);
        if( mcsighist[i]->GetMaximum() > max ) max = mcsighist[i]->GetMaximum();
	mcsighist[i]->Draw("samehist");
      }
    }
    if( log ) {
      datahist->SetMaximum( 90 * max );
      if (TString(histname).Contains("mt2bl")) datahist->SetMaximum(200*max);
    }
    else      datahist->SetMaximum( 1.4 * max );
    datahist->Draw("sameE1");
    datahist->Draw("sameaxis");
    
    if(!log) {
      datahist->GetYaxis()->SetRangeUser(0.,1.4*max);
      if (TString(histname).Contains("dphi")) datahist->GetYaxis()->SetRangeUser(0.,3.0*max);
      if (TString(histname).Contains("bbmass") && tsdir.Contains("bbmasslast") && drawingSignal
	  && (tsdir.Contains("met1") || tsdir.Contains("bbmass_nm1")) && !stacksig ) {
	datahist->GetYaxis()->SetRangeUser(0.,2.5*max);
      }
    }

    datamchists.datahist = datahist;
    
  }
  else{

    max = mctothist->GetMaximum() + mctothist->GetBinError(mctothist->GetMaximumBin());
    if( log ) mctothist->SetMaximum( 90 * max );
    else      mctothist->SetMaximum( 1.4 * max );

    mctothist->GetXaxis()->SetTitle(xtitle);
    mctothist->GetYaxis()->SetTitle(ytitle);
    mctothist->GetYaxis()->SetTitleSize(0.05);
    mctothist->Draw("hist");
    mcstack->Draw("hist same");
    mctothist->Draw("hist same axis");
    if (errband) bgerrgraph->Draw("same 2");
    if (nmcsig > 0) {
      for (int i = (int) mcsighist.size() - 1; i > -1 ; i--) {
	if (signorm > 0.) mcsighist[i]->Scale(signorm);
	if (stacksig) mcsighist[i]->Add(mctothist);
        if( mcsighist[i]->GetMaximum() > max ) max = mcsighist[i]->GetMaximum();
	mcsighist[i]->Draw("samehist");
      }
    }

    if(!log) mctothist->GetYaxis()->SetRangeUser(0.,1.4*max);
  }

  datamchists.mctothist = mctothist;

  if( drawLegend ){
    TLegend* myleg = 0;
    if (nmcsig > 0) myleg = getLegend( labels , overlayData, 0.57, 0.45, 0.88, 0.94, signorm, errband );
    else myleg = getLegend( labels , overlayData );
    myleg->Draw();
  }

  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.03);
  TString label = "CMS " + TString(cmslabel);
  text->DrawLatex(0.2,0.88,label);
  //text->DrawLatex(0.2,0.83,"0.98 fb^{-1} at #sqrt{s} = 7 TeV");
  text->DrawLatex(0.2,0.83,"#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");

  // if     ( TString(flavor).Contains("ee")  ) text->DrawLatex(0.2,0.78,"Events with ee");
  // else if( TString(flavor).Contains("mm")  ) text->DrawLatex(0.2,0.78,"Events with #mu#mu");
  // else if( TString(flavor).Contains("em")  ) text->DrawLatex(0.2,0.78,"Events with e#mu");
  // else if( TString(flavor).Contains("all") ) text->DrawLatex(0.2,0.78,"Events with ee/#mu#mu/e#mu");
  // else if( TString(flavor).Contains("sf")  ) text->DrawLatex(0.2,0.78,"Events with ee/#mu#mu");
  // else if( TString(flavor).Contains("3l")  ) text->DrawLatex(0.2,0.78,"Events with eee/ee#mu/#mu#mue/#mu#mu#mu");
  // else if( TString(flavor).Contains("e")  ) text->DrawLatex(0.2,0.78,"Events with e");
  // else if( TString(flavor).Contains("m")  ) text->DrawLatex(0.2,0.78,"Events with #mu");
  // else if( TString(flavor).Contains("sl")  ) text->DrawLatex(0.2,0.78,"Events with e/#mu");
  // else                                       text->DrawLatex(0.2,0.78,"Events with e/#mu");
  //else                                       text->DrawLatex(0.2,0.78,"Events with e#mu");
  //else                                       text->DrawLatex(0.2,0.78,"Events with #mu#mu");
  // label for specific control regions
  float yval = 0.78;
  if (tsdir.Contains("cr1_")) text->DrawLatex(0.2,yval,"M_{b#bar{b}} > 150 GeV");
  else if (tsdir.Contains("cr8_")) text->DrawLatex(0.2,yval,"M_{b#bar{b}} < 100 GeV");
  //  else if (tsdir.Contains("cr14_")) text->DrawLatex(0.2,yval,"M_{b#bar{b}} < 100 GeV or M_{b#bar{b}} > 150 GeV");
  else if (tsdir.Contains("_met100")) text->DrawLatex(0.2,yval,"E_{T}^{miss} > 100 GeV");
  else if (tsdir.Contains("_met125")) text->DrawLatex(0.2,yval,"E_{T}^{miss} > 125 GeV");
  else if (tsdir.Contains("_met150")) text->DrawLatex(0.2,yval,"E_{T}^{miss} > 150 GeV");
  else if (tsdir.Contains("bbmasslast") && tsdir.Contains("_bbmass_nm1")) text->DrawLatex(0.2,yval,"E_{T}^{miss} > 175 GeV");
  else if (tsdir.Contains("cr23_")) text->DrawLatex(0.2,yval,"CR-2l");
  else if (tsdir.Contains("cr5_")) text->DrawLatex(0.2,yval,"CR-0b");

  // // draw lines for signal region -- doesn't look so good...
  // if (fullhistname.Contains("sig_bbmasslast") && fullhistname.Contains("h_bbmass") ) {
  //   TLine* line_high = new TLine(150.,0.,150.,max*1.05);
  //   line_high->SetLineColor(kGray+2);
  //   line_high->SetLineWidth(2);
  //   line_high->SetLineStyle(2);
  //   line_high->Draw("same");
  //   TLine* line_low = new TLine(100.,0.,100.,max*1.05);
  //   line_low->SetLineColor(kGray+2);
  //   line_low->SetLineWidth(2);
  //   line_low->SetLineStyle(2);
  //   line_low->Draw("same");
  // }

  if( residual ){
    fullpad->cd();

    respad = new TPad("respad","respad",0,0.8,1,0.98);
    respad->SetRightMargin(0.05);
    respad->Draw();
    respad->cd();
    respad->SetTopMargin(0.05);

    gPad->SetGridy();

    //    if (false) {
    if ((TString(histname).Contains("sumpt") || TString(histname).Contains("dilpt")) 
        && !TString(histname).Contains("rapidity") ) {
      // int newnbins = 18;
      // Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,350,400,520};

      // int newnbins = 13;
      // Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,250,350,500};

      // binning for ttbar high stats sumpt plots
      int newnbins = 10;
      Double_t xbins[] = {0,20,40,60,80,100,150,200,250,350,500};

      // used for most plots of sumpt and dilpt
      //int newnbins = 14;
      //Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,250,350,500,800};

      // new z+jets binning for flatter systs
      // int newnbins = 15;
      // Double_t xbins[] = {0,20,40,60,80,100,120,140,160,180,200,250,300,400,500,800};

      // for monojet dilpt plot
      // int newnbins = 7;
      // Double_t xbins[] = {0,150,200,250,300,350,500,800};

      // int newnbins = 6;
      // Double_t xbins[] = {0,50,100,150,200,300,520};

      // int newnbins = 8;
      // Double_t xbins[] = {0,20,40,60,80,100,150,250,520};

      // int newnbins = 3;
      // Double_t xbins[] = {0,100,200,520};

      // int newnbins = 2;
      // Double_t xbins[] = {0,150,500};

      ratio = (TH1F*) datahist->Rebin(newnbins,Form("%s_ratio",datahist->GetName()),xbins);
      TH1F* mcdenom = (TH1F*) mctothist->Rebin(newnbins,Form("%s_denom",mctothist->GetName()),xbins);
      ratio->Divide(mcdenom);
    } else {
      ratio = (TH1F*) datahist->Clone(Form("%s_ratio",datahist->GetName()));
      ratio->Divide(mctothist);
    }

    ratio->GetYaxis()->SetTitleOffset(0.3);
    ratio->GetYaxis()->SetTitleSize(0.2);
    ratio->GetYaxis()->SetNdivisions(5);
    ratio->GetYaxis()->SetLabelSize(0.15);
    //ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetRangeUser(0.001,2.0);
    ratio->GetYaxis()->SetTitle("Data/SM  ");
    ratio->GetXaxis()->SetLabelSize(0);
    ratio->GetXaxis()->SetTitleSize(0);
    ratio->SetMarkerSize(0.7);

    TF1* fpol1 = new TF1("fpol1","pol1",datahist->GetXaxis()->GetXmin(),datahist->GetXaxis()->GetXmax());
    if( fit ) ratio->Fit(fpol1,"R");

    ratio->Draw();

    gratio = new TGraphErrors(ratio);
    for (int i=0; i<ratio->GetNbinsX(); ++i) {
      gratio->SetPointError(i,ratio->GetBinWidth(i+1)/2.,gratio->GetErrorY(i));
    }
    gratio->SetLineColor(kBlack);
    gratio->SetMarkerColor(kBlack);
    gratio->SetMarkerStyle(20);
    gratio->Draw("p0same");

    TLine line;
    line.SetLineWidth(1);
    //    line.DrawLine(ratio->GetXaxis()->GetXmin(),1,ratio->GetXaxis()->GetXmax(),1);
    line.DrawLine(xmin,1,xmax,1);

    datamchists.ratio = ratio;
  }

  datamchists.comppad = fullpad;
  //  return datamchists;
  //  return ratio;
  return gratio;
}

//____________________________________________________________________________
TGraphAsymmErrors* makeBand(TH1F* centhist, TH1F* uphist, TH1F* dnhist) {

  int nbins = centhist->GetNbinsX();

  // makes a graph with the correct central values
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nbins);

  // need to set correct upper and lower errors to create band

  // assume all bins are same width..
  float binwidth = centhist->GetBinWidth(1);

  for (int ibin=1; ibin < nbins+1; ++ibin) {
    float valx = centhist->GetBinCenter(ibin);
    float valy = centhist->GetBinContent(ibin);
    graph->SetPoint(ibin-1,valx,valy);
    float uperr = valy - uphist->GetBinContent(ibin);
    float dnerr = dnhist->GetBinContent(ibin) - valy;
    float err = max(uperr,dnerr);
    if (uperr < 0.) uperr = dnerr;
    else if (dnerr < 0.) dnerr = uperr;
    //    graph->SetPointError(ibin-1,binwidth/2.,binwidth/2.,err,err);
    graph->SetPointError(ibin-1,binwidth/2.,binwidth/2.,uperr,dnerr);
  }

  return graph;
}

//____________________________________________________________________________
TGraphAsymmErrors* makeBand(TGraphErrors* centgraph, TGraphErrors* upgraph, TGraphErrors* dngraph) {

  int nbins = centgraph->GetN();
  TGraphAsymmErrors* graph = new TGraphAsymmErrors(nbins);

  for (int ibin=0; ibin < nbins; ++ibin) {
    Double_t valx, valy;
    centgraph->GetPoint(ibin,valx,valy);
    graph->SetPoint(ibin,valx,valy);
    Double_t errx = centgraph->GetErrorX(ibin);
    Double_t upy, dny;
    upgraph->GetPoint(ibin,valx,upy);
    dngraph->GetPoint(ibin,valx,dny);
    float uperr = valy - upy;
    float dnerr = dny - valy;
    float err = max(uperr,dnerr);
    if (uperr < 0.) uperr = dnerr;
    else if (dnerr < 0.) dnerr = uperr;
    //    graph->SetPointError(ibin,errx,errx,err,err);
    graph->SetPointError(ibin,errx,errx,uperr,dnerr);
  }

  return graph;
}

// //____________________________________________________________________________
// TCanvas* compareNormalized(std::string histname, TFile* fsig, TFile* fsig2, TFile* fbg, int rebin, bool norm) {

//   TCanvas* c = new TCanvas(Form("c_%s",histname.c_str()),Form("c_%s",histname.c_str()));

//   TH1F* hsig = (TH1F*)((TH1F*)fsig->Get(histname.c_str()))->Clone(Form("%s_fsig",histname.c_str()));
//   TH1F* hsig2 = 0;
//   if (fsig2) hsig2 = (TH1F*)((TH1F*)fsig2->Get(histname.c_str()))->Clone(Form("%s_fsig2",histname.c_str()));
//   TH1F* hbg = (TH1F*)((TH1F*)fbg->Get(histname.c_str()))->Clone(Form("%s_fbg",histname.c_str()));

//   hsig->SetMarkerColor(kBlue);
//   hsig->SetLineColor(kBlue);
//   hsig->SetLineWidth(2);
//   if (hsig2) {
//     hsig2->SetMarkerColor(kGreen);
//     hsig2->SetLineColor(kGreen);
//     hsig2->SetLineWidth(2);
//   }
//   hbg->SetMarkerColor(kRed);
//   hbg->SetLineColor(kRed);
//   hbg->SetLineWidth(2);

//   if (rebin > 1) {
//     hsig->Rebin(rebin);
//     if (hsig2) hsig2->Rebin(rebin);
//     hbg->Rebin(rebin);
//   }

//   TH1F* hsig_norm = 0;
//   TH1F* hsig2_norm = 0;
//   TH1F* hbg_norm = 0;

//   if (norm) {
//     hsig_norm = (TH1F*)hsig->DrawNormalized("histe");
//     if (hsig2) hsig2_norm = (TH1F*)hsig2->DrawNormalized("histe same");
//     hbg_norm = (TH1F*)hbg->DrawNormalized("histe same");
//   } else {
//     hsig_norm = hsig;
//     if (hsig2) hsig2_norm = hsig2;
//     hbg_norm = hbg;
//     hsig_norm->Draw("histe");
//     if (hsig2) hsig2_norm->Draw("histe same");
//     hbg_norm->Draw("histe same");
//   }

//   if (hbg_norm->GetMaximum() > hsig_norm->GetMaximum()) {
//     hsig_norm->GetYaxis()->SetRangeUser(1E-4,1.1*hbg_norm->GetMaximum());
//   }

//   if (hsig2_norm && (hsig2_norm->GetMaximum() > hsig_norm->GetMaximum())) {
//     hsig_norm->GetYaxis()->SetRangeUser(1E-4,1.1*hsig2_norm->GetMaximum());
//   }

//   // std::string label_sig = "Wbb+jets";
//   // std::string label_sig2 = "TChiwh_150_1";
//   // std::string label_bg = "W+jets";
//   std::string label_sig = "TChiwh_250_1";
//   std::string label_sig2 = "TChiwh_150_1";
//   std::string label_bg = "all bkg";

//   TLegend *leg = new TLegend(0.66,0.77,0.93,0.9);
//   leg->SetFillColor(0);
//   leg->AddEntry(hsig_norm,label_sig.c_str(),"l");
//   if (hsig2_norm) leg->AddEntry(hsig2_norm,label_sig2.c_str(),"l");
//   leg->AddEntry(hbg_norm,label_bg.c_str(),"l");
//   leg->Draw("same");

//   gPad->Modified();

//   return c;
// }

//____________________________________________________________________________
TCanvas* compareNormalized(std::string histname, TFile* f1, std::string label1, TFile* f2, std::string label2, int rebin, bool norm, TFile* f3, std::string label3, TFile* f4, std::string label4) {

  TH1F* h1 = (TH1F*)f1->Get(histname.c_str());
  TH1F* h2 = (TH1F*)f2->Get(histname.c_str());
  TH1F* h3 = 0;
  TH1F* h4 = 0;
  if (f3) h3 = (TH1F*)f3->Get(histname.c_str());
  if (f4) h4 = (TH1F*)f4->Get(histname.c_str());

  TCanvas* c = compareNormalized(h1,label1,h2,label2,rebin,norm,h3,label3,h4,label4);

  return c;
}

//____________________________________________________________________________
TCanvas* compareNormalized(TH1F* h1, std::string label1, TH1F* h2, std::string label2, int rebin, bool norm, TH1F* h3, std::string label3, TH1F* h4, std::string label4) {

  TCanvas* c = new TCanvas(Form("c_%s",h1->GetName()),Form("c_%s",h1->GetName()));

  TH1F* h1_clone = (TH1F*)h1->Clone(Form("%s_clone",h1->GetName()));
  TH1F* h2_clone = (TH1F*)h2->Clone(Form("%s_clone",h2->GetName()));
  TH1F* h3_clone = 0;
  TH1F* h4_clone = 0;
  if (h3) h3_clone = (TH1F*)h3->Clone(Form("%s_clone",h3->GetName()));
  if (h4) h4_clone = (TH1F*)h4->Clone(Form("%s_clone",h4->GetName()));

  h1_clone->SetMarkerColor(kRed);
  h1_clone->SetLineColor(kRed);
  h1_clone->SetLineWidth(2);
  h2_clone->SetMarkerColor(kBlue);
  h2_clone->SetLineColor(kBlue);
  h2_clone->SetLineWidth(2);
  if (h3_clone) {
    h3_clone->SetMarkerColor(7);
    h3_clone->SetLineColor(7);
    h3_clone->SetLineWidth(2);
  }
  if (h4_clone) {
    h4_clone->SetMarkerColor(8);
    h4_clone->SetLineColor(8);
    h4_clone->SetLineWidth(2);
  }

  if (rebin > 1) {
    h1_clone->Rebin(rebin);
    if (h3_clone) h3_clone->Rebin(rebin);
    if (h4_clone) h4_clone->Rebin(rebin);
    h2_clone->Rebin(rebin);
  }

  TH1F* h1_norm = 0;
  TH1F* h2_norm = 0;
  TH1F* h3_norm = 0;
  TH1F* h4_norm = 0;

  if (norm) {
    h1_norm = (TH1F*)h1_clone->DrawNormalized("histe");
    h2_norm = (TH1F*)h2_clone->DrawNormalized("histe same");
    if (h3_clone) h3_norm = (TH1F*)h3_clone->DrawNormalized("histe same");
    if (h4_clone) h4_norm = (TH1F*)h4_clone->DrawNormalized("histe same");
  } else {
    h1_norm = h1_clone;
    h2_norm = h2_clone;
    if (h3_clone) h3_norm = h3_clone;
    if (h3_clone) h3_norm = h3_clone;
    h1_norm->Draw("histe");
    h2_norm->Draw("histe same");
    if (h3_clone) h3_norm->Draw("histe same");
    if (h4_clone) h4_norm->Draw("histe same");
  }

  if (h2_norm->GetMaximum() > h1_norm->GetMaximum()) {
    h1_norm->GetYaxis()->SetRangeUser(1E-4,1.1*h2_norm->GetMaximum());
  }

  if (h3_norm && (h3_norm->GetMaximum() > h1_norm->GetMaximum())) {
    h1_norm->GetYaxis()->SetRangeUser(1E-4,1.1*h3_norm->GetMaximum());
  }

  if (h4_norm && (h4_norm->GetMaximum() > h1_norm->GetMaximum())) {
    h1_norm->GetYaxis()->SetRangeUser(1E-4,1.1*h4_norm->GetMaximum());
  }

  TLegend *leg = new TLegend(0.57,0.69,0.93,0.9);
  leg->SetFillColor(0);
  leg->AddEntry(h1_norm,label1.c_str(),"l");
  leg->AddEntry(h2_norm,label2.c_str(),"l");
  if (h3_norm) leg->AddEntry(h3_norm,label3.c_str(),"l");
  if (h4_norm) leg->AddEntry(h4_norm,label4.c_str(),"l");
  leg->Draw("same");

  gPad->Modified();

  return c;
}

//______________________________________________________________________________
TH1F* cumulate (TH1F* in, bool increasing) {
  TH1F* h_out = new TH1F(in->GetName() + TString("tmp"), in->GetTitle(), in->GetNbinsX(),
			 in->GetBinLowEdge(1), in->GetBinLowEdge(in->GetNbinsX() + 1));
  h_out->Sumw2();
  h_out->SetFillColor(in->GetFillColor());
  h_out->SetFillStyle(in->GetFillStyle());
  h_out->SetLineStyle(in->GetLineStyle());
  h_out->SetLineColor(in->GetLineColor());
  h_out->SetLineWidth(in->GetLineWidth());
  double sum = 0;
  double err2 = 0;
  if (increasing) {
    for (int j = 0; j <= in->GetNbinsX() + 1; ++j) {
      sum += in->GetBinContent(j);
      err2 += in->GetBinError(j)*in->GetBinError(j);
      h_out->SetBinContent(j, sum);
      h_out->SetBinError(j, sqrt(err2));
    }
  } else {
    for (int j = in->GetNbinsX() + 1; j >= 0; --j) {
      sum += in->GetBinContent(j);
      err2 += in->GetBinError(j)*in->GetBinError(j);
      h_out->SetBinContent(j, sum);
      h_out->SetBinError(j, sqrt(err2));
    }
  }
  return h_out;
}

//______________________________________________________________________________
TGraphErrors* eff_rej (TH1F* signal, TH1F* background, bool normalize, bool increasing, bool print, float sf) {

  TH1F* sig = (TH1F*)signal->Clone("h_tmp_s");
  if (normalize) sig->Scale(1 / sig->Integral(0, sig->GetNbinsX() + 1));
  TH1F* bg = (TH1F*)background->Clone("h_tmp_bg");
  if (normalize) bg->Scale(1 / bg->Integral(0, bg->GetNbinsX() + 1));
  TH1F* sig_cum = cumulate(sig, increasing);
  TH1F* bg_cum = cumulate(bg, increasing);
  TGraphErrors* ret = new TGraphErrors(signal->GetNbinsX());
  for (int i = 1; i <= signal->GetNbinsX(); ++i) {
    const double x = sig_cum->GetBinCenter(i);
    const double s = sig_cum->GetBinContent(i) * sf;
    if (!normalize && (s < 3.0)) continue;
    const double b = bg_cum->GetBinContent(bg_cum->FindBin(x)) * sf;
    const double errs = sig_cum->GetBinError(i);
    const double errb = bg_cum->GetBinError(bg_cum->FindBin(x));
    ret->SetPoint(i - 1, 1.-b, s); // gotta love offsets
    ret->SetPointError(i - 1, errb, errs);
    if (print && s > 0) {
      std::cout << "x value: " << x << ", s: " << s << ", b: " << b;
      if (b > 0. && !normalize) std::cout << ", s/b: " << s/b << ", s/sqrt(b): " << s/sqrt(b)
					  << ", sig: " << sqrt(2*( (s+b) * log(1+(s/b)) - s) );
      std::cout << std::endl;
    }
  }
  return ret;
}

// ---------------------------------------------------
// function to make graph of S/sqrt(B)
// note: set increasing to false to cut var > XX
TGraph* s_over_rootb (TH1F* signal, TH1F* background, bool increasing, bool do_s_over_b, bool print, float sf) {

  TH1F* sig = (TH1F*)signal->Clone("h_tmp_s");
  TH1F* bg = (TH1F*)background->Clone("h_tmp_bg");
  TH1F* sig_cum = cumulate(sig, increasing);
  TH1F* bg_cum = cumulate(bg, increasing);
  TGraph* ret = new TGraphErrors(signal->GetNbinsX());
  for (int i = 1; i <= signal->GetNbinsX(); ++i) {
    const double x = sig_cum->GetBinCenter(i);
    const double s = sig_cum->GetBinContent(i) * sf;
    if (s < 3.0) continue;
    const double b = bg_cum->GetBinContent(bg_cum->FindBin(x)) * sf;
    //    const double errs = sig_cum->GetBinError(i);
    //    const double errb = bg_cum->GetBinError(bg_cum->FindBin(x));
    double s_over_rootb = 0.;
    if (b > 0.) {
      //      s_over_rootb = s/sqrt(b);
      s_over_rootb = s/sqrt(b + 0.25*0.25*b*b);
      if (do_s_over_b) s_over_rootb = s/b;
    }
    ret->SetPoint(i - 1, x, s_over_rootb); // gotta love offsets
    //    ret->SetPointError(i - 1, errb, errs);
    if (print) {
      std::cout << "x value: " << x << ", s: " << s << ", b: " << b;
      if (b > 0.) std::cout << ", s/b: " << s/b << ", s/sqrt(b): " << s/sqrt(b)
			    << ", s/sqrt(b+syst): " << s/sqrt(b + 0.25*0.25*b*b)
			    << ", sig: " << sqrt(2*( (s+b) * log(1+(s/b)) - s) );
      std::cout << std::endl;
    }
  }
  return ret;
}

//____________________________________________________________________
TLegend* init_legend(float x1, float y1, float x2, float y2) {
  // intializes a legend with white background, no border
  TLegend* l = new TLegend(x1,y1,x2,y2);
  l->SetLineColor(0);
  l->SetFillColor(0);
  l->SetShadowColor(0);
  return l;
}

//______________________________________________________________________________
TLegend* legendize(TCanvas* c, const TString& opt, const TString& label1, const TString& label2, const TString& label3 ) {

  TString labels[3] = {label1,label2,label3};
  // makes a basic legend by finding all TH1s and graphs in a canvas
  TList* clist = c->GetListOfPrimitives();
  TLegend* l = init_legend();
  // iterate over list of primitives
  int nhists = 0;
  TIter next(clist);
  TObject* obj = 0;
  while ((obj = next())) {
    // check for TH1 or TGraph
    if (TString(obj->ClassName()).Contains("TH1")) {
      TH1* h = (TH1*)obj;
      if (nhists <= 2) {
	l->AddEntry(h,labels[nhists],opt);
      } else {
	l->AddEntry(h,"",opt);
      }
      nhists++;
    } else if (TString(obj->ClassName()).Contains("TGraph")) {
      TGraph* h = (TGraph*)obj;
      if (nhists <= 2) {
	l->AddEntry(h,labels[nhists],opt);
      } else {
	l->AddEntry(h,"",opt);
      }
      nhists++;
    }
  }
  l->Draw("same");

  return l;
}

//______________________________________________________________________________
// returns the error on C = A*B (or C = A/B)
float err_mult(float A, float B, float errA, float errB, float C) {
  return sqrt(C*C*(pow(errA/A,2) + pow(errB/B,2)));
}

//______________________________________________________________________________
// returns the error on C = A/(A+B) 
//  note that if A and B are integers, simplifies to sqrt((C * (1-C)) / (A+B))
//  or thinking of an efficiency, sqrt((eff * (1-eff)) / N)
float err_binomial(float A, float B, float errA, float errB) {
  return (1/pow(A+B,2)) * sqrt(pow(B*errA,2) + pow(A*errB,2));
}

//______________________________________________________________________________
// integrates MET above various values and makes graph
TGraphErrors* makeMETGraph(TH1F* hdata, TH1F* hmc, float xoffset, bool errOnData) {

  // const int ncuts = 4;
  // const float vals[ncuts] = {50.,100.,150.,175.};
  const int ncuts = 6;
  const float vals[ncuts] = {50.,75.,100.,125.,150.,175.};
  const bool print = true;

  TGraphErrors* g = new TGraphErrors(ncuts);

  for (int i = 0; i < ncuts; ++i) {
    Double_t err_data = 0;
    Double_t err_mc = 0;
    float n_data = hdata->IntegralAndError(vals[i]+1,-1,err_data);
    float n_mc = hmc->IntegralAndError(vals[i]+1,-1,err_mc);
    float ratio = n_data/n_mc;
    float err = err_mult(n_data,n_mc,err_data,err_mc,ratio);
    g->SetPoint(i, i+0.5+xoffset, ratio);
    g->SetPointError(i, 0., err);
    if (print) {
      std::cout << "MET " << vals[i] << ": MC, data, ratio: " << Form("%.1f",n_mc) << " $\\pm$ " << Form("%.1f",err_mc);
      if (errOnData) {
	std::cout << " & " << Form("%.1f",n_data) << " $\\pm$ " << Form("%.1f",err_data) << " & " ;
      } else {
	std::cout << " & " << Form("%d",int(n_data)) << " & " ;
      }
      std::cout	<< Form("%.2f",ratio) << " $\\pm$ " << Form("%.2f",err) << std::endl;
    }
  }

  return g;
}

//______________________________________________________________________________
// returns mt tail fraction 
float getMTTailFrac(TFile* f, TString dir, Double_t& err, bool print = true) {

  const int cut = 101;

  TH1F* h = (TH1F*) f->Get(dir+"/h_lep1mt");
  if (!h) {
    std::cout << "ERROR: problem retreiving h_lep1mt from dir: " << dir << std::endl;
    return -1.;
  }

  Double_t tail_err, all_err;
  float tail = h->IntegralAndError(cut,-1,tail_err);
  float all = h->IntegralAndError(0,-1,all_err);

  err = err_mult(tail,all,tail_err,all_err,tail/all);
  if (print) std::cout << "tail fraction: " << Form("%.3f",tail/all) << " $\\pm$ " << Form("%.3f",err) << std::endl;
  return tail/all;
}

