#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TLatex.h"
#include <fstream>
#include "tdrStyle.C"
#include "CommonFunctions.C"
#include <vector>

using namespace std;

TH1F* hdefaultgenlep15;
TH1F* hscaleup;
TH1F* hscaledown;
TH1F* hdefault;
TH1F* htptrw;
TH2F* hdefaultgenlep152d;
TH2F* hscaleup2d;
//TH2F* hscaledown2d;
TH2F* hdefaultgenlep152drebinned;
TH2F* hscaleup2drebinned;
TH2F* hscaledown2drebinned;

void GetAfb(TH1F* h, Double_t &afb, Double_t  &afberr){
 
  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"");
  //event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"width");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"width");
  event_total = event_plus + event_minus;
  
  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
		     + event_minus*event_minus*event_plus_err*event_plus_err)/
		  (event_total*event_total*event_total*event_total));

}

//void comparisonplots(TString histname = "topSpinCorrGen", bool drawnorm = true, TString FName1 = "denominator_powheg_stat1/results/hist_noCuts.root", TString FName2 = "denominator_powhegscfix_stat1/results/hist_noCuts.root", TString FName3 = "denominator_mcatnlo/results/hist_noCuts.root"){
void comparescaleandtptvariations(TString histname = "ttpTGen", bool drawnorm = true, TString FName1 = "results/hist_noCuts.root", TString FName2 = "results/hist_noCuts_scaleup.root", TString FName3 = "results/hist_noCuts_scaledown.root", TString FName4 = "results/hist_noCuts_top_pt_rw.root", TString FName5 = "results/hist_noCuts_default.root"){
  setTDRStyle();

  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  hdefaultgenlep15 = (TH1F*)f_1->Get(Form("ttdil_h%s_allj_all", histname.Data())); 
  hdefaultgenlep152d = (TH2F*)f_1->Get(Form("ttdil_h%s_allj_all", histname.Data())); 
    
  std::cout << "Opening " << FName2.Data() << "\n";  
  TFile *f_2         = TFile::Open(FName2.Data());  
  hscaleup = (TH1F*)f_2->Get(Form("ttdil_h%s_allj_all", histname.Data()));
  hscaleup2d = (TH2F*)f_2->Get(Form("ttdil_h%s_allj_all", histname.Data())); 
  
  std::cout << "Opening " << FName3.Data() << "\n";  
  TFile *f_3         = TFile::Open(FName3.Data());  
  hscaledown = (TH1F*)f_3->Get(Form("ttdil_h%s_allj_all", histname.Data()));
  
  std::cout << "Opening " << FName4.Data() << "\n";  
  TFile *f_4         = TFile::Open(FName4.Data());  
  htptrw = (TH1F*)f_4->Get(Form("ttdil_h%s_allj_all", histname.Data()));

  std::cout << "Opening " << FName5.Data() << "\n";  
  TFile *f_5         = TFile::Open(FName5.Data());  
  hdefault = (TH1F*)f_5->Get(Form("ttdil_h%s_allj_all", histname.Data()));

  std::cout << "Opened " << Form("ttdil_h%s_allj_all", histname.Data()) << " and "<< Form("ttdil_h%s_allj_all", histname.Data()) <<"\n";

  //hscaledown->Print();
  
  Double_t pi = 3.141592653589793;
  bool printmean = false;


  hdefaultgenlep15->Scale(1./hdefaultgenlep15->Integral());
  hscaleup->Scale(1./hscaleup->Integral());
  hscaledown->Scale(1./hscaledown->Integral());
  hdefault->Scale(1./hdefault->Integral());
  htptrw->Scale(1./htptrw->Integral());


  hdefaultgenlep15->Rebin(8);
  hscaleup->Rebin(8);
  hscaledown->Rebin(8);
  hdefault->Rebin(8);
  htptrw->Rebin(8);

  hscaleup->Divide(hdefaultgenlep15);
  hscaledown->Divide(hdefaultgenlep15);
  htptrw->Divide(hdefault);




  double Asym1,Asym2,Asym3;
  double Asym1err,Asym2err,Asym3err;

    if(histname.Contains("topSpinCorr")){
      Asym1 = hscaledown->GetMean();
      Asym2 = htptrw->GetMean();
      Asym3 = hscaleup->GetMean();
      Asym1err = hscaledown->GetMeanError();
      Asym2err = htptrw->GetMeanError();
      Asym3err = hscaleup->GetMeanError();
      cout<<"mean: "<<9*Asym1<<" +/- "<<9*Asym1err<<", "<<9*Asym2<<" +/- "<<9*Asym2err<<", "<<9*Asym3<<" +/- "<<9*Asym3err<<endl;

      GetAfb(hscaledown,Asym1, Asym1err);
      GetAfb(htptrw,Asym2, Asym2err);
      GetAfb(hscaleup,Asym3, Asym3err);
      cout<<"asym: "<<4*Asym1<<" +/- "<<4*Asym1err<<", "<<4*Asym2<<" +/- "<<4*Asym2err<<", "<<4*Asym3<<" +/- "<<4*Asym3err<<endl;
    }
    else if(histname.Contains("Cos")){
      Asym1 = hscaledown->GetMean();
      Asym2 = htptrw->GetMean();
      Asym3 = hscaleup->GetMean();
      Asym1err = hscaledown->GetMeanError();
      Asym2err = htptrw->GetMeanError();
      Asym3err = hscaleup->GetMeanError();
      cout<<"mean: "<<3*Asym1<<" +/- "<<3*Asym1err<<", "<<3*Asym2<<" +/- "<<3*Asym2err<<", "<<3*Asym3<<" +/- "<<3*Asym3err<<endl;

      GetAfb(hscaledown,Asym1, Asym1err);
      GetAfb(htptrw,Asym2, Asym2err);
      GetAfb(hscaleup,Asym3, Asym3err);
      cout<<"asym: "<<2*Asym1<<" +/- "<<2*Asym1err<<", "<<2*Asym2<<" +/- "<<2*Asym2err<<", "<<2*Asym3<<" +/- "<<2*Asym3err<<endl;
    }
    else {
      Asym1 = hscaledown->GetMean();
      Asym2 = htptrw->GetMean();
      Asym3 = hscaleup->GetMean();
      Asym1err = hscaledown->GetMeanError();
      Asym2err = htptrw->GetMeanError();
      Asym3err = hscaleup->GetMeanError();
      cout<<"mean: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;

      //if( !(histname.Contains("ttMass") || histname.Contains("ttpT")) ){
        GetAfb(hscaledown,Asym1, Asym1err);
        GetAfb(htptrw,Asym2, Asym2err);
        GetAfb(hscaleup,Asym3, Asym3err);
        cout<<"asym: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;
      //}
      //else printmean = true;
    }




  TString accepthistname = "scale_and_toppt_ratios_";
  accepthistname += histname;

  
    
  htptrw->SetLineColor(kBlack);
  htptrw-> SetFillColor(0);
  htptrw->SetMarkerColor(kBlack);
  hscaleup->SetLineColor(kRed);
  hscaleup->SetMarkerColor(kRed);
  hscaleup-> SetFillColor(0);
  hscaledown->SetLineColor(kBlue);
  hscaledown->SetMarkerColor(kBlue);
  hscaledown-> SetFillColor(0);
  
  gStyle->SetPaintTextFormat("6.4f");

  TCanvas *c1 = new TCanvas();
  c1->cd();
  //gPad->SetLogy();
  
  hscaledown->SetMaximum(1.25*hscaledown->GetMaximum());
  if(hscaledown->GetMinimum() <0.15 *hscaledown->GetMaximum() ) hscaledown->SetMinimum(0.);  
  if(hscaledown->GetMinimum() > 0.) hscaledown->SetMinimum(0.75*hscaledown->GetMinimum() );  
  
  if(!drawnorm){
  	hscaledown->Draw("hist TEXT00E");
  }
  else{
  	hscaledown->Draw("hist");  	
  	htptrw->Draw("histsame");
  	hscaleup->Draw("histsame");
  }

  TLegend *leg = new TLegend(0.47,0.18,0.69,0.27);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  if(drawnorm) leg->AddEntry(htptrw, "topptrw/default","l");
  leg->AddEntry(hscaledown, "scaledown/default","l");
  //if(drawnorm) leg->AddEntry(htptrw, "powheg (stat1 leps)","l");
  //if(drawnorm) leg->AddEntry(hscaleup, "powheg SC fix (stat1 leps)","l");
  if(drawnorm) leg->AddEntry(hscaleup, "scaleup/default","l");
  
  leg->Draw("same");


  TPaveText *pt1 = new TPaveText(0.18, 0.77, 0.40, 0.92, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);
  
  TText *blah;
  blah = pt1->AddText("CMS Simulation, #sqrt{s}=8 TeV");
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  

  blah = pt1->AddText("");

  TString Asym2_temp = formatFloat(Asym2,"%6.4f");
  Asym2_temp.ReplaceAll(" " , "" );
  if(printmean) Asym2_temp = TString("   Mean: ") +  Asym2_temp;
  else Asym2_temp = TString("   Asym: ") +  Asym2_temp;
  blah = pt1->AddText(Asym2_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlack);  

  TString Asym1_temp = formatFloat(Asym1,"%6.4f");
  Asym1_temp.ReplaceAll(" " , "" );
  if(printmean) Asym1_temp = TString("   Mean: ") +  Asym1_temp;
  else Asym1_temp = TString("   Asym: ") +  Asym1_temp;
  blah = pt1->AddText(Asym1_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlue);  

  TString Asym3_temp = formatFloat(Asym3,"%6.4f");
  Asym3_temp.ReplaceAll(" " , "" );
  if(printmean) Asym3_temp = TString("   Mean: ") +  Asym3_temp;
  else Asym3_temp = TString("   Asym: ") +  Asym3_temp;
  blah = pt1->AddText(Asym3_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kRed);  


  
  pt1->Draw();


  c1->Print(Form("%s.pdf", accepthistname.Data()));


  f_1->Close();
  f_2->Close();
  f_3->Close();
  //output->Close();
  
}
