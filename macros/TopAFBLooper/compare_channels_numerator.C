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
#include "denominator/tdrStyle.C"
#include "denominator/CommonFunctions.C"
#include <vector>

using namespace std;

TH1F* hmuel;
TH1F* helmu;
TH1F* hdimu;
TH1F* hdiel;

TH1F* hmuel_numerator;
TH1F* helmu_numerator;
TH1F* hdimu_numerator;
TH1F* hdiel_numerator;


void GetAfb(TH1F* h, Double_t &afb, Double_t  &afberr){
 
  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  //event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"width");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"width");
  event_total = event_plus + event_minus;
  
  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
		     + event_minus*event_minus*event_plus_err*event_plus_err)/
		  (event_total*event_total*event_total*event_total));

}

void compare_channels_numerator(TString histname = "topSpinCorr", bool drawnorm = true) {
  setTDRStyle();

  TString observablename = "";

  if(histname=="lepChargeAsym") observablename="lep_charge_asymmetry";
  if(histname=="lepAzimAsym2") observablename="lep_azimuthal_asymmetry2";
  if(histname=="lepAzimAsym") observablename="lep_azimuthal_asymmetry";
  if(histname=="lepPlusCosTheta") observablename="lepPlus_costheta_cms";
  if(histname=="lepMinusCosTheta") observablename="lepMinus_costheta_cms";
  if(histname=="lepCosTheta") observablename="lep_costheta_cms";
  if(histname=="topSpinCorr") observablename="top_spin_correlation";
  if(histname=="rapiditydiffMarco") observablename="top_rapiditydiff_Marco";
  if(histname=="topCosTheta") observablename="top_costheta_cms";
  if(histname=="lepCosOpeningAngle") observablename="lep_cos_opening_angle";
  if(histname=="pseudorapiditydiff") observablename="top_pseudorapiditydiff_cms";
  if(histname=="rapiditydiff") observablename="top_rapiditydiff_cms";


  TString FName1 = "SIGoutput/ttdl_mcatnlo_smallTree_histos.root";

  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  

  hdiel_numerator = (TH1F*)f_1->Get(Form("h_numerator_%s_gen_diel", observablename.Data()));

  hmuel_numerator = (TH1F*)f_1->Get(Form("h_numerator_%s_gen_muel", observablename.Data())); 
  helmu_numerator = (TH1F*)f_1->Get(Form("h_numerator_%s_gen_elmu", observablename.Data())); 

  hdimu_numerator = (TH1F*)f_1->Get(Form("h_numerator_%s_gen_dimu", observablename.Data()));




  std::cout << "Opened " << Form("numerator_%s_gen_elmu", observablename.Data()) <<"\n";


  double KS42 = hdimu_numerator->KolmogorovTest(hmuel_numerator);
  double KS12 = hdiel_numerator->KolmogorovTest(hmuel_numerator);
  double KS32 = helmu_numerator->KolmogorovTest(hmuel_numerator);
  double KS41 = hdimu_numerator->KolmogorovTest(hdiel_numerator);

  cout<<"K-S mm,me: "<<KS42<<"; ee,me: "<<KS12<<"; em,me: "<<KS32<<"; mm,ee: "<<KS41<<endl;



  double Asym1,Asym2,Asym2b,Asym3;
  double Asym1err,Asym2err,Asym2berr,Asym3err;

  hdiel_numerator->Rebin(12);
  hmuel_numerator->Rebin(12);
  helmu_numerator->Rebin(12);
  hdimu_numerator->Rebin(12);

  TCanvas *c2 = new TCanvas();
  c2->cd();


  hmuel_numerator->SetLineColor(kBlue-7);
  helmu_numerator->SetLineColor(kBlue+2);
  hmuel_numerator-> SetFillColor(0);
  helmu_numerator-> SetFillColor(0);
  hmuel_numerator->SetMarkerColor(kBlue-7);
  helmu_numerator->SetMarkerColor(kBlue+2);
  hdimu_numerator->SetLineColor(kRed);
  hdimu_numerator->SetMarkerColor(kRed);
  hdimu_numerator-> SetFillColor(0);
  hdiel_numerator->SetLineColor(kBlack);
  hdiel_numerator->SetMarkerColor(kBlack);
  hdiel_numerator-> SetFillColor(0);

  hdiel_numerator->DrawNormalized("hist");   
  hmuel_numerator->DrawNormalized("histsame");
  helmu_numerator->DrawNormalized("histsame");
  hdimu_numerator->DrawNormalized("histsame");

      GetAfb(hdiel_numerator,Asym1, Asym1err);
      GetAfb(hmuel_numerator,Asym2, Asym2err);
      GetAfb(helmu_numerator,Asym2b, Asym2berr);
      GetAfb(hdimu_numerator,Asym3, Asym3err);
      cout<<"asym: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym2b<<" +/- "<<1*Asym2berr<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;






  TLegend *leg = new TLegend(0.47,0.18,0.69,0.27);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  leg->AddEntry(hdiel_numerator, "ee","l");
  if(drawnorm) leg->AddEntry(hmuel_numerator, "#mue","l");
  if(drawnorm) leg->AddEntry(helmu_numerator, "e#mu","l");
  if(drawnorm) leg->AddEntry(hdimu_numerator, "#mu#mu","l");
  
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

  TString Asym1_temp = formatFloat(Asym1,"%6.4f");  Asym1_temp.ReplaceAll(" " , "" );
  TString Asym1err_temp = formatFloat(Asym1err,"%6.4f");  Asym1err_temp.ReplaceAll(" " , "" );
  Asym1_temp = TString("   Asym: ") +  Asym1_temp + " #pm " + Asym1err_temp;
  blah = pt1->AddText(Asym1_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlack);  

  TString Asym2_temp = formatFloat(Asym2,"%6.4f");  Asym2_temp.ReplaceAll(" " , "" );
  TString Asym2err_temp = formatFloat(Asym2err,"%6.4f");  Asym2err_temp.ReplaceAll(" " , "" );
  Asym2_temp = TString("   Asym: ") +  Asym2_temp + " #pm " + Asym2err_temp;
  blah = pt1->AddText(Asym2_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlue-7);  

  TString Asym2b_temp = formatFloat(Asym2b,"%6.4f");  Asym2b_temp.ReplaceAll(" " , "" );
  TString Asym2berr_temp = formatFloat(Asym2berr,"%6.4f");  Asym2berr_temp.ReplaceAll(" " , "" );
  Asym2b_temp = TString("   Asym: ") +  Asym2b_temp + " #pm " + Asym2berr_temp;
  blah = pt1->AddText(Asym2b_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlue+2);  

  TString Asym3_temp = formatFloat(Asym3,"%6.4f");  Asym3_temp.ReplaceAll(" " , "" );
  TString Asym3err_temp = formatFloat(Asym3err,"%6.4f");  Asym3err_temp.ReplaceAll(" " , "" );
  Asym3_temp = TString("   Asym: ") +  Asym3_temp + " #pm " + Asym3err_temp;
  blah = pt1->AddText(Asym3_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kRed); 

  pt1->Draw(); 



  TPaveText *pt2 = new TPaveText(0.60, 0.77, 0.82, 0.92, "brNDC");
  pt2->SetName("pt2name");
  pt2->SetBorderSize(0);
  pt2->SetFillStyle(0);
  
  TText *blah2;

  blah2 = pt2->AddText("");

  TString KS42_temp = formatFloat(KS42,"%6.4f");
  KS42_temp.ReplaceAll(" " , "" );
  KS42_temp = TString("   KS #mu#mu,#mue: ") +  KS42_temp;
  blah2 = pt2->AddText(KS42_temp.Data());
  blah2->SetTextSize(0.032);
  blah2->SetTextAlign(11);  
  blah2->SetTextColor(kBlack);  

  TString KS12_temp = formatFloat(KS12,"%6.4f");
  KS12_temp.ReplaceAll(" " , "" );
  KS12_temp = TString("   KS ee,#mue: ") +  KS12_temp;
  blah2 = pt2->AddText(KS12_temp.Data());
  blah2->SetTextSize(0.032);
  blah2->SetTextAlign(11);  
  blah2->SetTextColor(kBlack);  

  TString KS32_temp = formatFloat(KS32,"%6.4f");
  KS32_temp.ReplaceAll(" " , "" );
  KS32_temp = TString("   KS e#mu,#mue: ") +  KS32_temp;
  blah2 = pt2->AddText(KS32_temp.Data());
  blah2->SetTextSize(0.032);
  blah2->SetTextAlign(11);  
  blah2->SetTextColor(kBlack);  

  TString KS41_temp = formatFloat(KS41,"%6.4f");
  KS41_temp.ReplaceAll(" " , "" );
  KS41_temp = TString("   KS #mu#mu,ee: ") +  KS41_temp;
  blah2 = pt2->AddText(KS41_temp.Data());
  blah2->SetTextSize(0.032);
  blah2->SetTextAlign(11);  
  blah2->SetTextColor(kBlack);  

  
  pt2->Draw();





  leg->Draw("same");

  c2->Print(Form("%s_numerator.pdf", histname.Data()));




  f_1->Close();
  
}
