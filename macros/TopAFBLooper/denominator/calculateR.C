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



TH1F* hmueg_numerator_1;
TH1F* hdimu_numerator_1;
TH1F* hdiel_numerator_1;
TH1F* hall_numerator_1;
TH1F* hmueg_denominator_1;
TH1F* hdimu_denominator_1;
TH1F* hdiel_denominator_1;
TH1F* hall_denominator_1;

TH1F* hmueg_numerator_2;
TH1F* hdimu_numerator_2;
TH1F* hdiel_numerator_2;
TH1F* hall_numerator_2;
TH1F* hmueg_denominator_2;
TH1F* hdimu_denominator_2;
TH1F* hdiel_denominator_2;
TH1F* hall_denominator_2;





TH2F* hmueg2d;
TH2F* hdimu2d;
//TH2F* hdiel2d;
TH2F* hmueg2drebinned;
TH2F* hdimu2drebinned;
TH2F* hdiel2drebinned;

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

void calculateR(TString histname = "topSpinCorr", bool drawnorm = true) {
  setTDRStyle();
  TH1::SetDefaultSumw2();

  TString FName1 = Form("acceptance_default/mcnlo/accept_%s.root", histname.Data());
  TString FName2 = Form("acceptance_nospincorrelation/mcnlo/accept_%s.root", histname.Data());


  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  

  hdiel_numerator_1 = (TH1F*)f_1->Get(Form("numerator_%s_diel", histname.Data()));
  double hdiel_numerator_integral_1 = hdiel_numerator_1->Integral();

  hmueg_numerator_1 = (TH1F*)f_1->Get(Form("numerator_%s_mueg", histname.Data())); 
  double hmueg_numerator_integral_1 = hmueg_numerator_1->Integral();

  hdimu_numerator_1 = (TH1F*)f_1->Get(Form("numerator_%s_dimu", histname.Data()));
  double hdimu_numerator_integral_1 = hdimu_numerator_1->Integral();

  hall_numerator_1 = (TH1F*)f_1->Get(Form("numerator_%s_all", histname.Data()));
  double hall_numerator_integral_1 = hall_numerator_1->Integral();

  hdiel_denominator_1 = (TH1F*)f_1->Get(Form("denominator_%s_diel", histname.Data()));
  double hdiel_denominator_integral_1 = hdiel_denominator_1->Integral();
  double hdiel_ratio_1 = hdiel_numerator_integral_1/hdiel_denominator_integral_1;

  hmueg_denominator_1 = (TH1F*)f_1->Get(Form("denominator_%s_mueg", histname.Data())); 
  double hmueg_denominator_integral_1 = hmueg_denominator_1->Integral();
  double hmueg_ratio_1 = hmueg_numerator_integral_1/hmueg_denominator_integral_1;

  hdimu_denominator_1 = (TH1F*)f_1->Get(Form("denominator_%s_dimu", histname.Data()));
  double hdimu_denominator_integral_1 = hdimu_denominator_1->Integral();
  double hdimu_ratio_1 = hdimu_numerator_integral_1/hdimu_denominator_integral_1;

  hall_denominator_1 = (TH1F*)f_1->Get(Form("denominator_%s_all", histname.Data()));
  double hall_denominator_integral_1 = hall_denominator_1->Integral();
  double hall_ratio_1 = hall_numerator_integral_1/hall_denominator_integral_1;


  std::cout << "Opening " << FName2.Data() << "\n";
  TFile *f_2         = TFile::Open(FName2.Data());  

  hdiel_numerator_2 = (TH1F*)f_2->Get(Form("numerator_%s_diel", histname.Data()));
  double hdiel_numerator_integral_2 = hdiel_numerator_2->Integral();

  hmueg_numerator_2 = (TH1F*)f_2->Get(Form("numerator_%s_mueg", histname.Data())); 
  double hmueg_numerator_integral_2 = hmueg_numerator_2->Integral();

  hdimu_numerator_2 = (TH1F*)f_2->Get(Form("numerator_%s_dimu", histname.Data()));
  double hdimu_numerator_integral_2 = hdimu_numerator_2->Integral();

  hall_numerator_2 = (TH1F*)f_2->Get(Form("numerator_%s_all", histname.Data()));
  double hall_numerator_integral_2 = hall_numerator_2->Integral();

  hdiel_denominator_2 = (TH1F*)f_2->Get(Form("denominator_%s_diel", histname.Data()));
  double hdiel_denominator_integral_2 = hdiel_denominator_2->Integral();
  double hdiel_ratio_2 = hdiel_numerator_integral_2/hdiel_denominator_integral_2;

  hmueg_denominator_2 = (TH1F*)f_2->Get(Form("denominator_%s_mueg", histname.Data())); 
  double hmueg_denominator_integral_2 = hmueg_denominator_2->Integral();
  double hmueg_ratio_2 = hmueg_numerator_integral_2/hmueg_denominator_integral_2;

  hdimu_denominator_2 = (TH1F*)f_2->Get(Form("denominator_%s_dimu", histname.Data()));
  double hdimu_denominator_integral_2 = hdimu_denominator_2->Integral();
  double hdimu_ratio_2 = hdimu_numerator_integral_2/hdimu_denominator_integral_2;

  hall_denominator_2 = (TH1F*)f_2->Get(Form("denominator_%s_all", histname.Data()));
  double hall_denominator_integral_2 = hall_denominator_2->Integral();
  double hall_ratio_2 = hall_numerator_integral_2/hall_denominator_integral_2;

  cout<<hall_ratio_2/hall_ratio_1<<" "<<hdiel_ratio_2/hdiel_ratio_1<<" "<<hmueg_ratio_2/hmueg_ratio_1<<" "<<hdimu_ratio_2/hdimu_ratio_1<<endl;


  TH1D* hall_acceptance_1 =  (TH1D*) hall_numerator_1->Clone( "hall_acceptance_1" );
  hall_acceptance_1->Divide(hall_denominator_1);

  TH1D* hall_acceptance_2 =  (TH1D*) hall_numerator_2->Clone( "hall_acceptance_2" );
  hall_acceptance_2->Divide(hall_denominator_2);
  hall_acceptance_2->Divide(hall_acceptance_1);


  TH1D* hdiel_acceptance_1 =  (TH1D*) hdiel_numerator_1->Clone( "hdiel_acceptance_1" );
  hdiel_acceptance_1->Divide(hdiel_denominator_1);

  TH1D* hdiel_acceptance_2 =  (TH1D*) hdiel_numerator_2->Clone( "hdiel_acceptance_2" );
  hdiel_acceptance_2->Divide(hdiel_denominator_2);
  hdiel_acceptance_2->Divide(hdiel_acceptance_1);


  TH1D* hmueg_acceptance_1 =  (TH1D*) hmueg_numerator_1->Clone( "hmueg_acceptance_1" );
  hmueg_acceptance_1->Divide(hmueg_denominator_1);

  TH1D* hmueg_acceptance_2 =  (TH1D*) hmueg_numerator_2->Clone( "hmueg_acceptance_2" );
  hmueg_acceptance_2->Divide(hmueg_denominator_2);
  hmueg_acceptance_2->Divide(hmueg_acceptance_1);


  TH1D* hdimu_acceptance_1 =  (TH1D*) hdimu_numerator_1->Clone( "hdimu_acceptance_1" );
  hdimu_acceptance_1->Divide(hdimu_denominator_1);

  TH1D* hdimu_acceptance_2 =  (TH1D*) hdimu_numerator_2->Clone( "hdimu_acceptance_2" );
  hdimu_acceptance_2->Divide(hdimu_denominator_2);
  hdimu_acceptance_2->Divide(hdimu_acceptance_1);



  TCanvas *c1 = new TCanvas();
  c1->cd();

  hmueg_acceptance_2->SetLineColor(kBlue);
  hmueg_acceptance_2-> SetFillColor(0);
  //hmueg_acceptance_2->SetMinimum(0.998);
  //hmueg_acceptance_2->SetMaximum(1.002);
  hmueg_acceptance_2->SetMarkerColor(kBlue);
  hdimu_acceptance_2->SetLineColor(kRed);
  hdimu_acceptance_2->SetMarkerColor(kRed);
  hdimu_acceptance_2-> SetFillColor(0);
  //hdimu_acceptance_2->SetMinimum(0.998);
  //hdimu_acceptance_2->SetMaximum(1.002);
  hall_acceptance_2->SetLineColor(kBlack);
  hall_acceptance_2->SetMarkerColor(kBlack);
  hall_acceptance_2-> SetFillColor(0);
  hall_acceptance_2->SetMinimum(0.9);
  hall_acceptance_2->SetMaximum(1.1);
  hall_acceptance_2->Draw("hist");   
  hmueg_acceptance_2->Draw("histsame");
  //hdimu_acceptance_2->Draw("histsame");  

  c1->Print(Form("%s_uncorr_corr_acc_ratio.pdf", histname.Data()));


/*


  double Asym1,Asym2,Asym3;
  double Asym1err,Asym2err,Asym3err;

      Asym1 = hdiel->GetMean();
      Asym2 = hmueg->GetMean();
      Asym3 = hdimu->GetMean();
      Asym1err = hdiel->GetMeanError();
      Asym2err = hmueg->GetMeanError();
      Asym3err = hdimu->GetMeanError();
      cout<<"mean: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;


      GetAfb(hdiel_denominator,Asym1, Asym1err);
      GetAfb(hmueg_denominator,Asym2, Asym2err);
      GetAfb(hdimu_denominator,Asym3, Asym3err);
      cout<<"asym_denom: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;

      GetAfb(h_numerator_acccorr_diel,Asym1, Asym1err);
      GetAfb(h_numerator_acccorr_mueg,Asym2, Asym2err);
      GetAfb(h_numerator_acccorr_dimu,Asym3, Asym3err);
      cout<<"asym_acorr: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;


      GetAfb(hdiel,Asym1, Asym1err);
      GetAfb(hmueg,Asym2, Asym2err);
      GetAfb(hdimu,Asym3, Asym3err);
      cout<<"asym: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;


      double AsymD,AsymC;
      double AsymDerr,AsymCerr;
      GetAfb(hall_denominator,AsymD, AsymDerr);
      GetAfb(h_numerator_acccorr_all,AsymC, AsymCerr);
      cout<<"combined asym test: "<<1*AsymD<<" +/- "<<1*AsymDerr<<", "<<1*AsymC<<" +/- "<<1*AsymCerr<<endl;



  double KS32 = hdimu->KolmogorovTest(hmueg);
  double KS12 = hdiel->KolmogorovTest(hmueg);
  double KS31 = hdimu->KolmogorovTest(hdiel);

  cout<<"K-S mm,em: "<<KS32<<"; ee,em: "<<KS12<<"; mm,ee: "<<KS31<<endl;

  
  TString accepthistname = "compare_channels_acceptance_";
  accepthistname += histname;
  
  
    
  hmueg->SetLineColor(kBlue);
  hmueg-> SetFillColor(0);
  hmueg->SetMarkerColor(kBlue);
  hdimu->SetLineColor(kRed);
  hdimu->SetMarkerColor(kRed);
  hdimu-> SetFillColor(0);
  hdiel->SetLineColor(kBlack);
  hdiel->SetMarkerColor(kBlack);
  hdiel-> SetFillColor(0);
  
  gStyle->SetPaintTextFormat("6.4f");

  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  hdiel->SetMaximum(1.25*hdiel->GetMaximum());
  if(hdiel->GetMinimum() <0.15 *hdiel->GetMaximum() ) hdiel->SetMinimum(0.);  
  if(hdiel->GetMinimum() > 0.) hdiel->SetMinimum(0.75*hdiel->GetMinimum() );  
  
  if(!drawnorm){
  	hdiel->Draw("hist TEXT00E");
  }
  else{
  	//hdiel->Scale(1.,"width");
  	//hmueg->Scale(1.,"width");
  	//hdimu->Scale(1.,"width");
  	  	
  	//TH1 *frame=new TH1F("frame","",1000,0.,3.1415926535);frame->SetMaximum(0.0024);frame->SetMinimum(0.12);frame->Draw();
         //TH1 *frame=new TH1F("frame","",1000,-1.,1.);frame->SetMaximum(0.45);frame->SetMinimum(0.);frame->Draw();
  	//hdiel->SetMaximum(0.4);
  	//hdiel->SetMinimum(0.); 
  	//hmueg->SetMaximum(0.4);
  	//hmueg->SetMinimum(0.);
  	//hdimu->SetMaximum(0.4);
  	//hdimu->SetMinimum(0.); 
  	//hdiel->SetMinimum(0.5);
  	//hdiel->SetMaximum(1.5);
  	//hdiel->GetXaxis()->SetTitle("");
  	hdiel->DrawNormalized("histE");  	
  	hmueg->DrawNormalized("histEsame");
  	hdimu->DrawNormalized("histEsame");

    //hdiel_acceptance->DrawNormalized("histE");   
    //hmueg_acceptance->DrawNormalized("histEsame");
    //hdimu_acceptance->DrawNormalized("histEsame");
  }

  TLegend *leg = new TLegend(0.47,0.18,0.69,0.27);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  leg->AddEntry(hdiel, "ee","l");
  //if(drawnorm) leg->AddEntry(hmueg, "powheg (stat1 leps)","l");
  //if(drawnorm) leg->AddEntry(hdimu, "powheg SC fix (stat1 leps)","l");
  if(drawnorm) leg->AddEntry(hmueg, "e#mu","l");
  if(drawnorm) leg->AddEntry(hdimu, "#mu#mu","l");
  
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
  blah->SetTextColor(kBlue);  

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

  TString KS32_temp = formatFloat(KS32,"%6.4f");
  KS32_temp.ReplaceAll(" " , "" );
  KS32_temp = TString("   KS #mu#mu,e#mu: ") +  KS32_temp;
  blah2 = pt2->AddText(KS32_temp.Data());
  blah2->SetTextSize(0.032);
  blah2->SetTextAlign(11);  
  blah2->SetTextColor(kBlack);  

  TString KS12_temp = formatFloat(KS12,"%6.4f");
  KS12_temp.ReplaceAll(" " , "" );
  KS12_temp = TString("   KS ee,e#mu: ") +  KS12_temp;
  blah2 = pt2->AddText(KS12_temp.Data());
  blah2->SetTextSize(0.032);
  blah2->SetTextAlign(11);  
  blah2->SetTextColor(kBlack);  

  TString KS31_temp = formatFloat(KS31,"%6.4f");
  KS31_temp.ReplaceAll(" " , "" );
  KS31_temp = TString("   KS #mu#mu,ee: ") +  KS31_temp;
  blah2 = pt2->AddText(KS31_temp.Data());
  blah2->SetTextSize(0.032);
  blah2->SetTextAlign(11);  
  blah2->SetTextColor(kBlack);  

  
  pt2->Draw();


  c1->Print(Form("%s_gen.pdf", accepthistname.Data()));

  TCanvas *c2 = new TCanvas();
  c2->cd();

  hdiel_acceptance->Divide(hdiel_acceptance,hdiel,1., 1.);
  hmueg_acceptance->Divide(hmueg_acceptance,hmueg,1., 1.);
  hdimu_acceptance->Divide(hdimu_acceptance,hdimu,1., 1.);

  hmueg_acceptance->SetLineColor(kBlue);
  hmueg_acceptance-> SetFillColor(0);
  hmueg_acceptance->SetMinimum(0.998);
  hmueg_acceptance->SetMaximum(1.002);
  hmueg_acceptance->SetMarkerColor(kBlue);
  hdimu_acceptance->SetLineColor(kRed);
  hdimu_acceptance->SetMarkerColor(kRed);
  hdimu_acceptance-> SetFillColor(0);
  hdimu_acceptance->SetMinimum(0.998);
  hdimu_acceptance->SetMaximum(1.002);
  hdiel_acceptance->SetLineColor(kBlack);
  hdiel_acceptance->SetMarkerColor(kBlack);
  hdiel_acceptance-> SetFillColor(0);
  hdiel_acceptance->SetMinimum(0.998);
  hdiel_acceptance->SetMaximum(1.002);

  hdiel_acceptance->Draw("hist");   
  hmueg_acceptance->Draw("histsame");
  hdimu_acceptance->Draw("histsame");

      GetAfb(hdiel_acceptance,Asym1, Asym1err);
      GetAfb(hmueg_acceptance,Asym2, Asym2err);
      GetAfb(hdimu_acceptance,Asym3, Asym3err);
      cout<<"asym: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;

  TPaveText *pt1rat = new TPaveText(0.18, 0.77, 0.40, 0.92, "brNDC");
  pt1rat->SetName("pt1ratname");
  pt1rat->SetBorderSize(0);
  pt1rat->SetFillStyle(0);
  
  TText *blahrat;
  blahrat = pt1rat->AddText("CMS Simulation, #sqrt{s}=8 TeV");
  blahrat->SetTextSize(0.032);
  blahrat->SetTextAlign(11);  

  blahrat = pt1rat->AddText("");

  Asym1_temp = formatFloat(Asym1,"%6.4f");  Asym1_temp.ReplaceAll(" " , "" );
  Asym1err_temp = formatFloat(Asym1err,"%6.4f");  Asym1err_temp.ReplaceAll(" " , "" );
  Asym1_temp = TString("   Asym: ") +  Asym1_temp + " #pm " + Asym1err_temp;
  blahrat = pt1rat->AddText(Asym1_temp.Data());
  blahrat->SetTextSize(0.032);
  blahrat->SetTextAlign(11);  
  blahrat->SetTextColor(kBlack);  

  Asym2_temp = formatFloat(Asym2,"%6.4f");  Asym2_temp.ReplaceAll(" " , "" );
  Asym2err_temp = formatFloat(Asym2err,"%6.4f");  Asym2err_temp.ReplaceAll(" " , "" );
  Asym2_temp = TString("   Asym: ") +  Asym2_temp + " #pm " + Asym2err_temp;
  blahrat = pt1rat->AddText(Asym2_temp.Data());
  blahrat->SetTextSize(0.032);
  blahrat->SetTextAlign(11);  
  blahrat->SetTextColor(kBlue);  

  Asym3_temp = formatFloat(Asym3,"%6.4f");  Asym3_temp.ReplaceAll(" " , "" );
  Asym3err_temp = formatFloat(Asym3err,"%6.4f");  Asym3err_temp.ReplaceAll(" " , "" );
  Asym3_temp = TString("   Asym: ") +  Asym3_temp + " #pm " + Asym3err_temp;
  blahrat = pt1rat->AddText(Asym3_temp.Data());
  blahrat->SetTextSize(0.032);
  blahrat->SetTextAlign(11);  
  blahrat->SetTextColor(kRed); 

  pt1rat->Draw(); 





  leg->Draw("same");

  c2->Print(Form("%s_gen_ratio.pdf", accepthistname.Data()));


  //hdiel->Write();
  //hdimu->Write();
  //hmueg->Write();
  
  //hdiel2d->Write();
  //hdimu2d->Write();
  //hmueg2d->Write();

  //hdiel2drebinned->Write();
  //hdimu2drebinned->Write();
  //hmueg2drebinned->Write();
*/
  f_1->Close();
  f_2->Close();
  //f_3->Close();
  //output->Close();
  
}
