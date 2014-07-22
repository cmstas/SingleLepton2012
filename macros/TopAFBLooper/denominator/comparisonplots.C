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

TH1F* hnumerator;
TH1F* hdenominator;
TH1F* hacceptance;
TH2F* hnumerator2d;
TH2F* hdenominator2d;
//TH2F* hacceptance2d;
TH2F* hnumerator2drebinned;
TH2F* hdenominator2drebinned;
TH2F* hacceptance2drebinned;

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

//void comparisonplots(TString histname = "topSpinCorrGen", bool drawnorm = true, TString FName1 = "denominator_powheg_stat1/results/hist_noCuts.root", TString FName2 = "denominator_powhegscfix_stat1/results/hist_noCuts.root", TString FName3 = "denominator_mcatnlo/results/hist_noCuts.root"){
void comparisonplots(TString histname = "topSpinCorrGen", bool drawnorm = true, TString FName1 = "denominator_powheg/results/hist_noCuts.root", TString FName2 = "denominator_powhegscfix/results/hist_noCuts.root", TString FName3 = "denominator_mcatnlo/results/hist_noCuts.root"){
  setTDRStyle();

  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  hnumerator = (TH1F*)f_1->Get(Form("ttdil_h%s_allj_all", histname.Data())); 
  hnumerator2d = (TH2F*)f_1->Get(Form("ttdil_h%s_allj_all", histname.Data())); 
    
  std::cout << "Opening " << FName2.Data() << "\n";  
  TFile *f_2         = TFile::Open(FName2.Data());  
  hdenominator = (TH1F*)f_2->Get(Form("ttdil_h%s_allj_all", histname.Data()));
  hdenominator2d = (TH2F*)f_2->Get(Form("ttdil_h%s_allj_all", histname.Data())); 
  
  std::cout << "Opening " << FName3.Data() << "\n";  
  TFile *f_3         = TFile::Open(FName3.Data());  
  hacceptance = (TH1F*)f_3->Get(Form("ttdil_h%s_allj_all", histname.Data()));
  //hdenominator2d = (TH2F*)f_3->Get(Form("ttdil_h%s_allj_all", histname.Data()));
  
  std::cout << "Opened " << Form("ttdil_h%s_allj_all", histname.Data()) << " and "<< Form("ttdil_h%s_allj_all", histname.Data()) <<"\n";

  //hacceptance->Print();
  
  Double_t pi = 3.141592653589793;


/*  
  if(histname.Contains("lepChargeAsym") ||  histname.Contains("rapiditydiff")) {

  	hnumerator = (TH1F*) hnumerator->Rebin(6,Form("numerator_%s", histname.Data()),bins1);
  	hdenominator = (TH1F*) hdenominator->Rebin(6,Form("denominator_%s", histname.Data()),bins1);

  	hnumerator2drebinned = new TH2F(Form("numerator_%s_mtt", histname.Data()),Form("numerator_%s_mtt", histname.Data()),2,bins1forMtt,3, binsMtt);
  	TAxis *xaxis = hnumerator2d->GetXaxis();
  	TAxis *yaxis = hnumerator2d->GetYaxis();
  	for (int j=1;j<=yaxis->GetNbins();j++) {
  		for (int i=1;i<=xaxis->GetNbins();i++) {
  			hnumerator2drebinned->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d->GetBinContent(i,j));
  		}
  	}
  	
  	hdenominator2drebinned = new TH2F(Form("denominator_%s_mtt", histname.Data()),Form("denominator_%s_mtt", histname.Data()),2,bins1forMtt,3, binsMtt);
  	TAxis *xaxisd = hdenominator2d->GetXaxis();
  	TAxis *yaxisd = hdenominator2d->GetYaxis();
  	for (int j=1;j<=yaxisd->GetNbins();j++) {
  		for (int i=1;i<=xaxisd->GetNbins();i++) {
  			hdenominator2drebinned->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d->GetBinContent(i,j));
  		}
  	}
  	
  }
  
  else  if(histname.Contains("lepAzimAsym2") ) {

  	hnumerator = (TH1F*) hnumerator->Rebin(6,Form("numerator_%s", histname.Data()),bins3);
  	hdenominator = (TH1F*) hdenominator->Rebin(6,Form("denominator_%s", histname.Data()),bins3);

  	hnumerator2drebinned = new TH2F(Form("numerator_%s_mtt", histname.Data()),Form("numerator_%s_mtt", histname.Data()),2,bins3forMtt,3, binsMtt);
  	TAxis *xaxis = hnumerator2d->GetXaxis();
  	TAxis *yaxis = hnumerator2d->GetYaxis();
  	for (int j=1;j<=yaxis->GetNbins();j++) {
  		for (int i=1;i<=xaxis->GetNbins();i++) {
  			hnumerator2drebinned->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d->GetBinContent(i,j));
  		}
  	}
  	
  	hdenominator2drebinned = new TH2F(Form("denominator_%s_mtt", histname.Data()),Form("denominator_%s_mtt", histname.Data()),2,bins3forMtt,3, binsMtt);
  	TAxis *xaxisd = hdenominator2d->GetXaxis();
  	TAxis *yaxisd = hdenominator2d->GetYaxis();
  	for (int j=1;j<=yaxisd->GetNbins();j++) {
  		for (int i=1;i<=xaxisd->GetNbins();i++) {
  			hdenominator2drebinned->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d->GetBinContent(i,j));
  		}
  	}
  	
  }




  
  if(true) {
 	
  	hnumerator = (TH1F*) hnumerator->Rebin(6,Form("numerator_%s", histname.Data()),bins2);
  	hdenominator = (TH1F*) hdenominator->Rebin(6,Form("denominator_%s", histname.Data()),bins2);
  	
    	hnumerator2drebinned = new TH2F(Form("numerator_%s_mtt", histname.Data()),Form("numerator_%s_mtt", histname.Data()),2,bins2forMtt,3, binsMtt);
  	TAxis *xaxis = hnumerator2d->GetXaxis();
  	TAxis *yaxis = hnumerator2d->GetYaxis();
  	for (int j=1;j<=yaxis->GetNbins();j++) {
  		for (int i=1;i<=xaxis->GetNbins();i++) {
  			hnumerator2drebinned->Fill(xaxis->GetBinCenter(i),yaxis->GetBinCenter(j), hnumerator2d->GetBinContent(i,j));
  		}
  	}
  	
  	hdenominator2drebinned = new TH2F(Form("denominator_%s_mtt", histname.Data()),Form("denominator_%s_mtt", histname.Data()),2,bins2forMtt,3, binsMtt);
  	TAxis *xaxisd = hdenominator2d->GetXaxis();
  	TAxis *yaxisd = hdenominator2d->GetYaxis();
  	for (int j=1;j<=yaxisd->GetNbins();j++) {
  		for (int i=1;i<=xaxisd->GetNbins();i++) {
  			hdenominator2drebinned->Fill(xaxisd->GetBinCenter(i),yaxisd->GetBinCenter(j), hdenominator2d->GetBinContent(i,j));
  		}
  	}
    	
  }
  */
  



  double Asym1,Asym2,Asym3;
  double Asym1err,Asym2err,Asym3err;

    if(histname.Contains("topSpinCorr")){
      Asym1 = hacceptance->GetMean();
      Asym2 = hnumerator->GetMean();
      Asym3 = hdenominator->GetMean();
      Asym1err = hacceptance->GetMeanError();
      Asym2err = hnumerator->GetMeanError();
      Asym3err = hdenominator->GetMeanError();
      cout<<"mean: "<<9*Asym1<<" +/- "<<9*Asym1err<<", "<<9*Asym2<<" +/- "<<9*Asym2err<<", "<<9*Asym3<<" +/- "<<9*Asym3err<<endl;

      GetAfb(hacceptance,Asym1, Asym1err);
      GetAfb(hnumerator,Asym2, Asym2err);
      GetAfb(hdenominator,Asym3, Asym3err);
      cout<<"asym: "<<4*Asym1<<" +/- "<<4*Asym1err<<", "<<4*Asym2<<" +/- "<<4*Asym2err<<", "<<4*Asym3<<" +/- "<<4*Asym3err<<endl;
    }
    else if(histname.Contains("Cos")){
      Asym1 = hacceptance->GetMean();
      Asym2 = hnumerator->GetMean();
      Asym3 = hdenominator->GetMean();
      Asym1err = hacceptance->GetMeanError();
      Asym2err = hnumerator->GetMeanError();
      Asym3err = hdenominator->GetMeanError();
      cout<<"mean: "<<3*Asym1<<" +/- "<<3*Asym1err<<", "<<3*Asym2<<" +/- "<<3*Asym2err<<", "<<3*Asym3<<" +/- "<<3*Asym3err<<endl;

      GetAfb(hacceptance,Asym1, Asym1err);
      GetAfb(hnumerator,Asym2, Asym2err);
      GetAfb(hdenominator,Asym3, Asym3err);
      cout<<"asym: "<<2*Asym1<<" +/- "<<2*Asym1err<<", "<<2*Asym2<<" +/- "<<2*Asym2err<<", "<<2*Asym3<<" +/- "<<2*Asym3err<<endl;
    }
    else {
      Asym1 = hacceptance->GetMean();
      Asym2 = hnumerator->GetMean();
      Asym3 = hdenominator->GetMean();
      Asym1err = hacceptance->GetMeanError();
      Asym2err = hnumerator->GetMeanError();
      Asym3err = hdenominator->GetMeanError();
      cout<<"mean: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;

      GetAfb(hacceptance,Asym1, Asym1err);
      GetAfb(hnumerator,Asym2, Asym2err);
      GetAfb(hdenominator,Asym3, Asym3err);
      cout<<"asym: "<<1*Asym1<<" +/- "<<1*Asym1err<<", "<<1*Asym2<<" +/- "<<1*Asym2err<<", "<<1*Asym3<<" +/- "<<1*Asym3err<<endl;
    }


  hnumerator->Rebin(8);
  hdenominator->Rebin(8);
  hacceptance->Rebin(8);
  
  TString accepthistname = "SCfixcompare_";
  accepthistname += histname;
  
  //TFile *output = new TFile(Form("%s.root", accepthistname.Data()), "RECREATE");  
  
  //hacceptance =  (TH1F*) hnumerator->Clone(accepthistname.Data());
  //hacceptance->SetTitle(accepthistname.Data());
  //hacceptance->Reset();
  //hacceptance->Divide(hnumerator,hdenominator,1., 1.);
  
  //hacceptance2d =  (TH2F*) hnumerator2d->Clone( Form("%s_mtt_notrebinned", accepthistname.Data()) );
  //hacceptance2d->Reset();
  //hacceptance2d->SetTitle(Form("%s_mtt_notrebinned", accepthistname.Data()));
  //hacceptance2d->Divide(hnumerator2d,hdenominator2d,1., 1.);
  
  //hacceptance2drebinned =  (TH2F*) hnumerator2drebinned->Clone( Form("%s_mtt", accepthistname.Data()) );
  //hacceptance2drebinned->Reset();
  //hacceptance2drebinned->SetTitle(Form("%s_mtt", accepthistname.Data()));
  //hacceptance2drebinned->Divide(hnumerator2drebinned,hdenominator2drebinned,1., 1.);
  
  
    
  hnumerator->SetLineColor(kBlue);
  hnumerator-> SetFillColor(0);
  hnumerator->SetMarkerColor(kBlue);
  hdenominator->SetLineColor(kRed);
  hdenominator->SetMarkerColor(kRed);
  hdenominator-> SetFillColor(0);
  hacceptance->SetLineColor(kBlack);
  hacceptance->SetMarkerColor(kBlack);
  hacceptance-> SetFillColor(0);
  
  gStyle->SetPaintTextFormat("6.4f");

  TCanvas *c1 = new TCanvas();
  c1->cd();
  
  hacceptance->SetMaximum(1.25*hacceptance->GetMaximum());
  if(hacceptance->GetMinimum() <0.15 *hacceptance->GetMaximum() ) hacceptance->SetMinimum(0.);  
  if(hacceptance->GetMinimum() > 0.) hacceptance->SetMinimum(0.75*hacceptance->GetMinimum() );  
  
  if(!drawnorm){
  	hacceptance->Draw("hist TEXT00E");
  }
  else{
  	//hacceptance->Scale(1.,"width");
  	//hnumerator->Scale(1.,"width");
  	//hdenominator->Scale(1.,"width");
  	  	
  	//TH1 *frame=new TH1F("frame","",1000,0.,3.1415926535);frame->SetMaximum(0.0024);frame->SetMinimum(0.12);frame->Draw();
         //TH1 *frame=new TH1F("frame","",1000,-1.,1.);frame->SetMaximum(0.45);frame->SetMinimum(0.);frame->Draw();
  	//hacceptance->SetMaximum(0.4);
  	//hacceptance->SetMinimum(0.); 
  	//hnumerator->SetMaximum(0.4);
  	//hnumerator->SetMinimum(0.);
  	//hdenominator->SetMaximum(0.4);
  	//hdenominator->SetMinimum(0.); 
  	//hacceptance->SetMinimum(0.5);
  	//hacceptance->SetMaximum(1.5);
  	//hacceptance->GetXaxis()->SetTitle("");
  	hacceptance->DrawNormalized("hist");  	
  	hnumerator->DrawNormalized("histsame");
  	hdenominator->DrawNormalized("histsame");
  }

  TLegend *leg = new TLegend(0.47,0.18,0.69,0.27);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.032);
  leg->SetFillStyle(0);
  leg->AddEntry(hacceptance, "mc@nlo","l");
  //if(drawnorm) leg->AddEntry(hnumerator, "powheg (stat1 leps)","l");
  //if(drawnorm) leg->AddEntry(hdenominator, "powheg SC fix (stat1 leps)","l");
  if(drawnorm) leg->AddEntry(hnumerator, "powheg","l");
  if(drawnorm) leg->AddEntry(hdenominator, "powheg SC fix","l");
  
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

  TString Asym1_temp = formatFloat(Asym1,"%6.4f");
  Asym1_temp.ReplaceAll(" " , "" );
  Asym1_temp = TString("   Asym: ") +  Asym1_temp;
  blah = pt1->AddText(Asym1_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlack);  

  TString Asym2_temp = formatFloat(Asym2,"%6.4f");
  Asym2_temp.ReplaceAll(" " , "" );
  Asym2_temp = TString("   Asym: ") +  Asym2_temp;
  blah = pt1->AddText(Asym2_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kBlue);  

  TString Asym3_temp = formatFloat(Asym3,"%6.4f");
  Asym3_temp.ReplaceAll(" " , "" );
  Asym3_temp = TString("   Asym: ") +  Asym3_temp;
  blah = pt1->AddText(Asym3_temp.Data());
  blah->SetTextSize(0.032);
  blah->SetTextAlign(11);  
  blah->SetTextColor(kRed);  


  
  pt1->Draw();


  c1->Print(Form("%s.pdf", accepthistname.Data()));
  //hacceptance->Write();
  //hdenominator->Write();
  //hnumerator->Write();
  
  //hacceptance2d->Write();
  //hdenominator2d->Write();
  //hnumerator2d->Write();

  //hacceptance2drebinned->Write();
  //hdenominator2drebinned->Write();
  //hnumerator2drebinned->Write();

  f_1->Close();
  f_2->Close();
  f_3->Close();
  //output->Close();
  
}
