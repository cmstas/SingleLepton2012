#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "TCanvas.h"
#include <iostream>
#include "TFile.h"
#include "TDirectory.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "THStack.h"
#include "TLatex.h"
#include "TMath.h"
#include <fstream>
#include "tdrStyle.C"
#include "CommonFunctions.C"
#include "acceptanceplots.h"
#include <vector>

using namespace std;

TH1D* hnumerator;
TH1D* hdenominator;
TH1D* hacceptance;
TH1D* hacceptance_copy;
TH1D* hacceptance_statup;
TH1D* hacceptance_statdown;
TH1D* hnumerator_finebins;
TH1D* hdenominator_finebins[41];
TH1D* hdiff_finebins[41];


TH2D* hnumerator2d_mtt;
TH2D* hdenominator2d_mtt;
TH2D* hnumerator2drebinned_mtt;
TH2D* hdenominator2drebinned_mtt;
TH2D* hacceptance2drebinned_mtt;

TH2D* hnumerator2d_ttpt;
TH2D* hdenominator2d_ttpt;
TH2D* hnumerator2drebinned_ttpt;
TH2D* hdenominator2drebinned_ttpt;
TH2D* hacceptance2drebinned_ttpt;

TH2D* hnumerator2d_ttrapidity2;
TH2D* hdenominator2d_ttrapidity2;
TH2D* hnumerator2drebinned_ttrapidity2;
TH2D* hdenominator2drebinned_ttrapidity2;
TH2D* hacceptance2drebinned_ttrapidity2;


void GetAfb(TH1D* h, Double_t &afb, Double_t  &afberr){
 
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



void denominatorPDFvars(TString histname = "lepAzimAsym2", bool drawnorm = false, TString FName2 = "results/hist_noCuts.root"){

  int nbinsx = -999;
  if(histname == "lepChargeAsym" || histname == "lepAzimAsym2" || histname == "lepAzimAsym") nbinsx = 12;
  else nbinsx = 6;
  if(histname == "ttMass" || histname == "ttpT" || histname == "ttrap") nbinsx = 3;

  Double_t bins[nbinsx+1];

  if(histname == "lepChargeAsym") memcpy(bins,bins_lepChargeAsym,13*8);
  //if(histname == "lepChargeAsym") memcpy(bins,bins_lepChargeAsym,7*8);
  if(histname == "lepAzimAsym2") memcpy(bins,bins_lepAzimAsym2,13*8);
  //if(histname == "lepAzimAsym2") memcpy(bins,bins_lepAzimAsym2,7*8);
  if(histname == "lepAzimAsym") memcpy(bins,bins_lepAzimAsym,13*8);
  //if(histname == "lepAzimAsym") memcpy(bins,bins_lepAzimAsym,7*8);
  if(histname == "topCosTheta") memcpy(bins,bins_topCosTheta,7*8);
  if(histname == "pseudorapiditydiff") memcpy(bins,bins_pseudorapiditydiff,7*8);
  if(histname == "rapiditydiff") memcpy(bins,bins_rapiditydiff,7*8);
  if(histname == "rapiditydiffMarco") memcpy(bins,bins_rapiditydiffMarco,7*8);
  if(histname == "lepCosTheta" || histname == "lepPlusCosTheta" || histname == "lepMinusCosTheta") memcpy(bins,bins_lepCosTheta,7*8);
  if(histname == "topSpinCorr") memcpy(bins,bins_topSpinCorr,7*8);
  if(histname == "lepCosOpeningAngle") memcpy(bins,bins_lepCosOpeningAngle,7*8);

  //if(histname == "ttMass") memcpy(bins,ybinsmtt,4*8);
  //if(histname == "ttpT") memcpy(bins,ybinsttpt,4*8);
  //if(histname == "ttrap") memcpy(bins,ybinsttrapidity2,4*8);
 

/*
  if(histname == "lepChargeAsym") memcpy(binsfor2D,bins_lepChargeAsym_for2D,3*8);
  if(histname == "lepAzimAsym") memcpy(binsfor2D,bins_lepAzimAsym_for2D,3*8);
  if(histname == "lepAzimAsym2") memcpy(binsfor2D,bins_lepAzimAsym2_for2D,3*8);
  if(histname == "topCosTheta") memcpy(binsfor2D,bins_topCosTheta_for2D,3*8);
  if(histname == "pseudorapiditydiff") memcpy(binsfor2D,bins_pseudorapiditydiff_for2D,3*8);
  if(histname == "rapiditydiff") memcpy(binsfor2D,bins_rapiditydiff_for2D,3*8);
  if(histname == "rapiditydiffMarco") memcpy(binsfor2D,bins_rapiditydiffMarco_for2D,3*8);
  if(histname == "lepCosTheta" || histname == "lepPlusCosTheta" || histname == "lepMinusCosTheta") memcpy(binsfor2D,bins_lepCosTheta_for2D,3*8);
  if(histname == "topSpinCorr") memcpy(binsfor2D,bins_topSpinCorr_for2D,3*8);
  if(histname == "lepCosOpeningAngle") memcpy(binsfor2D,bins_lepCosOpeningAngle_for2D,3*8);
*/


  const int nChannels = 4;

  char suffixdenominator[nChannels][4]  = {"ee", "mm", "em", "all"};
  char suffixnumerator[nChannels][8]  = {"diel", "dimu", "mueg", "all"};

  setTDRStyle();

  TString accepthistname = "denominatorPDFvars_";
  accepthistname += histname;


  for (int ic = 3; ic < nChannels; ++ic)
  {

  std::cout << "Opening " << FName2.Data() << "\n";  
  TFile *f_2         = TFile::Open(FName2.Data());  

  for (int PDFset = 0; PDFset < 41; ++PDFset)
  {

  if(histname != "ttrap") hdenominator_finebins[PDFset] = (TH1D*)f_2->Get(Form("ttdil_h%sGen_PDFset%i_%s", histname.Data(), PDFset, suffixdenominator[ic]));
  else hdenominator_finebins[PDFset] = (TH1D*) ((TH2D*)f_2->Get(Form("ttdil_hlepAzimAsym2ttRapidity2Gen2d_PDFset%i_%s", PDFset, suffixdenominator[ic])))->ProjectionY();

  std::cout << "Opened " << Form("ttdil_h%sGen_PDFset%i_%s", histname.Data(), PDFset, suffixdenominator[ic]) <<"\n";

  if(histname == "ttMass" || histname == "ttpT" || histname == "ttrap")  hdenominator_finebins[PDFset] = (TH1D*) hdenominator_finebins[PDFset]->Rebin(6,Form("denominator_%s_%s", histname.Data(), suffixnumerator[ic]));
  else hdenominator_finebins[PDFset] = (TH1D*) hdenominator_finebins[PDFset]->Rebin(nbinsx,Form("denominator_%s_%s", histname.Data(), suffixnumerator[ic]),bins);
  hdenominator_finebins[PDFset]->Scale(1./hdenominator_finebins[PDFset]->Integral(),"width");

  hdiff_finebins[PDFset] = (TH1D*) hdenominator_finebins[PDFset]->Clone(Form("diff_%s_%s", histname.Data(), suffixnumerator[ic]));
  if(PDFset>0) hdiff_finebins[PDFset] -> Divide(hdiff_finebins[0]);

/*
  if(PDFset%2==0) {
    hdenominator_finebins[PDFset]->SetLineColor(kRed);
    hdenominator_finebins[PDFset]->SetMarkerColor(kRed);
    hdenominator_finebins[PDFset]-> SetFillColor(0);
  }
  else {
    hdenominator_finebins[PDFset]->SetLineColor(kGreen+1);
    hdenominator_finebins[PDFset]->SetMarkerColor(kGreen+1);
    hdenominator_finebins[PDFset]-> SetFillColor(0);
  }
*/

    hdenominator_finebins[PDFset]->SetLineColor( 1+ ((PDFset+1)/2)%9 );
    hdenominator_finebins[PDFset]->SetMarkerColor( 1+ ((PDFset+1)/2)%9 );
    hdenominator_finebins[PDFset]-> SetFillColor(0);
    hdenominator_finebins[PDFset]->SetLineWidth(1);


    hdiff_finebins[PDFset]->SetLineColor( 1+ ((PDFset+1)/2)%9 );
    hdiff_finebins[PDFset]->SetMarkerColor( 1+ ((PDFset+1)/2)%9 );
    hdiff_finebins[PDFset]-> SetFillColor(0);
    hdiff_finebins[PDFset]->SetLineWidth(1);


  }

  hdiff_finebins[0] -> Divide(hdiff_finebins[0]);

  hdenominator_finebins[0]->SetLineColor(kBlack);
  hdenominator_finebins[0]->SetLineWidth(2);
  hdenominator_finebins[0]->SetMarkerColor(kBlack);
  hdenominator_finebins[0]-> SetFillColor(0);

  //hdenominator_finebins = (TH1D*) hdenominator_finebins->Clone(Form("denominator_finebins_%s_%s", histname.Data(), suffixnumerator[ic]));


  gStyle->SetPaintTextFormat("6.4f");

  TCanvas *c1 = new TCanvas(Form("%s_%s_canvas", accepthistname.Data(), suffixnumerator[ic]), Form("%s_%s_canvas", accepthistname.Data(), suffixnumerator[ic]), 500, 725);


        TPad *p1, *p2;
        TLine *line;

            p1 = new TPad("p1", "dist", 0.0, 0.31, 1, 1.);
            p2 = new TPad("p2", "diff", 0.0, 0.0, 1, 0.31);
            p1->Draw();
            p2->Draw();
            p1->cd();

            p1->SetBottomMargin(0.15);
            //p2->SetTopMargin(0.1);
            p2->SetBottomMargin(0.31);


  hdenominator_finebins[0]->SetMaximum(1.25*hdenominator_finebins[0]->GetMaximum());
  if(hdenominator_finebins[0]->GetMinimum() <0.15 *hdenominator_finebins[0]->GetMaximum() ) hdenominator_finebins[0]->SetMinimum(0.);  
  if(hdenominator_finebins[0]->GetMinimum() > 0.) hdenominator_finebins[0]->SetMinimum(0.75*hdenominator_finebins[0]->GetMinimum() );  

  hdenominator_finebins[0]->GetYaxis()->SetTitle("Number of events");
  hdenominator_finebins[0]->GetYaxis()->SetDecimals(kTRUE);
  hdenominator_finebins[0]->SetTitleSize(0.06, "XYZ");
  hdenominator_finebins[0]->SetLabelSize(0.05, "XYZ");
  //if(!histname.Contains("lepAzimAsym2") ) hdenominator_finebins[0]->GetXaxis()->SetNdivisions(504,0);
  //else hdenominator_finebins[0]->GetXaxis()->SetNdivisions(506);
  hdenominator_finebins[0]->GetYaxis()->SetTitleOffset(1.3);


  if(histname.Contains("lepChargeAsym") ) {
    hdenominator_finebins[0]->GetXaxis()->SetTitle("#Delta|#eta_{l}|");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(#Delta|#eta_{l}|)");
    hdiff_finebins[0]->GetXaxis()->SetTitle("#Delta|#eta_{l}|");
  }
  if(histname.Contains("lepAzimAsym2") ) {
    //hdenominator_finebins[0]->GetXaxis()->SetTitle("#Delta#phi_{l+l-}");
    hdenominator_finebins[0]->GetXaxis()->SetTitle("#Delta#phi_{l#lower[-0.4]{+}l#lower[-0.48]{-}}");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(#Delta#phi_{l#lower[-0.4]{+}l#lower[-0.48]{-}})");
    hdiff_finebins[0]->GetXaxis()->SetTitle("#Delta#phi_{l#lower[-0.4]{+}l#lower[-0.48]{-}}");
  }
  if(histname.Contains("lepCosTheta") ) {
    //hdenominator_finebins[0]->GetXaxis()->SetTitle("cos(#theta_{l})");
    hdenominator_finebins[0]->GetXaxis()->SetTitle("cos(^{}#theta_{l}#kern[-0.35]{*})");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(cos(^{}#theta_{l}#kern[-0.35]{*}))");
    hdiff_finebins[0]->GetXaxis()->SetTitle("cos(^{}#theta_{l}#kern[-0.35]{*})");
  }
  if(histname.Contains("lepPlusCosTheta") ) {
    hdenominator_finebins[0]->GetXaxis()->SetTitle("cos(#theta_{l+})");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(cos(#theta_{l+}))");
    hdiff_finebins[0]->GetXaxis()->SetTitle("cos(#theta_{l+})");
  }
  if(histname.Contains("lepMinusCosTheta") ) {
    hdenominator_finebins[0]->GetXaxis()->SetTitle("cos(#theta_{l-})");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(cos(#theta_{l-}))");
    hdiff_finebins[0]->GetXaxis()->SetTitle("cos(#theta_{l-})");
  }
  if(histname.Contains("topSpinCorr") ) {
    //hdenominator_finebins[0]->GetXaxis()->SetTitle("cos(#theta_{l+})cos(#theta_{l-})");
    hdenominator_finebins[0]->GetXaxis()->SetTitle("cos(^{}#theta_{l#lower[-0.4]{+}}#kern[-1.38]{*}) cos(^{}#theta_{l#lower[-0.48]{-}}#kern[-1.0]{*})");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(cos(^{}#theta_{l#lower[-0.4]{+}}#kern[-1.38]{*}) cos(^{}#theta_{l#lower[-0.48]{-}}#kern[-1.0]{*}))");
    hdiff_finebins[0]->GetXaxis()->SetTitle("cos(^{}#theta_{l#lower[-0.4]{+}}#kern[-1.38]{*}) cos(^{}#theta_{l#lower[-0.48]{-}}#kern[-1.0]{*})");
  }
  if(histname.Contains("lepCosOpeningAngle") ) {
    hdenominator_finebins[0]->GetXaxis()->SetTitle("cos(#phi)");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(cos(#phi))");
    hdiff_finebins[0]->GetXaxis()->SetTitle("cos(#phi)");
  }
  if(histname.Contains("rapiditydiffMarco") ) {
    hdenominator_finebins[0]->GetXaxis()->SetTitle("#Delta|y_{t}|");
    hdenominator_finebins[0]->GetYaxis()->SetTitle("1/#sigma d#sigma/d(#Delta|y_{t}|)");
    hdiff_finebins[0]->GetXaxis()->SetTitle("#Delta|y_{t}|");
  }

  hdenominator_finebins[0]->Draw("hist");
  for (int PDFset = 1; PDFset < 41; ++PDFset) {
    hdenominator_finebins[PDFset]->Draw("hist same");
  }

  if(!drawnorm){
/*
    THStack *hs = new THStack("hs_statband", "Stat band");
    hacceptance_statdown->SetLineColor(kWhite);
    hacceptance_statdown->SetFillColor(kWhite);
    hacceptance_statdown->SetFillStyle(0);
    hs->Add(hacceptance_statdown);
    tdrStyle->SetHatchesSpacing(0.6);
    hacceptance_statup->SetFillStyle(3335);
    hacceptance_statup->SetLineColor(kWhite);
    hacceptance_statup->SetFillColor(15);
    hs->Add(hacceptance_statup);

    hs->Draw();

    hs->SetMaximum(1.25*hs->GetMaximum()/1.05);
    if(hs->GetMinimum() <0.15 *hs->GetMaximum() ) hs->SetMinimum(0.);  
    else hs->SetMinimum(0.75*hs->GetMinimum() );  

    hs->GetXaxis()->SetTitle(hacceptance->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetDecimals(kTRUE);
    hs->GetYaxis()->SetTitle("Acceptance #times Efficiency");
    hs->GetYaxis()->SetDecimals(kTRUE);
    hs->GetXaxis()->SetTitleSize(0.06);
    hs->GetXaxis()->SetLabelSize(0.05);
    hs->GetYaxis()->SetTitleSize(0.06);
    hs->GetYaxis()->SetLabelSize(0.05);
    if(!histname.Contains("lepAzimAsym2") ) hs->GetXaxis()->SetNdivisions(504,0);
    else hs->GetXaxis()->SetNdivisions(506);
    hs->GetYaxis()->SetNdivisions(508);
    hs->GetYaxis()->SetTitleOffset(1.4);
    hs->GetXaxis()->SetTitleOffset(1.0);
*/
    //hacceptance->SetMarkerSize(1.5);
    //hacceptance->Draw("hist TEXT30E");
    //hacceptance->Draw("hist E");


  }

  TLegend *leg = new TLegend(0.48,0.78,0.90,0.92);
  leg->SetFillStyle(0);
  leg->SetLineColor(0);
  leg->SetTextSize(0.036);
  leg->SetFillStyle(0);
  leg->AddEntry(hdenominator, "MC@NLO parton level","l"); 
  leg->Draw("same");

  TPaveText *pt1 = new TPaveText(0.155, 0.94, 0.41, 0.98, "brNDC");
  pt1->SetName("pt1name");
  pt1->SetBorderSize(0);
  pt1->SetFillStyle(0);

  TText *blah;
  blah = pt1->AddText("CMS Simulation, #sqrt{s} = 8 TeV");
  blah->SetTextSize(0.036);
  blah->SetTextAlign(11);   

  pt1->Draw();


     p2->cd();

     line = new TLine(hdenominator_finebins[0]->GetXaxis()->GetBinLowEdge(1), 1.0,  hdenominator_finebins[0]->GetXaxis()->GetBinUpEdge(hdenominator_finebins[0]->GetNbinsX()), 1.0);

            hdiff_finebins[0]->SetMinimum(0.921);
            hdiff_finebins[0]->SetMaximum(1.079);



            hdiff_finebins[0]->GetXaxis()->SetTitle( hdenominator_finebins[0]->GetXaxis()->GetTitle() );
            hdiff_finebins[0]->GetXaxis()->SetTitleSize( hdenominator_finebins[0]->GetXaxis()->GetTitleSize()*0.69/0.31);
            hdiff_finebins[0]->GetXaxis()->SetLabelSize( hdenominator_finebins[0]->GetXaxis()->GetLabelSize()*0.69/0.31);
            //hdiff_finebins[0]->GetXaxis()->SetLabelOffset(-0.88);

            hdiff_finebins[0]->GetYaxis()->SetTitle("#DeltaPDF_{i}/default");
            hdiff_finebins[0]->GetYaxis()->SetNdivisions(805);
            //hdiff_finebins[0]->GetYaxis()->SetTitleFont(hData_unfolded->GetYaxis()->GetTitleFont());
            hdiff_finebins[0]->GetYaxis()->SetTitleOffset(0.7);
            hdiff_finebins[0]->GetYaxis()->SetTitleSize(0.120);
            hdiff_finebins[0]->GetYaxis()->SetLabelSize( hdenominator_finebins[0]->GetYaxis()->GetLabelSize()*0.69/0.31);
            //hdiff_finebins[0]->GetYaxis()->SetLabelFont(hData_unfolded->GetYaxis()->GetLabelFont());



  hdiff_finebins[0]->Draw("hist");
  for (int PDFset = 1; PDFset < 41; ++PDFset) {
    hdiff_finebins[PDFset]->Draw("hist same");
  }

            //line->Draw();
            c1->Modified();
            c1->Update();


            p1->cd();   


  c1->Print(Form("%s_%s.pdf", accepthistname.Data(), suffixnumerator[ic]));

  f_2->Close();


  }


}
