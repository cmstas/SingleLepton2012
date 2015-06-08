#include "TH1D.h"
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

TH1D* histo1;
TH1D* histo2;
TH1D* histo3;
TH1D* datahisto;

TH2D* histo2D;
TH1D* histo2D_ProjectionX[13];



void resolutionplots(TString histname = "ttMasspull",  int drawlog =0, double rangelow = -1, double rangehigh = 3, double rangeylow = 0, double rangeyhigh = 0, int rebin = 1,TString FName1 = "SIGoutput/ttdl_mcatnlo_histos.root"){
  setTDRStyle();

  TString xaxistitle;


  if(histname == "lep_charge_asymmetry") xaxistitle = "#Delta(#Delta|#eta_{l}|)";
  if(histname == "lep_azimuthal_asymmetry") xaxistitle = "#Delta(#Delta#phi_{ll})";
  if(histname == "lep_azimuthal_asymmetry2") xaxistitle = "#Delta(#Delta#phi_{ll})";
  if(histname == "top_rapiditydiff_cms") xaxistitle = "#Delta rapiditydiff";
  if(histname == "top_pseudorapiditydiff_cms") xaxistitle = "#Delta pseudorapiditydiff";
  if(histname == "top_rapiditydiff_Marco") xaxistitle = "#Delta(#Delta|y_{t}|)";
  if(histname == "top_costheta_cms") xaxistitle = "#Deltacos(#theta_{t})";
  if(histname == "lepPlus_costheta_cms") xaxistitle = "#Deltacos(#theta*)";
  if(histname == "lepMinus_costheta_cms") xaxistitle = "#Deltacos(#theta*)";
  if(histname == "lep_costheta_cms") xaxistitle = "#Deltacos(#theta*)";
  if(histname == "top_spin_correlation") xaxistitle = "#Deltac1c2";
  if(histname == "lep_cos_opening_angle") xaxistitle = "#Deltacos(#phi)";
  if(histname == "tt_mass") xaxistitle = "#DeltaM_{t#bart}";
  if(histname == "ttRapidity2") xaxistitle = "#Deltay_{t#bart}";
  if(histname == "absttRapidity2") xaxistitle = "#Delta|y|_{t#bart}";
  if(histname == "tt_pT") xaxistitle = "#Deltap_{T}^{t#bart}";
  if(histname == "top1_pt") xaxistitle = "#Delta p_{T}^{t}";
  if(histname == "top2_pt") xaxistitle = "#Delta p_{T}^{#bart}";



  bool is2D = rangeyhigh == 0 ? 0 : 1; 
  std::cout << "Opening " << FName1.Data() << "\n";
  TFile *f_1         = TFile::Open(FName1.Data());  
  histo1 = (TH1D*)f_1->Get(Form("h_sig_%s_Delta_all", histname.Data()));
  std::cout << "hist " << histo1->GetName() << " with entries " << histo1->GetEntries() << std::endl;

  histo2D = (TH2D*)f_1->Get(Form("h_sig_%s_Delta_2D_all", histname.Data()));
  std::cout << "hist " << histo2D->GetName() << " with entries " << histo2D->GetEntries() << std::endl;

  gStyle->SetOptStat(1001001100);
  if(drawlog && !is2D) gStyle->SetOptLogy(1);
  else gStyle->SetOptLogy(0);

  if(drawlog && is2D) gStyle->SetOptLogz(1);
  else gStyle->SetOptLogz(0);

  if(histname.Contains("lep_charge_asymmetry") || histname.Contains("lep_azimuthal_asymmetry")) gStyle->SetOptLogy(1);


  int nbins = histo2D->GetNbinsY();
  int ndiv = 0;
  
  if(nbins==3) ndiv = 2;
  if(nbins==6) ndiv = 3;
  if(nbins==12) ndiv = 4;
  cout<<nbins<<" "<<ndiv<<endl;

  TCanvas *c1 = new TCanvas();
  c1->Divide(ndiv,ndiv);

  for (int i = 0; i < nbins+1; ++i)
  {
    if(i==0) histo2D_ProjectionX[i] = histo2D->ProjectionX();
    else histo2D_ProjectionX[i] = histo2D->ProjectionX(Form("px%i", i), i, i, "");

    rangelow = histo2D_ProjectionX[i]->GetXaxis()->GetXmin(); 
    rangehigh = histo2D_ProjectionX[i]->GetXaxis()->GetXmax();

    //if(!is2D) histo2D_ProjectionX[i]->Rebin(rebin);
    histo2D_ProjectionX[i]->GetXaxis()->SetRangeUser(rangelow,rangehigh);
    if(is2D) histo2D_ProjectionX[i]->GetYaxis()->SetRangeUser(rangeylow,rangeyhigh);
    histo2D_ProjectionX[i]->SetLineColor(kBlue);
    histo2D_ProjectionX[i]-> SetFillColor(0);
    if( !( histname.Contains("lep_charge_asymmetry") || histname.Contains("lep_azimuthal_asymmetry") ) ) histo2D_ProjectionX[i]-> SetMinimum(0.);

    if(i==0) histo2D_ProjectionX[i]->GetXaxis()->SetTitle(xaxistitle+Form(" (inclusive)",i));
    else histo2D_ProjectionX[i]->GetXaxis()->SetTitle(xaxistitle+Form(" (bin %i)",i));
    histo2D_ProjectionX[i]->GetXaxis()->SetTitleOffset(1.1);
    histo2D_ProjectionX[i]->GetXaxis()->SetTitleSize(0.05);
    histo2D_ProjectionX[i]->GetXaxis()->SetLabelSize(0.05);
    histo2D_ProjectionX[i]->GetYaxis()->SetTitle("Events/bin");
    histo2D_ProjectionX[i]->GetYaxis()->SetTitleOffset(1.6);
    histo2D_ProjectionX[i]->GetYaxis()->SetTitleSize(0.05);
    histo2D_ProjectionX[i]->GetYaxis()->SetLabelSize(0.05);

    if(nbins==3 || i==0) c1->cd(i+1);
    else c1->cd(i+ndiv);

    if(!is2D) histo2D_ProjectionX[i]->Draw("hist");
    else histo2D_ProjectionX[i]->Draw("COLZ");

    //if(!is2D) histo1->Fit("gaus","","",(rangelow+rangehigh)/2. + (rangelow-rangehigh)/16., (rangelow+rangehigh)/2. - (rangelow-rangehigh)/16.);

    /*
    TLegend *leg = new TLegend(0.74,0.76,0.92,0.92);
    leg->SetFillStyle(0);
    leg->SetLineColor(0);
    leg->SetTextSize(0.035);
    leg->SetFillStyle(0);
    leg->AddEntry(histo1, "ttdil","l");

    leg->Draw("same");
    */

  /*
    TPaveText *pt1 = new TPaveText(0.20, 0.76, 0.40, 0.91, "brNDC");
    pt1->SetName("pt1name");
    pt1->SetBorderSize(0);
    pt1->SetFillStyle(0);
    
    TText *blah;
    blah = pt1->AddText("CMS Preliminary, 5.0 fb^{-1} at  #sqrt{s}=7 TeV");
    blah->SetTextSize(0.032);
    blah->SetTextAlign(11);  

    pt1->Draw();
  */
  }

  if(!drawlog) c1->Print(Form("resolution_%s.pdf", histname.Data()));
  else c1->Print(Form("resolution_%s_log.pdf", histname.Data()));

  f_1->Close();
  
}
