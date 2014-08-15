{

  gROOT->ProcessLine(".L ~/cms2/SingleLepton2012/scripts/whmet/plotUtilities.C+");

  float ylow = 0.;
  float yhigh = 2.5;

  TCanvas* c = new TCanvas();
  c->SetGrid();

  int xbins = 6;
  TH2F* haxis = new TH2F("haxis",";;Data/SM", xbins, 0, xbins, 40, ylow, yhigh);
  haxis->GetXaxis()->SetBinLabel(1,"> 50");
  haxis->GetXaxis()->SetBinLabel(2,"> 75");
  haxis->GetXaxis()->SetBinLabel(3,"> 100");
  haxis->GetXaxis()->SetBinLabel(4,"> 125");
  haxis->GetXaxis()->SetBinLabel(5,"> 150");
  haxis->GetXaxis()->SetBinLabel(6,"> 175");
  haxis->GetXaxis()->SetTitle("E_{T}^{miss} Threshold [GeV]");

  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_crs/";
  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_crs_blinded/";

  TFile* fd = new TFile(inbase+"data_histos.root");
  TFile* fmc = new TFile(inbase+"allbg_histos.root");

  TH1F* hd_mt = (TH1F*)fd->Get("cr23_mt2bl_nm1/h_met");
  TH1F* hmc_mt = (TH1F*)fmc->Get("cr23_mt2bl_nm1/h_met");

  std::cout << "MT > 100:" << std::endl;
  TGraphErrors* g_mt = makeMETGraph(hd_mt,hmc_mt,-0.05);
  g_mt->SetName("g_mt");
  g_mt->SetLineColor(kRed);
  g_mt->SetMarkerColor(kRed);

  TH1F* hd_mt2bl = (TH1F*)fd->Get("cr23_mt2blfirst/h_met");
  TH1F* hmc_mt2bl = (TH1F*)fmc->Get("cr23_mt2blfirst/h_met");

  std::cout << "MT2bl > 200:" << std::endl;
  TGraphErrors* g_mt2bl = makeMETGraph(hd_mt2bl,hmc_mt2bl,0.05);
  g_mt2bl->SetName("g_mt2bl");
  g_mt2bl->SetLineColor(kBlue);
  g_mt2bl->SetMarkerColor(kBlue);
  g_mt2bl->SetMarkerStyle(21);

  TGraph* dummy = new TGraph();
  dummy->SetLineColor(kMagenta);
  dummy->SetLineWidth(3);
  dummy->SetLineStyle(1);

  TLegend* leg = init_legend(0.55,0.72,0.96,0.92);
  leg->AddEntry(g_mt,"M_{T} > 100 GeV","lp");
  leg->AddEntry(g_mt2bl,"M_{T2}^{bl} > 200 GeV","lp");
  leg->AddEntry(dummy,"Scale factor #pm uncert.","l");

  float sf = 1.0;
  float syst = 0.4;
  float xmin = 2.;

  TLine* line_cent = new TLine(xmin,sf,xbins,sf);
  line_cent->SetLineColor(kMagenta);
  line_cent->SetLineWidth(3);
  line_cent->SetLineStyle(1);
  TLine* line_high = new TLine(xmin,sf+syst,xbins,sf+syst);
  line_high->SetLineColor(kMagenta);
  line_high->SetLineWidth(3);
  line_high->SetLineStyle(2);
  TLine* line_low = new TLine(xmin,sf-syst,xbins,sf-syst);
  line_low->SetLineColor(kMagenta);
  line_low->SetLineWidth(3);
  line_low->SetLineStyle(2);

  haxis->Draw();
  haxis->GetXaxis()->SetRangeUser(xmin,xbins);
  line_cent->Draw("same");
  line_high->Draw("same");
  line_low->Draw("same");
  g_mt->Draw("ep same");
  g_mt2bl->Draw("ep same");
  leg->Draw("same");

  float sf_center = 1.903;
  float y_offset = 0.04;
  float sf_xmin = 3.97;
  float sf_xmax = 4.33;
  TLine *line_exp = new TLine();
  line_exp->SetLineColor(kMagenta);
  line_exp->SetLineWidth(3);
  line_exp->SetLineStyle(2);
  line_exp->DrawLine(sf_xmin,sf_center+y_offset,sf_xmax,sf_center+y_offset);
  line_exp->DrawLine(sf_xmin,sf_center-y_offset,sf_xmax,sf_center-y_offset);

  TLatex *text = new TLatex();
  text->SetNDC();

  text->SetTextFont(42);
  text->SetTextSize(0.038);
  text->SetTextAlign(20);
  //  TString label = "CMS          #sqrt{s} = 8 TeV         L = 19.5 fb^{-1}";
  TString label = "CMS Unpublished    #sqrt{s} = 8 TeV    L = 19.5 fb^{-1}";
  text->DrawLatex(0.56,0.95,label);

  // TLatex *text = new TLatex();
  // text->SetNDC();
  // text->SetTextSize(0.03);
  // //  text->DrawLatex(0.2,0.88,"CMS Preliminary");
  // text->DrawLatex(0.2,0.88,"CMS Unpublished");
  // //text->DrawLatex(0.2,0.83,"0.98 fb^{-1} at #sqrt{s} = 7 TeV");
  // text->DrawLatex(0.2,0.83,"#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = 19.5 fb^{-1}");
  // text->DrawLatex(0.2,0.78,"CR-2l");

  TString savedir = "/media/data/olivito/cms2/SingleLepton2012/plots/May27_14/";
  TString plotname = "datamc_crdl_mt2bl_mt";
  c->SaveAs(savedir + "/" + plotname + ".eps");
  c->SaveAs(savedir + "/" + plotname + ".pdf");
  c->SaveAs(savedir + "/" + plotname + ".png");

}
