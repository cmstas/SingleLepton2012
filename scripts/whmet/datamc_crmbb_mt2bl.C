{

  gROOT->ProcessLine(".L ~/cms2/SingleLepton2012/scripts/plotUtilities.C+");

  float ylow = 0.;
  float yhigh = 2.;

  TCanvas* c = new TCanvas();
  c->SetGrid();

  int xbins = 6;
  TH2F* haxis = new TH2F("haxis",";;Data/MC", xbins, 0, xbins, 40, ylow, yhigh);
  haxis->GetXaxis()->SetBinLabel(1,"> 50");
  haxis->GetXaxis()->SetBinLabel(2,"> 75");
  haxis->GetXaxis()->SetBinLabel(3,"> 100");
  haxis->GetXaxis()->SetBinLabel(4,"> 125");
  haxis->GetXaxis()->SetBinLabel(5,"> 150");
  haxis->GetXaxis()->SetBinLabel(6,"> 175");
  haxis->GetXaxis()->SetTitle("MET Threshold [GeV]");

  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_cr14_mtfirst/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_crs_mtfirst/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_wbbmtrw30/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_cr14nomt/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_topptrw2_tightpumva/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_sig_nvtx_split/";
  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_crs_wbbnloxsec/";

  TFile* fd = new TFile(inbase+"data_histos.root");
  TFile* f1 = new TFile(inbase+"bg_1lb_histos.root");
  //  TFile* fo = new TFile(inbase+"bg_others_histos.root");
  TFile* fo = new TFile(inbase+"bg_others_dd_histos.root");

  TString histname = "cr14_mtfirst/h_met";
  //  TString histname = "cr1_metlast_mtfirst/h_met";
  TH1F* hd_mt = (TH1F*)fd->Get(histname);
  TH1F* h1_mt = (TH1F*)f1->Get(histname);
  TH1F* ho_mt = (TH1F*)fo->Get(histname);

  TH1F* hsub_mt = hd_mt->Clone("hsub_mt");
  hsub_mt->Add(ho_mt,-1.);

  std::cout << "mt > 100:" << std::endl;
  // const int ncuts = 4;
  // const float vals[ncuts] = {50.,100.,150.,175.};
  // TGraphErrors* g = new TGraphErrors(ncuts);
  TGraphErrors* g_mt = makeMETGraph(hsub_mt,h1_mt,-0.05);
  g_mt->SetName("g_mt");
  g_mt->SetLineColor(kRed);
  g_mt->SetMarkerColor(kRed);

  TString histname2 = "cr14_mt_nm1/h_met";
  //  TString histname2 = "cr1_metlast_mt_nm1/h_met";
  TH1F* hd_mt2bl = (TH1F*)fd->Get(histname2);
  TH1F* h1_mt2bl = (TH1F*)f1->Get(histname2);
  TH1F* ho_mt2bl = (TH1F*)fo->Get(histname2);

  TH1F* hsub_mt2bl = hd_mt2bl->Clone("hsub_mt2bl");
  hsub_mt2bl->Add(ho_mt2bl,-1.);

  std::cout << "mt2bl > 200:" << std::endl;
  TGraphErrors* g_mt2bl = makeMETGraph(hsub_mt2bl,h1_mt2bl,-0.05,true);
  g_mt2bl->SetName("g_mt2bl");
  g_mt2bl->SetLineColor(kBlue);
  g_mt2bl->SetMarkerColor(kBlue);

  // for (int i = 0; i < ncuts; ++i) {
  //   Double_t err_data = 0;
  //   Double_t err_mc_1lb = 0;
  //   Double_t err_mc_others = 0;
  //   float n_data = hd->IntegralAndError(vals[i],-1,err_data);
  //   float n_mc_1lb = h1->IntegralAndError(vals[i],-1,err_mc_1lb);
  //   float n_mc_others = ho->IntegralAndError(vals[i],-1,err_mc_others);
  //   float n_sub = n_data - n_mc_others;
  //   float err_sub = sqrt(err_data**2 + err_mc_others**2);
  //   float ratio = n_sub/n_mc_1lb;
  //   float err = err_mult(n_sub,n_mc_1lb,err_sub,err_mc_1lb,ratio);
  //   g->SetPoint(i, i+0.5, ratio);
  //   g->SetPointError(i, 0., err);
  //   const bool print = true;
  //   if (print) {
  //     std::cout << "MET " << vals[i] << ": MC, data, ratio: " << Form("%.1f",n_mc_1lb) << " $\\pm$ " << Form("%.1f",err_mc_1lb)
  // 		<< " & " << n_sub << " & " << Form("%.2f",ratio) << " $\\pm$ " << Form("%.2f",err) << std::endl;
  //   }
  // }

  TLegend* leg = init_legend(0.37,0.72,0.75,0.92);
  // leg->AddEntry(g_mt,"mt > 100","lp");
  leg->AddEntry(g_mt2bl,"M_{T2}^{bl} > 200 GeV","lp");

  //  float syst = 0.25;
  float syst = 0.0;
  TLine* line_high = new TLine(0.,1.+syst,xbins,1.+syst);
  line_high->SetLineColor(kMagenta);
  line_high->SetLineWidth(3);
  line_high->SetLineStyle(2);
  TLine* line_low = new TLine(0.,1.-syst,xbins,1.-syst);
  line_low->SetLineColor(kMagenta);
  line_low->SetLineWidth(3);
  line_low->SetLineStyle(2);

  haxis->Draw();
  //  line_high->Draw("same");
  line_low->Draw("same");
  //  g_mt->Draw("ep same");
  g_mt2bl->Draw("ep same");
  leg->Draw("same");

}
