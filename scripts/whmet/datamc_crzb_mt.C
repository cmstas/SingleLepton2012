{

  gROOT->ProcessLine(".L ~/cms2/SingleLepton2012/scripts/plotUtilities.C+");

  float ylow = 0.;
  float yhigh = 2.;

  TCanvas* c = new TCanvas();
  c->SetGrid();

  // TH2F* haxis = new TH2F("haxis",";;Data/MC", 4, 0, 4, 40, ylow, yhigh);
  // haxis->GetXaxis()->SetBinLabel(1,"MET > 50");
  // haxis->GetXaxis()->SetBinLabel(2,"MET > 100");
  // haxis->GetXaxis()->SetBinLabel(3,"MET > 150");
  // haxis->GetXaxis()->SetBinLabel(4,"MET > 175");

  int xbins = 6;
  TH2F* haxis = new TH2F("haxis",";;Data/MC", xbins, 0, xbins, 40, ylow, yhigh);
  haxis->GetXaxis()->SetBinLabel(1,"> 50");
  haxis->GetXaxis()->SetBinLabel(2,"> 75");
  haxis->GetXaxis()->SetBinLabel(3,"> 100");
  haxis->GetXaxis()->SetBinLabel(4,"> 125");
  haxis->GetXaxis()->SetBinLabel(5,"> 150");
  haxis->GetXaxis()->SetBinLabel(6,"> 175");
  haxis->GetXaxis()->SetTitle("MET Threshold [GeV]");

  //  std::string inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V21_genmtvars/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_isotau_nozveto_cr13/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_cr5_mtfirst/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_wbbmtrw30/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_masslastnomt/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_topptrw2_tightpumva/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_crs/";
  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_sig_nvtx_split/";

  TFile* fd = new TFile(inbase+"data_histos.root");
  TFile* fmc = new TFile(inbase+"allbg_histos.root");

  // TH1F* hd_mt = (TH1F*)fd->Get("cr23_met_nm1/h_met");
  // TH1F* hmc_mt = (TH1F*)fmc->Get("cr23_met_nm1/h_met");
  // TH1F* hd_mt = (TH1F*)fd->Get("cr5_metlast_mtfirst/h_met");
  // TH1F* hmc_mt = (TH1F*)fmc->Get("cr5_metlast_mtfirst/h_met");
  TH1F* hd_mt = (TH1F*)fd->Get("cr5_metlast_met_nm1/h_met");
  TH1F* hmc_mt = (TH1F*)fmc->Get("cr5_metlast_met_nm1/h_met");

  std::cout << "MT > 100:" << std::endl;
  TGraphErrors* g_mt = makeMETGraph(hd_mt,hmc_mt,-0.05);
  g_mt->SetName("g_mt");
  g_mt->SetLineColor(kRed);
  g_mt->SetMarkerColor(kRed);

  TH1F* hd_mt2bl = (TH1F*)fd->Get("cr23_mt2blfirst/h_met");
  TH1F* hmc_mt2bl = (TH1F*)fmc->Get("cr23_mt2blfirst/h_met");

  TLegend* leg = init_legend(0.37,0.72,0.75,0.92);
  leg->AddEntry(g_mt,"M_{T} > 100 GeV","lp");

  //  float syst = 0.4;
  float syst = 0.1;
  TLine* line_high = new TLine(0.,1.+syst,xbins,1.+syst);
  line_high->SetLineColor(kMagenta);
  line_high->SetLineWidth(3);
  line_high->SetLineStyle(2);
  TLine* line_low = new TLine(0.,1.-syst,xbins,1.-syst);
  line_low->SetLineColor(kMagenta);
  line_low->SetLineWidth(3);
  line_low->SetLineStyle(2);

  haxis->Draw();
  line_high->Draw("same");
  //  line_low->Draw("same");
  g_mt->Draw("ep same");
  leg->Draw("same");



}
