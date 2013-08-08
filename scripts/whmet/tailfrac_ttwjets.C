{

  gROOT->ProcessLine(".L ~/cms2/SingleLepton2012/scripts/plotUtilities.C+");

  float ylow = 0.;
  float yhigh = 4.;

  TCanvas* c = new TCanvas();
  c->SetGrid();

  // TH2F* haxis = new TH2F("haxis",";;Data/MC", 4, 0, 4, 40, ylow, yhigh);
  // haxis->GetXaxis()->SetBinLabel(1,"MET > 50");
  // haxis->GetXaxis()->SetBinLabel(2,"MET > 100");
  // haxis->GetXaxis()->SetBinLabel(3,"MET > 150");
  // haxis->GetXaxis()->SetBinLabel(4,"MET > 175");

  // int xbins = 8;
  // TH2F* haxis = new TH2F("haxis",";;Ratio", xbins, 0, xbins, 40, ylow, yhigh);
  // haxis->GetXaxis()->SetBinLabel(1,"> 50");
  // haxis->GetXaxis()->SetBinLabel(2,"> 75");
  // haxis->GetXaxis()->SetBinLabel(3,"> 100");
  // haxis->GetXaxis()->SetBinLabel(4,"> 125");
  // haxis->GetXaxis()->SetBinLabel(5,"> 137");
  // haxis->GetXaxis()->SetBinLabel(6,"> 150");
  // haxis->GetXaxis()->SetBinLabel(7,"> 162");
  // haxis->GetXaxis()->SetBinLabel(8,"> 175");
  // haxis->GetXaxis()->SetTitle("MET Threshold [GeV]");

  int xbins = 6;
  TH2F* haxis = new TH2F("haxis",";;Ratio", xbins, 0, xbins, 40, ylow, yhigh);
  haxis->GetXaxis()->SetBinLabel(1,"50-75");
  haxis->GetXaxis()->SetBinLabel(2,"75-100");
  haxis->GetXaxis()->SetBinLabel(3,"100-125");
  haxis->GetXaxis()->SetBinLabel(4,"125-150");
  haxis->GetXaxis()->SetBinLabel(5,"150-175");
  haxis->GetXaxis()->SetBinLabel(6,"> 175");
  haxis->GetXaxis()->SetTitle("MET Bin [GeV]");

  //  std::string inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V21_genmtvars/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_isotau_nozveto_cr13/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_cr5_mtfirst/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_wbbmtrw30/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_masslastfinermet/";
  //  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_masslastexcmet/";
  TString inbase = "/media/data/olivito/cms2/SingleLepton2012/macros/WHLooper/output/V24_tailfrachists/";

  TFile* ftt = new TFile(inbase+"ttbar_mg_1l_histos.root");
  TFile* fwj = new TFile(inbase+"wjets_nobb_histos.root");

  TGraphErrors* g_mt = new TGraphErrors(xbins);

  std::vector<TString> dirsvec;
  //  dirsvec.push_back("mt_nm1");
  dirsvec.push_back("nomt_met50");
  dirsvec.push_back("nomt_met75");
  dirsvec.push_back("nomt_met100");
  dirsvec.push_back("nomt_met125");
  //  dirsvec.push_back("nomt_met137");
  dirsvec.push_back("nomt_met150");
  //  dirsvec.push_back("nomt_met162");
  dirsvec.push_back("nomt_met175");

  // then add prefixes to dir names
  // use getMTTailFrac etc
  for (unsigned int i=0; i < dirsvec.size(); ++i) {
    TString ttdir = "sig_bbmasslast_"+dirsvec[i];
    TString wjdir = "cr5_bbmasslast_"+dirsvec[i];
    Double_t err_tt,err_wj;
    float f_tt = getMTTailFrac(ftt,ttdir,err_tt,false);
    float f_wj = getMTTailFrac(fwj,wjdir,err_wj,false);
    float ratio = f_wj/f_tt;
    float ratio_err = err_mult(f_wj,f_tt,err_wj,err_tt,ratio);
    g_mt->SetPoint(i,i +0.5,ratio);
    g_mt->SetPointError(i,0.,ratio_err);
    std::cout << Form("%.3f",f_wj) << " $\\pm$ " << Form("%.3f",err_wj) << "  &  "
	      << Form("%.3f",f_tt) << " $\\pm$ " << Form("%.3f",err_tt) << "  &  "
	      << Form("%.1f",ratio) << " $\\pm$ " << Form("%.1f",ratio_err) << std::endl;
  }

  g_mt->SetName("g_mt");
  g_mt->SetLineColor(kRed);
  g_mt->SetMarkerColor(kRed);

  TLegend* leg = init_legend(0.37,0.72,0.75,0.92);
  leg->AddEntry(g_mt,"M_{T} > 100 GeV","lp");
  //  leg->AddEntry(g_mt2bl,"M_{T2}^{bl} > 200 GeV","lp");

  //  float syst = 0.4;
  float syst = 0.5;
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
  //  g_mt2bl->Draw("ep same");
  leg->Draw("same");



}
