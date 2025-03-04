#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <vector>
#include <string>
#include "TMath.h"
#include "TDirectory.h"
#include "topAFB_looper.h"
void topAFB_looper::bookHistos(const char *prefix, int nchannel, int nhists) {
  //  Book histograms...
  //  Naming Convention:
  //  Prefix comes from the sample and it is passed to the scanning function
  //  Suffix is "ee" "em" "em" "all" which depends on the final state
  //  For example: histogram named tt_hnJet_ee would be the Njet distribution
  //  for the ee final state in the ttbar sample.
  
  // MAKE SURE TO CAL SUMW2 FOR EACH 1D HISTOGRAM BEFORE FILLING!!!!!!
  
  cout << "Begin book histos..." << endl;

  TDirectory *rootdir = gDirectory->GetDirectory("Rint:");
  rootdir->cd();
    
  char suffixall[4][4]  = {"ee", "mm", "em", "all"};
  //char njetCh[4][5]     = {"0j", "1j", "2j", "allj"};

  //Define binning. These numbers are copied from StopTreeLooper.cc.
  int nbins = 240; //240 bins in the asymmetry variable makes more sense than 80 since it can be divided into 12 and we will be rebinning to 6 or 12 bins for the acceptance histograms

  int nbinsmtt = 120;
  int nbinsttpt = 300;
  int nbinsttrapidity2 = 300;

  double minmtt = 0.0;
  double minttpt = 0.0;
  double minttrapidity2 = 0.0;

  double maxmtt = 1200.0;
  double maxttpt = 300.0;
  double maxttrapidity2 = 1.5;


  double pi = TMath::Pi();

  
  //the plots not differ by number of jets
  for (int i=0; i<4; i++) {

    for (int j = 0; j < 41; j++) {
      char suffix[7];
      sprintf(suffix, "PDFset%i_%s", j, suffixall[i]);
/*
      hnJet[i][j] = new TH1D(Form("%s_hnJet_%s",prefix,suffix),Form("%s_nJet_%s",prefix,suffix),9,0,9);	
      hnJet[i][j]->GetXaxis()->SetTitle("Number of jets");
      hnJet[i][j]->Sumw2();

      hnBtagJet[i][j] = new TH1D(Form("%s_hnBtagJet_%s",prefix,suffix),Form("%s_nBtagJet_%s",prefix,suffix),6,0,6);	
      hnBtagJet[i][j]->GetXaxis()->SetTitle("Number of b tagged jets");
      hnBtagJet[i][j]->Sumw2();

      hnVtx[i][j] = new TH1D(Form("%s_hnVtx_%s",prefix,suffix),Form("%s_nVtx_%s",prefix,suffix),21,0.,21.);	
      hnVtx[i][j]->GetXaxis()->SetTitle("Number of vertices");
      hnVtx[i][j]->Sumw2();


      hAMWTweightnojetsmear[i][j] = new TH1D(Form("%s_hAMWTweightnojetsmear_%s",prefix,suffix),Form("%s_AMWTweightnojetsmear_%s",prefix,suffix),120,-0.2,1.0);
      hAMWTweightnojetsmear[i][j]->GetXaxis()->SetTitle("AMWT weight (no jet smearing)");
      hAMWTweightnojetsmear[i][j]->Sumw2();


      
      hlepChargeAsym_2d[i][j] = new TH2D(Form("%s_hlepChargeAsym2d_%s",prefix,suffix),Form("%s_lepChargeAsym2d_%s",prefix,suffix),80,-4,4, 80,-4,4);
      hlepChargeAsym_2d[i][j]->GetXaxis()->SetTitle("Charge_Asymmetry_lep_gen");
      hlepChargeAsym_2d[i][j]->GetYaxis()->SetTitle("Charge_Asymmetry_lep_rec");
      hlepChargeAsym_2d[i][j]->Sumw2();

      hlepAzimAsym_2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2d_%s",prefix,suffix),Form("%s_lepAzimAsym2d_%s",prefix,suffix),80,-1,1, 80,-1,1);
      hlepAzimAsym_2d[i][j]->GetXaxis()->SetTitle("Azimuthal_Asymmetry_lep_gen");
      hlepAzimAsym_2d[i][j]->GetYaxis()->SetTitle("Azimuthal_Asymmetry_lep_rec");
      hlepAzimAsym_2d[i][j]->Sumw2();

      htopSpinCorr_2d[i][j] = new TH2D(Form("%s_htopSpinCorr2d_%s",prefix,suffix),Form("%s_topSpinCorr2d_%s",prefix,suffix),80,-1,1,80,-1,1);
      htopSpinCorr_2d[i][j]->GetXaxis()->SetTitle("Spin_Correlation_top_gen");
      htopSpinCorr_2d[i][j]->GetYaxis()->SetTitle("Spin_Correlation_top_rec");
      htopSpinCorr_2d[i][j]->Sumw2();

      hlepCosOpeningAngle_2d[i][j] = new TH2D(Form("%s_hlepCosOpeningAngle2d_%s",prefix,suffix),Form("%s_lepCosOpeningAngle2d_%s",prefix,suffix),80,-1,1,80,-1,1);
      hlepCosOpeningAngle_2d[i][j]->GetXaxis()->SetTitle("cos(#phi) gen");
      hlepCosOpeningAngle_2d[i][j]->GetYaxis()->SetTitle("cos(#phi) reco");
      hlepCosOpeningAngle_2d[i][j]->Sumw2();
      
      htopCosTheta_2d[i][j] = new TH2D(Form("%s_htopCosTheta2d_%s",prefix,suffix),Form("%s_topCosTheta2d_%s",prefix,suffix),80,-1,1, 80,-1,1);
      htopCosTheta_2d[i][j]->GetXaxis()->SetTitle("Cos(theta)_top_gen");
      htopCosTheta_2d[i][j]->GetYaxis()->SetTitle("Cos(theta)_top_rec");
      htopCosTheta_2d[i][j]->Sumw2();

      hlepCosTheta_2d[i][j] = new TH2D(Form("%s_hlepCosTheta2d_%s",prefix,suffix),Form("%s_lepCosTheta2d_%s",prefix,suffix),80,-1,1, 80,-1,1);
      hlepCosTheta_2d[i][j]->GetXaxis()->SetTitle("Cos(theta)_lep_gen");
      hlepCosTheta_2d[i][j]->GetYaxis()->SetTitle("Cos(theta)_lep_rec");
      hlepCosTheta_2d[i][j]->Sumw2();
      */
      
      
      hlepChargeAsym_gen[i][j] = new TH1D(Form("%s_hlepChargeAsymGen_%s",prefix,suffix),Form("%s_lepChargeAsymGen_%s",prefix,suffix),nbins,-2,2);
      hlepChargeAsym_gen[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_gen[i][j]->Sumw2();

      hlepAzimAsym_gen[i][j] = new TH1D(Form("%s_hlepAzimAsymGen_%s",prefix,suffix),Form("%s_lepAzimAsymGen_%s",prefix,suffix),nbins,-pi,pi);
      hlepAzimAsym_gen[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym_gen[i][j]->Sumw2();

      hlepAzimAsym2_gen[i][j] = new TH1D(Form("%s_hlepAzimAsym2Gen_%s",prefix,suffix),Form("%s_lepAzimAsym2Gen_%s",prefix,suffix),nbins,0,pi);
      hlepAzimAsym2_gen[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_gen[i][j]->Sumw2();

      htopSpinCorr_gen[i][j] = new TH1D(Form("%s_htopSpinCorrGen_%s",prefix,suffix),Form("%s_topSpinCorrGen_%s",prefix,suffix),nbins,-1,1);
      htopSpinCorr_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_gen[i][j]->Sumw2();

      hlepCosOpeningAngle_gen[i][j] = new TH1D(Form("%s_hlepCosOpeningAngleGen_%s",prefix,suffix),Form("%s_lepCosOpeningAngleGen_%s",prefix,suffix),nbins,-1,1);
      hlepCosOpeningAngle_gen[i][j]->GetXaxis()->SetTitle("cos(#phi)");
      hlepCosOpeningAngle_gen[i][j]->Sumw2();
      
      htopCosTheta_gen[i][j] = new TH1D(Form("%s_htopCosThetaGen_%s",prefix,suffix),Form("%s_topCosThetaGen_%s",prefix,suffix),nbins,-1,1);
      htopCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_gen[i][j]->Sumw2();

      hlepCosTheta_gen[i][j] = new TH1D(Form("%s_hlepCosThetaGen_%s",prefix,suffix),Form("%s_lepCosThetaGen_%s",prefix,suffix),nbins,-1,1);
      hlepCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_gen[i][j]->Sumw2();

      hlepPlusCosTheta_gen[i][j] = new TH1D(Form("%s_hlepPlusCosThetaGen_%s",prefix,suffix),Form("%s_lepPlusCosThetaGen_%s",prefix,suffix),nbins,-1,1);
      hlepPlusCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_gen[i][j]->Sumw2();
      
      hlepMinusCosTheta_gen[i][j] = new TH1D(Form("%s_hlepMinusCosThetaGen_%s",prefix,suffix),Form("%s_lepMinusCosThetaGen_%s",prefix,suffix),nbins,-1,1);
      hlepMinusCosTheta_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_gen[i][j]->Sumw2();
            
      hpseudorapiditydiff_gen[i][j] = new TH1D(Form("%s_hpseudorapiditydiffGen_%s",prefix,suffix),Form("%s_pseudorapiditydiffGen_%s",prefix,suffix),nbins,-2,2);
      hpseudorapiditydiff_gen[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_gen[i][j]->Sumw2();
      
      hrapiditydiff_gen[i][j] = new TH1D(Form("%s_hrapiditydiffGen_%s",prefix,suffix),Form("%s_rapiditydiffGen_%s",prefix,suffix),nbins,-2,2);
      hrapiditydiff_gen[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_gen[i][j]->Sumw2();
      
      hrapiditydiffMarco_gen[i][j] = new TH1D(Form("%s_hrapiditydiffMarcoGen_%s",prefix,suffix),Form("%s_rapiditydiffMarcoGen_%s",prefix,suffix),nbins,-2,2);
      hrapiditydiffMarco_gen[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_gen[i][j]->Sumw2();      
      
      
      
      
      //2d histos for asyms vs mass at gen level
      hlepChargeAsym_gen2d[i][j] = new TH2D(Form("%s_hlepChargeAsymGen2d_%s",prefix,suffix),Form("%s_lepChargeAsymGen2d_%s",prefix,suffix),nbins,-2,2, nbinsmtt, minmtt, maxmtt);
      hlepChargeAsym_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepChargeAsym_gen2d[i][j]->Sumw2();

      hlepAzimAsym_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsymGen2d_%s",prefix,suffix),Form("%s_lepAzimAsymGen2d_%s",prefix,suffix),nbins,-pi,pi, nbinsmtt, minmtt, maxmtt);
      hlepAzimAsym_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym_gen2d[i][j]->Sumw2();

      hlepAzimAsym2_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2Gen2d_%s",prefix,suffix),Form("%s_lepAzimAsym2Gen2d_%s",prefix,suffix),nbins,0,pi, nbinsmtt, minmtt, maxmtt);
      hlepAzimAsym2_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym2_gen2d[i][j]->Sumw2();

      htopSpinCorr_gen2d[i][j] = new TH2D(Form("%s_htopSpinCorrGen2d_%s",prefix,suffix),Form("%s_topSpinCorrGen2d_%s",prefix,suffix),nbins,-1,1, nbinsmtt, minmtt, maxmtt);
      htopSpinCorr_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      htopSpinCorr_gen2d[i][j]->Sumw2();

      hlepCosOpeningAngle_gen2d[i][j] = new TH2D(Form("%s_hlepCosOpeningAngleGen2d_%s",prefix,suffix),Form("%s_lepCosOpeningAngleGen2d_%s",prefix,suffix),nbins,-1,1, nbinsmtt, minmtt, maxmtt);
      hlepCosOpeningAngle_gen2d[i][j]->GetXaxis()->SetTitle("cos(#phi)");
      hlepCosOpeningAngle_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepCosOpeningAngle_gen2d[i][j]->Sumw2();
      
      htopCosTheta_gen2d[i][j] = new TH2D(Form("%s_htopCosThetaGen2d_%s",prefix,suffix),Form("%s_topCosThetaGen2d_%s",prefix,suffix),nbins,-1,1, nbinsmtt, minmtt, maxmtt);
      htopCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      htopCosTheta_gen2d[i][j]->Sumw2();

      hlepCosTheta_gen2d[i][j] = new TH2D(Form("%s_hlepCosThetaGen2d_%s",prefix,suffix),Form("%s_lepCosThetaGen2d_%s",prefix,suffix),nbins,-1,1, nbinsmtt, minmtt, maxmtt);
      hlepCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepCosTheta_gen2d[i][j]->Sumw2();

      hlepPlusCosTheta_gen2d[i][j] = new TH2D(Form("%s_hlepPlusCosThetaGen2d_%s",prefix,suffix),Form("%s_lepPlusCosThetaGen2d_%s",prefix,suffix),nbins,-1,1, nbinsmtt, minmtt, maxmtt);
      hlepPlusCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepPlusCosTheta_gen2d[i][j]->Sumw2();
      
      hlepMinusCosTheta_gen2d[i][j] = new TH2D(Form("%s_hlepMinusCosThetaGen2d_%s",prefix,suffix),Form("%s_lepMinusCosThetaGen2d_%s",prefix,suffix),nbins,-1,1, nbinsmtt, minmtt, maxmtt);
      hlepMinusCosTheta_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hlepMinusCosTheta_gen2d[i][j]->Sumw2();
            
      hpseudorapiditydiff_gen2d[i][j] = new TH2D(Form("%s_hpseudorapiditydiffGen2d_%s",prefix,suffix),Form("%s_pseudorapiditydiffGen2d_%s",prefix,suffix),nbins,-2,2, nbinsmtt, minmtt, maxmtt);
      hpseudorapiditydiff_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hpseudorapiditydiff_gen2d[i][j]->Sumw2();
      
      hrapiditydiff_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffGen2d_%s",prefix,suffix),Form("%s_rapiditydiffGen2d_%s",prefix,suffix),nbins,-2,2, nbinsmtt, minmtt, maxmtt);
      hrapiditydiff_gen2d[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiff_gen2d[i][j]->Sumw2();
      
      hrapiditydiffMarco_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffMarcoGen2d_%s",prefix,suffix),Form("%s_rapiditydiffMarcoGen2d_%s",prefix,suffix),nbins,-2,2, nbinsmtt, minmtt, maxmtt);
      hrapiditydiffMarco_gen2d[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_gen2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiffMarco_gen2d[i][j]->Sumw2();



      hlepChargeAsym_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepChargeAsymttpTGen2d_%s",prefix,suffix),Form("%s_lepChargeAsymttpTGen2d_%s",prefix,suffix),nbins,-2,2, nbinsttpt, minttpt, maxttpt);
      hlepChargeAsym_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepChargeAsym_ttpT_gen2d[i][j]->Sumw2();

      hlepAzimAsym_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsymttpTGen2d_%s",prefix,suffix),Form("%s_lepAzimAsymttpTGen2d_%s",prefix,suffix),nbins,-pi,pi, nbinsttpt, minttpt, maxttpt);
      hlepAzimAsym_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym_ttpT_gen2d[i][j]->Sumw2();

      hlepAzimAsym2_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2ttpTGen2d_%s",prefix,suffix),Form("%s_lepAzimAsym2ttpTGen2d_%s",prefix,suffix),nbins,0,pi, nbinsttpt, minttpt, maxttpt);
      hlepAzimAsym2_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepAzimAsym2_ttpT_gen2d[i][j]->Sumw2();

      htopSpinCorr_ttpT_gen2d[i][j] = new TH2D(Form("%s_htopSpinCorrttpTGen2d_%s",prefix,suffix),Form("%s_topSpinCorrttpTGen2d_%s",prefix,suffix),nbins,-1,1, nbinsttpt, minttpt, maxttpt);
      htopSpinCorr_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      htopSpinCorr_ttpT_gen2d[i][j]->Sumw2();

      hlepCosOpeningAngle_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepCosOpeningAnglettpTGen2d_%s",prefix,suffix),Form("%s_lepCosOpeningAnglettpTGen2d_%s",prefix,suffix),nbins,-1,1, nbinsttpt, minttpt, maxttpt);
      hlepCosOpeningAngle_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#phi)");
      hlepCosOpeningAngle_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepCosOpeningAngle_ttpT_gen2d[i][j]->Sumw2();
      
      htopCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_htopCosThetattpTGen2d_%s",prefix,suffix),Form("%s_topCosThetattpTGen2d_%s",prefix,suffix),nbins,-1,1, nbinsttpt, minttpt, maxttpt);
      htopCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      htopCosTheta_ttpT_gen2d[i][j]->Sumw2();

      hlepCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepCosThetattpTGen2d_%s",prefix,suffix),Form("%s_lepCosThetattpTGen2d_%s",prefix,suffix),nbins,-1,1, nbinsttpt, minttpt, maxttpt);
      hlepCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepCosTheta_ttpT_gen2d[i][j]->Sumw2();

      hlepPlusCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepPlusCosThetattpTGen2d_%s",prefix,suffix),Form("%s_lepPlusCosThetattpTGen2d_%s",prefix,suffix),nbins,-1,1, nbinsttpt, minttpt, maxttpt);
      hlepPlusCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepPlusCosTheta_ttpT_gen2d[i][j]->Sumw2();
      
      hlepMinusCosTheta_ttpT_gen2d[i][j] = new TH2D(Form("%s_hlepMinusCosThetattpTGen2d_%s",prefix,suffix),Form("%s_lepMinusCosThetattpTGen2d_%s",prefix,suffix),nbins,-1,1, nbinsttpt, minttpt, maxttpt);
      hlepMinusCosTheta_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hlepMinusCosTheta_ttpT_gen2d[i][j]->Sumw2();
            
      hpseudorapiditydiff_ttpT_gen2d[i][j] = new TH2D(Form("%s_hpseudorapiditydiffttpTGen2d_%s",prefix,suffix),Form("%s_pseudorapiditydiffttpTGen2d_%s",prefix,suffix),nbins,-2,2, nbinsttpt, minttpt, maxttpt);
      hpseudorapiditydiff_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hpseudorapiditydiff_ttpT_gen2d[i][j]->Sumw2();
      
      hrapiditydiff_ttpT_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffttpTGen2d_%s",prefix,suffix),Form("%s_rapiditydiffttpTGen2d_%s",prefix,suffix),nbins,-2,2, nbinsttpt, minttpt, maxttpt);
      hrapiditydiff_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiff_ttpT_gen2d[i][j]->Sumw2();
      
      hrapiditydiffMarco_ttpT_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffMarcottpTGen2d_%s",prefix,suffix),Form("%s_rapiditydiffMarcottpTGen2d_%s",prefix,suffix),nbins,-2,2, nbinsttpt, minttpt, maxttpt);
      hrapiditydiffMarco_ttpT_gen2d[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_ttpT_gen2d[i][j]->GetYaxis()->SetTitle("pT_{t#bar{t}} (GeV/c^{2})");
      hrapiditydiffMarco_ttpT_gen2d[i][j]->Sumw2();


    
      hlepChargeAsym_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepChargeAsymttRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepChargeAsymttRapidity2Gen2d_%s",prefix,suffix),nbins,-2,2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hlepChargeAsym_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}|");
      hlepChargeAsym_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepChargeAsym_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepAzimAsym_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsymttRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepAzimAsymttRapidity2Gen2d_%s",prefix,suffix),nbins,-pi,pi, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hlepAzimAsym_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepAzimAsym_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepAzimAsym2_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepAzimAsym2ttRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepAzimAsym2ttRapidity2Gen2d_%s",prefix,suffix),nbins,0,pi, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hlepAzimAsym2_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepAzimAsym2_ttRapidity2_gen2d[i][j]->Sumw2();

      htopSpinCorr_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_htopSpinCorrttRapidity2Gen2d_%s",prefix,suffix),Form("%s_topSpinCorrttRapidity2Gen2d_%s",prefix,suffix),nbins,-1,1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      htopSpinCorr_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      htopSpinCorr_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepCosOpeningAngle_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepCosOpeningAnglettRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepCosOpeningAnglettRapidity2Gen2d_%s",prefix,suffix),nbins,-1,1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hlepCosOpeningAngle_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#phi)");
      hlepCosOpeningAngle_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepCosOpeningAngle_ttRapidity2_gen2d[i][j]->Sumw2();
      
      htopCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_htopCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_topCosThetattRapidity2Gen2d_%s",prefix,suffix),nbins,-1,1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      htopCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      htopCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepCosThetattRapidity2Gen2d_%s",prefix,suffix),nbins,-1,1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hlepCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();

      hlepPlusCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepPlusCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepPlusCosThetattRapidity2Gen2d_%s",prefix,suffix),nbins,-1,1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hlepPlusCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepPlusCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();
      
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hlepMinusCosThetattRapidity2Gen2d_%s",prefix,suffix),Form("%s_lepMinusCosThetattRapidity2Gen2d_%s",prefix,suffix),nbins,-1,1, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hlepMinusCosTheta_ttRapidity2_gen2d[i][j]->Sumw2();
            
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hpseudorapiditydiffttRapidity2Gen2d_%s",prefix,suffix),Form("%s_pseudorapiditydiffttRapidity2Gen2d_%s",prefix,suffix),nbins,-2,2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hpseudorapiditydiff_ttRapidity2_gen2d[i][j]->Sumw2();
      
      hrapiditydiff_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffttRapidity2Gen2d_%s",prefix,suffix),Form("%s_rapiditydiffttRapidity2Gen2d_%s",prefix,suffix),nbins,-2,2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hrapiditydiff_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hrapiditydiff_ttRapidity2_gen2d[i][j]->Sumw2();
      
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j] = new TH2D(Form("%s_hrapiditydiffMarcottRapidity2Gen2d_%s",prefix,suffix),Form("%s_rapiditydiffMarcottRapidity2Gen2d_%s",prefix,suffix),nbins,-2,2, nbinsttrapidity2, minttrapidity2, maxttrapidity2);
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j]->GetYaxis()->SetTitle("|y_{t#bar{t}}|");
      hrapiditydiffMarco_ttRapidity2_gen2d[i][j]->Sumw2();



            
      /*
      //Reco-Gen hists for asyms
      hlepChargeAsymGenDiff[i][j] = new TH1D(Form("%s_hlepChargeAsymGenDiff_%s",prefix,suffix),Form("%s_lepChargeAsymGenDiff_%s",prefix,suffix),80,-4,4);
      hlepChargeAsymGenDiff[i][j]->GetXaxis()->SetTitle("|#eta_{l^{+}}| - |#eta_{l^{-}}| (reco-gen)");
      hlepChargeAsymGenDiff[i][j]->Sumw2();

      hlepAzimAsymGenDiff[i][j] = new TH1D(Form("%s_hlepAzimAsymGenDiff_%s",prefix,suffix),Form("%s_lepAzimAsymGenDiff_%s",prefix,suffix),80,-2,2);
      hlepAzimAsymGenDiff[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{l^{+}l^{-}}) (reco-gen)");
      hlepAzimAsymGenDiff[i][j]->Sumw2();

      hlepAzimAsym2GenDiff[i][j] = new TH1D(Form("%s_hlepAzimAsym2GenDiff_%s",prefix,suffix),Form("%s_lepAzimAsym2GenDiff_%s",prefix,suffix),80,-pi,pi);
      hlepAzimAsym2GenDiff[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}} (reco-gen)");
      hlepAzimAsym2GenDiff[i][j]->Sumw2();

      htopSpinCorrGenDiff[i][j] = new TH1D(Form("%s_htopSpinCorrGenDiff_%s",prefix,suffix),Form("%s_topSpinCorrGenDiff_%s",prefix,suffix),80,-2,2);
      htopSpinCorrGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} ) (reco-gen)");
      htopSpinCorrGenDiff[i][j]->Sumw2();

      hlepCosOpeningAngleGenDiff[i][j] = new TH1D(Form("%s_hlepCosOpeningAngleGenDiff_%s",prefix,suffix),Form("%s_lepCosOpeningAngleGenDiff_%s",prefix,suffix),80,-2,2);
      hlepCosOpeningAngleGenDiff[i][j]->GetXaxis()->SetTitle("cos(#phi)");
      hlepCosOpeningAngleGenDiff[i][j]->Sumw2();
      
      htopCosThetaGenDiff[i][j] = new TH1D(Form("%s_htopCosThetaGenDiff_%s",prefix,suffix),Form("%s_topCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      htopCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}}) (reco-gen)");
      htopCosThetaGenDiff[i][j]->Sumw2();

      hlepCosThetaGenDiff[i][j] = new TH1D(Form("%s_hlepCosThetaGenDiff_%s",prefix,suffix),Form("%s_lepCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      hlepCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t}) (reco-gen)");
      hlepCosThetaGenDiff[i][j]->Sumw2();

      hlepPlusCosThetaGenDiff[i][j] = new TH1D(Form("%s_hlepPlusCosThetaGenDiff_%s",prefix,suffix),Form("%s_lepPlusCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      hlepPlusCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) (reco-gen)");
      hlepPlusCosThetaGenDiff[i][j]->Sumw2();
      
      hlepMinusCosThetaGenDiff[i][j] = new TH1D(Form("%s_hlepMinusCosThetaGenDiff_%s",prefix,suffix),Form("%s_lepMinusCosThetaGenDiff_%s",prefix,suffix),80,-2,2);
      hlepMinusCosThetaGenDiff[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}}) (reco-gen)");
      hlepMinusCosThetaGenDiff[i][j]->Sumw2();
            
      hpseudorapiditydiffGenDiff[i][j] = new TH1D(Form("%s_hpseudorapiditydiffGenDiff_%s",prefix,suffix),Form("%s_pseudorapiditydiffGenDiff_%s",prefix,suffix),80,-4,4);
      hpseudorapiditydiffGenDiff[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}| (reco-gen)");
      hpseudorapiditydiffGenDiff[i][j]->Sumw2();
      
      hrapiditydiffGenDiff[i][j] = new TH1D(Form("%s_hrapiditydiffGenDiff_%s",prefix,suffix),Form("%s_rapiditydiffGenDiff_%s",prefix,suffix),80,-4,4);
      hrapiditydiffGenDiff[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}}) (reco-gen)");
      hrapiditydiffGenDiff[i][j]->Sumw2();
      
      hrapiditydiffMarcoGenDiff[i][j] = new TH1D(Form("%s_hrapiditydiffMarcoGenDiff_%s",prefix,suffix),Form("%s_rapiditydiffMarcoGenDiff_%s",prefix,suffix),80,-4,4);
      hrapiditydiffMarcoGenDiff[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}| (reco-gen)");
      hrapiditydiffMarcoGenDiff[i][j]->Sumw2();
            
      


      
 
      hlepChargeAsym[i][j] = new TH1D(Form("%s_hlepChargeAsym_%s",prefix,suffix),Form("%s_lepChargeAsym_%s",prefix,suffix),80,-4,4);
      hlepChargeAsym[i][j]->GetXaxis()->SetTitle(" |#eta_{l^{+}}| - |#eta_{l^{-}}| ");
      hlepChargeAsym[i][j]->Sumw2();
 
      hlepRapDiff[i][j] = new TH1D(Form("%s_hlepRapDiff_%s",prefix,suffix),Form("%s_lepRapDiff_%s",prefix,suffix),100,-4,4);
      hlepRapDiff[i][j]->GetXaxis()->SetTitle(" #eta_{l^{+}} - #eta_{l^{-}} ");
      hlepRapDiff[i][j]->Sumw2();

      hlepAngleBetween[i][j] = new TH1D(Form("%s_hlepAngleBetween_%s",prefix,suffix),Form("%s_lepAngleBetween_%s",prefix,suffix),80,-1,1);
      hlepAngleBetween[i][j]->GetXaxis()->SetTitle("cos( #alpha_{l^{+}l^{-}}^{lab}) ");
      hlepAngleBetween[i][j]->Sumw2();
      
      hlepAngleBetweenCMS[i][j] = new TH1D(Form("%s_hlepAngleBetweenCMS_%s",prefix,suffix),Form("%s_lepAngleBetweenCMS_%s",prefix,suffix),80,-1,1);
      hlepAngleBetweenCMS[i][j]->GetXaxis()->SetTitle("cos( #alpha_{l^{+}l^{-}}^{t#bar{t}}) ");
      hlepAngleBetweenCMS[i][j]->Sumw2();


      hpseudorapiditydiff2[i][j] = new TH1D(Form("%s_hpseudorapiditydiff2_%s",prefix,suffix),Form("%s_pseudorapiditydiff2_%s",prefix,suffix),100,-6,6);
      hpseudorapiditydiff2[i][j]->GetXaxis()->SetTitle("#eta_{t} - #eta_{#bar{t}}");
      hpseudorapiditydiff2[i][j]->Sumw2();
      
      hrapiditydiff2[i][j] = new TH1D(Form("%s_hrapiditydiff2_%s",prefix,suffix),Form("%s_rapiditydiff2_%s",prefix,suffix),100,-4,4);
      hrapiditydiff2[i][j]->GetXaxis()->SetTitle("y_{t}-y_{#bar{t}} ");
      hrapiditydiff2[i][j]->Sumw2();

      hlepPlusCosTheta[i][j] = new TH1D(Form("%s_hlepPlusCosTheta_%s",prefix,suffix),Form("%s_lepPlusCosTheta_%s",prefix,suffix),80,-1,1);
      hlepPlusCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t})");
      hlepPlusCosTheta[i][j]->Sumw2();
      
      hlepMinusCosTheta[i][j] = new TH1D(Form("%s_hlepMinusCosTheta_%s",prefix,suffix),Form("%s_lepMinusCosTheta_%s",prefix,suffix),80,-1,1);
      hlepMinusCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#bar{t}})");
      hlepMinusCosTheta[i][j]->Sumw2();


      hjetAzimAsym[i][j] = new TH1D(Form("%s_hjetAzimAsym_%s",prefix,suffix),Form("%s_jetAzimAsym_%s",prefix,suffix),80,-1,1);
      hjetAzimAsym[i][j]->GetXaxis()->SetTitle("cos(#Delta #phi_{j1j2})");
      hjetAzimAsym[i][j]->Sumw2();
      
      hjetRapDiff[i][j] = new TH1D(Form("%s_hjetRapDiff_%s",prefix,suffix),Form("%s_jetRapDiff_%s",prefix,suffix),100,-4,4);
      hjetRapDiff[i][j]->GetXaxis()->SetTitle(" #eta_{j1} - #eta_{j2} ");
      hjetRapDiff[i][j]->Sumw2();
      
      hjetAngleBetween[i][j] = new TH1D(Form("%s_hjetAngleBetween_%s",prefix,suffix),Form("%s_jetAngleBetween_%s",prefix,suffix),80,-1,1);
      hjetAngleBetween[i][j]->GetXaxis()->SetTitle("cos( #alpha_{j1j2}^{lab}) ");
      hjetAngleBetween[i][j]->Sumw2();
      
      hjetAngleBetweenCMS[i][j] = new TH1D(Form("%s_hjetAngleBetweenCMS_%s",prefix,suffix),Form("%s_jetAngleBetweenCMS_%s",prefix,suffix),80,-1,1);
      hjetAngleBetweenCMS[i][j]->GetXaxis()->SetTitle("cos( #alpha_{j1j2}^{t#bar{t}}) ");
      hjetAngleBetweenCMS[i][j]->Sumw2();

      hlepPhi[i][j] = new TH1D(Form("%s_hlepPhi_%s",prefix,suffix),Form("%s_lepPhi_%s",prefix,suffix),80,-pi,pi);
      hlepPhi[i][j]->GetXaxis()->SetTitle("Lepton #phi");
      hlepPhi[i][j]->Sumw2();
            
      hlepPlusPhi[i][j] = new TH1D(Form("%s_hlepPlusPhi_%s",prefix,suffix),Form("%s_lepPlusPhi_%s",prefix,suffix),80,-pi,pi);
      hlepPlusPhi[i][j]->GetXaxis()->SetTitle("#phi_{l^{+}}");
      hlepPlusPhi[i][j]->Sumw2();
            
      hlepMinusPhi[i][j] = new TH1D(Form("%s_hlepMinusPhi_%s",prefix,suffix),Form("%s_lepMinusPhi_%s",prefix,suffix),80,-pi,pi);
      hlepMinusPhi[i][j]->GetXaxis()->SetTitle("#phi_{l^{-}}");
      hlepMinusPhi[i][j]->Sumw2();

      hjetPhi[i][j] = new TH1D(Form("%s_hjetPhi_%s",prefix,suffix),Form("%s_jetPhi_%s",prefix,suffix),80,-pi,pi);
      hjetPhi[i][j]->GetXaxis()->SetTitle("Jet #phi");
      hjetPhi[i][j]->Sumw2();

      hlepPlusEta[i][j] = new TH1D(Form("%s_hlepPlusEta_%s",prefix,suffix),Form("%s_lepPlusEta_%s",prefix,suffix),60,-3.,3.);
      hlepPlusEta[i][j]->GetXaxis()->SetTitle("#eta_{l^{+}}");
      hlepPlusEta[i][j]->Sumw2();

      hlepMinusEta[i][j] = new TH1D(Form("%s_hlepMinusEta_%s",prefix,suffix),Form("%s_lepMinusEta_%s",prefix,suffix),60,-3.,3.);
      hlepMinusEta[i][j]->GetXaxis()->SetTitle("#eta_{l^{-}}");
      hlepMinusEta[i][j]->Sumw2();

      hlepPlusPt[i][j] = new TH1D(Form("%s_hlepPlusPt_%s",prefix,suffix),Form("%s_lepPlusPt_%s",prefix,suffix),60,0.,240.);
      hlepPlusPt[i][j]->GetXaxis()->SetTitle("l^{+} p_{T} (GeV/c)");
      hlepPlusPt[i][j]->Sumw2();
      
      hlepMinusPt[i][j] = new TH1D(Form("%s_hlepMinusPt_%s",prefix,suffix),Form("%s_lepMinusPt_%s",prefix,suffix),60,0.,240.);
      hlepMinusPt[i][j]->GetXaxis()->SetTitle("l^{-} p_{T} (GeV/c)");
      hlepMinusPt[i][j]->Sumw2();








      hlepAzimAsym[i][j] = new TH1D(Form("%s_hlepAzimAsym_%s",prefix,suffix),Form("%s_lepAzimAsym_%s",prefix,suffix),80,-pi,pi);
      hlepAzimAsym[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym[i][j]->Sumw2();

      hlepAzimAsym2[i][j] = new TH1D(Form("%s_hlepAzimAsym2_%s",prefix,suffix),Form("%s_lepAzimAsym2_%s",prefix,suffix),80,0,pi);
      hlepAzimAsym2[i][j]->GetXaxis()->SetTitle("#Delta #phi_{l^{+}l^{-}}");
      hlepAzimAsym2[i][j]->Sumw2();

      htopSpinCorr[i][j] = new TH1D(Form("%s_htopSpinCorr_%s",prefix,suffix),Form("%s_topSpinCorr_%s",prefix,suffix),80,-1,1);
      htopSpinCorr[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{t}) #times cos(#theta_{l^{-}}^{#bar{t}} )");
      htopSpinCorr[i][j]->Sumw2();

      hlepCosOpeningAngle[i][j] = new TH1D(Form("%s_hlepCosOpeningAngle_%s",prefix,suffix),Form("%s_lepCosOpeningAngle_%s",prefix,suffix),80,-1,1);
      hlepCosOpeningAngle[i][j]->GetXaxis()->SetTitle("cos(#phi)");
      hlepCosOpeningAngle[i][j]->Sumw2();
      
      htopCosTheta[i][j] = new TH1D(Form("%s_htopCosTheta_%s",prefix,suffix),Form("%s_topCosTheta_%s",prefix,suffix),80,-1,1);
      htopCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{t}^{t#bar{t}})");
      htopCosTheta[i][j]->Sumw2();

      hpseudorapiditydiff[i][j] = new TH1D(Form("%s_hpseudorapiditydiff_%s",prefix,suffix),Form("%s_pseudorapiditydiff_%s",prefix,suffix),80,-4,4);
      hpseudorapiditydiff[i][j]->GetXaxis()->SetTitle("|#eta_{t}| - |#eta_{#bar{t}}|");
      hpseudorapiditydiff[i][j]->Sumw2();
      
      hrapiditydiff[i][j] = new TH1D(Form("%s_hrapiditydiff_%s",prefix,suffix),Form("%s_rapiditydiff_%s",prefix,suffix),80,-4,4);
      hrapiditydiff[i][j]->GetXaxis()->SetTitle("(y_{t}-y_{#bar{t}}) #times (y_{t}+y_{#bar{t}})");
      hrapiditydiff[i][j]->Sumw2();
      
      hrapiditydiffMarco[i][j] = new TH1D(Form("%s_hrapiditydiffMarco_%s",prefix,suffix),Form("%s_rapiditydiffMarco_%s",prefix,suffix),80,-4,4);
      hrapiditydiffMarco[i][j]->GetXaxis()->SetTitle("|y_{t}|-|y_{#bar{t}}|");
      hrapiditydiffMarco[i][j]->Sumw2();

      hlepCosTheta[i][j] = new TH1D(Form("%s_hlepCosTheta_%s",prefix,suffix),Form("%s_lepCosTheta_%s",prefix,suffix),80,-1,1);
      hlepCosTheta[i][j]->GetXaxis()->SetTitle("cos(#theta_{l}^{t})");
      hlepCosTheta[i][j]->Sumw2();

      httpT[i][j] = new TH1D(Form("%s_httpT_%s",prefix,suffix),Form("%s_ttpT_%s",prefix,suffix),101,-4.,400.);
      httpT[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} estimate (GeV/c)");
      httpT[i][j]->Sumw2();

      httpT_2d[i][j] = new TH2D(Form("%s_httpT2d_%s",prefix,suffix),Form("%s_ttpT2d_%s",prefix,suffix),101,-4.,400., 101,-4.,400.);
      httpT_2d[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} gen (GeV/c)");
      httpT_2d[i][j]->GetYaxis()->SetTitle("t#bar{t} p_{T} estimate (GeV/c)");
      httpT_2d[i][j]->Sumw2();

      htopP_2d[i][j] = new TH2D(Form("%s_htopP2d_%s",prefix,suffix),Form("%s_topP2d_%s",prefix,suffix),100,0.,800., 100,0.,800.);
      htopP_2d[i][j]->GetXaxis()->SetTitle("top CM momentum gen (GeV/c)");
      htopP_2d[i][j]->GetYaxis()->SetTitle("top CM momentum estimate (GeV/c)");
      htopP_2d[i][j]->Sumw2();

      httMass[i][j] = new TH1D(Form("%s_httMass_%s",prefix,suffix),Form("%s_ttMass_%s",prefix,suffix),100,100.,1100.);
      httMass[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} estimate (GeV/c^{2})");
      httMass[i][j]->Sumw2();
      
      httMass_pull[i][j] = new TH1D(Form("%s_httMasspull_%s",prefix,suffix),Form("%s_ttMasspull_%s",prefix,suffix),200,-3,3);
      httMass_pull[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} (reco-gen)/gen");
      httMass_pull[i][j]->Sumw2();

      httMass_diff[i][j] = new TH1D(Form("%s_httMassdiff_%s",prefix,suffix),Form("%s_ttMassdiff_%s",prefix,suffix),100,-300,300);
      httMass_diff[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} (reco-gen)");
      httMass_diff[i][j]->Sumw2();

      htopMass_diff_plus[i][j] = new TH1D(Form("%s_htopMassdiff_plus_%s",prefix,suffix),Form("%s_topMassdiff_plus_%s",prefix,suffix),150,-150,150);
      htopMass_diff_plus[i][j]->GetXaxis()->SetTitle("M_{t} (reco-gen)");
      htopMass_diff_plus[i][j]->Sumw2();

      htopMass_diff_minus[i][j] = new TH1D(Form("%s_htopMassdiff_minus_%s",prefix,suffix),Form("%s_topMassdiff_minus_%s",prefix,suffix),150,-150,150);
      htopMass_diff_minus[i][j]->GetXaxis()->SetTitle("M_{#bar{t}} (reco-gen)");
      htopMass_diff_minus[i][j]->Sumw2();

      htopPCM_diff_plus[i][j] = new TH1D(Form("%s_htopPCMdiff_plus_%s",prefix,suffix),Form("%s_topPCMdiff_plus_%s",prefix,suffix),80,-300,300);
      htopPCM_diff_plus[i][j]->GetXaxis()->SetTitle("top momentum in CM (reco-gen)");
      htopPCM_diff_plus[i][j]->Sumw2();

      htopPCM_diff_minus[i][j] = new TH1D(Form("%s_htopPCMdiff_minus_%s",prefix,suffix),Form("%s_topPCMdiff_minus_%s",prefix,suffix),80,-300,300);
      htopPCM_diff_minus[i][j]->GetXaxis()->SetTitle("#bar{t} momentum in CM (reco-gen)");
      htopPCM_diff_minus[i][j]->Sumw2();

      httMass_2d[i][j] = new TH2D(Form("%s_httMass2d_%s",prefix,suffix),Form("%s_ttMass2d_%s",prefix,suffix),100,100.,1100., 100,100.,1100.);
      httMass_2d[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} gen (GeV/c^{2})");
      httMass_2d[i][j]->GetYaxis()->SetTitle("M_{t#bar{t}} estimate (GeV/c^{2})");
      httMass_2d[i][j]->Sumw2();

      htopMass[i][j] = new TH1D(Form("%s_htopMass_%s",prefix,suffix),Form("%s_topMass_%s",prefix,suffix),102,98,302);
      htopMass[i][j]->GetXaxis()->SetTitle("M_{t} estimate (GeV/c^{2})");
      htopMass[i][j]->Sumw2();
*/


      httpT_gen[i][j] = new TH1D(Form("%s_httpTGen_%s",prefix,suffix),Form("%s_ttpTGen_%s",prefix,suffix),101,-4.,400.);
      httpT_gen[i][j]->GetXaxis()->SetTitle("t#bar{t} p_{T} gen (GeV/c)");
      httpT_gen[i][j]->Sumw2();

      httMass_gen[i][j] = new TH1D(Form("%s_httMassGen_%s",prefix,suffix),Form("%s_ttMassGen_%s",prefix,suffix),100,100.,1100.);
      httMass_gen[i][j]->GetXaxis()->SetTitle("M_{t#bar{t}} gen (GeV/c^{2})");
      httMass_gen[i][j]->Sumw2();

      htopMass_plus_gen[i][j] = new TH1D(Form("%s_htopMass_plusGen_%s",prefix,suffix),Form("%s_topMass_plusGen_%s",prefix,suffix),100,100.,250.);
      htopMass_plus_gen[i][j]->GetXaxis()->SetTitle("M_{t} gen (GeV/c^{2})");
      htopMass_plus_gen[i][j]->Sumw2();

      htopMass_minus_gen[i][j] = new TH1D(Form("%s_htopMass_minusGen_%s",prefix,suffix),Form("%s_topMass_minusGen_%s",prefix,suffix),100,100.,250.);
      htopMass_minus_gen[i][j]->GetXaxis()->SetTitle("M_{#bar{t}} gen (GeV/c^{2})");
      htopMass_minus_gen[i][j]->Sumw2();

/*
      httRapidity[i][j] = new TH1D(Form("%s_httRapidity_%s",prefix,suffix),Form("%s_ttRapidity_%s",prefix,suffix),120,-6,6);
      httRapidity[i][j]->GetXaxis()->SetTitle("y_{t} + y_{#bar{t}} estimate");
      httRapidity[i][j]->Sumw2();
      
      httRapidity2[i][j] = new TH1D(Form("%s_httRapidity2_%s",prefix,suffix),Form("%s_ttRapidity2_%s",prefix,suffix),120,-6,6);
      httRapidity2[i][j]->GetXaxis()->SetTitle("y_{t#bar{t}} estimate");
      httRapidity2[i][j]->Sumw2();
	
      hmassltb[i][j] =  new TH1D(Form("%s_hmassltb_%s",prefix,suffix),Form("%s_massltb_%s",prefix,suffix),75,0.,510.);
      hmassltb[i][j]->GetXaxis()->SetTitle("M_{l1b1} (GeV/c^{2})");
      hmassltb[i][j]->Sumw2();

      hmassllb[i][j] =  new TH1D(Form("%s_hmassllb_%s",prefix,suffix),Form("%s_massllb_%s",prefix,suffix),75,0.,510.);
      hmassllb[i][j]->GetXaxis()->SetTitle("M_{l2b2} (GeV/c^{2})");
      hmassllb[i][j]->Sumw2();
      
      hmassltb1Dmasscut[i][j] =  new TH1D(Form("%s_hmassltb1Dmasscut_%s",prefix,suffix),Form("%s_massltb1Dmasscut_%s",prefix,suffix),75,0.,510.);
      hmassltb1Dmasscut[i][j]->GetXaxis()->SetTitle("M_{l1b1} (GeV/c^{2}) for M_{l2b2} > 170 GeV/c^{2}");
      hmassltb1Dmasscut[i][j]->Sumw2();

      hmassllb1Dmasscut[i][j] =  new TH1D(Form("%s_hmassllb1Dmasscut_%s",prefix,suffix),Form("%s_massllb1Dmasscut_%s",prefix,suffix),75,0.,510.);
      hmassllb1Dmasscut[i][j]->GetXaxis()->SetTitle("M_{l2b2} (GeV/c^{2}) for M_{l1b1} > 170 GeV/c^{2}");
      hmassllb1Dmasscut[i][j]->Sumw2();

      htheSumJetPt[i][j] =  new TH1D(Form("%s_theSumJetPt_%s",prefix,suffix),Form("%s_theSumJetPt_%s",prefix,suffix),100,0.,600.);
      htheSumJetPt[i][j]->GetXaxis()->SetTitle("Sum of jet p_{T} (GeV/c)");
      htheSumJetPt[i][j]->Sumw2();

      htheSumBtagJetPt[i][j] =  new TH1D(Form("%s_theSumBtagJetPt_%s",prefix,suffix),Form("%s_theSumBtagJetPt_%s",prefix,suffix),100,0.,600.);
      htheSumBtagJetPt[i][j]->GetXaxis()->SetTitle("Sum of b tagged jet p_{T} (GeV/c)");
      htheSumBtagJetPt[i][j]->Sumw2();

      hthefirstJetPt[i][j] =  new TH1D(Form("%s_thefirstJetPt_%s",prefix,suffix),Form("%s_thefirstJetPt_%s",prefix,suffix),60,0.,360.);
      hthefirstJetPt[i][j]->GetXaxis()->SetTitle("First jet p_{T} (GeV/c)");
      hthefirstJetPt[i][j]->Sumw2();

      hthesecondJetPt[i][j] =  new TH1D(Form("%s_thesecondJetPt_%s",prefix,suffix),Form("%s_thesecondJetPt_%s",prefix,suffix),60,0.,360.);
      hthesecondJetPt[i][j]->GetXaxis()->SetTitle("Second jet p_{T} (GeV/c)");
      hthesecondJetPt[i][j]->Sumw2();
      
      htheleadinglepPt[i][j] = new TH1D(Form("%s_htheleadinglepPt_%s",prefix,suffix),Form("%s_theleadinglepPt_%s",prefix,suffix),60,0.,240.);
	  htheleadinglepPt[i][j]->GetXaxis()->SetTitle("Leading lepton p_{T} (GeV/c)");
      htheleadinglepPt[i][j]->Sumw2();

      hthesecondlepPt[i][j] = new TH1D(Form("%s_hthesecondlepPt_%s",prefix,suffix),Form("%s_thesecondlepPt_%s",prefix,suffix),60,0.,240.);
	  hthesecondlepPt[i][j]->GetXaxis()->SetTitle("Second lepton p_{T} (GeV/c)");
      hthesecondlepPt[i][j]->Sumw2();

      htheSumLepPt[i][j] =  new TH1D(Form("%s_theSumLepPt_%s",prefix,suffix),Form("%s_theSumLepPt_%s",prefix,suffix),120,0.,480.);
      htheSumLepPt[i][j]->GetXaxis()->SetTitle("Sum of lepton p_{T} (GeV/c)");
      htheSumLepPt[i][j]->Sumw2();

      hlepEta[i][j] = new TH1D(Form("%s_hlepEta_%s",prefix,suffix),Form("%s_lepEta_%s",prefix,suffix),60,-3.,3.);
      hlepEta[i][j]->GetXaxis()->SetTitle("Lepton #eta");
      hlepEta[i][j]->Sumw2();

      hlepPt[i][j] = new TH1D(Form("%s_hlepPt_%s",prefix,suffix),Form("%s_lepPt_%s",prefix,suffix),60,0.,240.);
      hlepPt[i][j]->GetXaxis()->SetTitle("Lepton p_{T} (GeV/c)");
      hlepPt[i][j]->Sumw2();
      
      hjetPt[i][j] = new TH1D(Form("%s_hjetPt_%s",prefix,suffix),Form("%s_jetPt_%s",prefix,suffix),60,0.,360.);
      hjetPt[i][j]->GetXaxis()->SetTitle("Jet p_{T} (GeV/c)");
      hjetPt[i][j]->Sumw2();

      hjetEta[i][j] = new TH1D(Form("%s_hjetEta_%s",prefix,suffix),Form("%s_jetEta_%s",prefix,suffix),60,-3.,3.);
      hjetEta[i][j]->GetXaxis()->SetTitle("Jet #eta");
      hjetEta[i][j]->Sumw2();
      
     
      hMET[i][j] = new TH1D(Form("%s_hMET_%s",prefix,suffix),Form("%s_MET_%s",prefix,suffix),100,0.,400.);
      hMET[i][j]->GetXaxis()->SetTitle("E_{T}^{miss} (GeV)");
      hMET[i][j]->Sumw2();

      hdMET[i][j] = new TH1D(Form("%s_hdMET_%s",prefix,suffix),Form("%s_dMET_%s",prefix,suffix),100,0.,400.);
      hdMET[i][j]->GetXaxis()->SetTitle("#DeltaMET (GeV)");
      hdMET[i][j]->Sumw2();

       
      hmasslb_2d[i][j] = new TH2D(Form("%s_hmasslb2d_%s",prefix,suffix),  Form("%s_masslb2d_%s" ,prefix,suffix),100,0,500,100,0,500);
      hmasslb_2d[i][j]->GetXaxis()->SetTitle("M_{l2b2}(GeV/c^{2})");
      hmasslb_2d[i][j]->GetYaxis()->SetTitle("M_{l1b1}(GeV/c^{2})");
      hmasslb_2d[i][j]->Sumw2();



      // generator level distributions

      httMassGluongenp[i][j] = new TH1D(Form("%s_httMassGluongenp_%s",prefix,suffix),Form("%s_ttMassGluongenp_%s",prefix,suffix),200, 0 , 1000);
      httMassGluongenp[i][j]->GetXaxis()->SetTitle("ttMassGluongenp ");
      httMassGluongenp[i][j]->Sumw2();
      
      httMassQuarkgenp[i][j] = new TH1D(Form("%s_httMassQuarkgenp_%s",prefix,suffix),Form("%s_ttMassQuarkgenp_%s",prefix,suffix),200, 0 , 1000);
      httMassQuarkgenp[i][j]->GetXaxis()->SetTitle("ttMassQuarkgenp ");
      httMassQuarkgenp[i][j]->Sumw2();

      httRapidityGluongenp[i][j] = new TH1D(Form("%s_httRapidityGluongenp_%s",prefix,suffix),Form("%s_ttRapidityGluongenp_%s",prefix,suffix),200, -10, 10);
      httRapidityGluongenp[i][j]->GetXaxis()->SetTitle("ttRapidityGluongenp ");
      httRapidityGluongenp[i][j]->Sumw2();
      
      httRapidityQuarkgenp[i][j] = new TH1D(Form("%s_httRapidityQuarkgenp_%s",prefix,suffix),Form("%s_ttRapidityQuarkgenp_%s",prefix,suffix),200,-10, 10);
      httRapidityQuarkgenp[i][j]->GetXaxis()->SetTitle("ttRapidityQuarkgenp ");
      httRapidityQuarkgenp[i][j]->Sumw2();
*/


      //daughter lepton angle in tau rest frame to check if MC is correctly using the tau polarisation
      hlepPlusCosThetaTau_gen[i][j] = new TH1D(Form("%s_hlepPlusCosThetaTauGen_%s",prefix,suffix),Form("%s_lepPlusCosThetaTauGen_%s",prefix,suffix),80,-1,1);
      hlepPlusCosThetaTau_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{+}}^{#tau})");
      hlepPlusCosThetaTau_gen[i][j]->Sumw2();
      
      hlepMinusCosThetaTau_gen[i][j] = new TH1D(Form("%s_hlepMinusCosThetaTauGen_%s",prefix,suffix),Form("%s_lepMinusCosThetaTauGen_%s",prefix,suffix),80,-1,1);
      hlepMinusCosThetaTau_gen[i][j]->GetXaxis()->SetTitle("cos(#theta_{l^{-}}^{#tau})");
      hlepMinusCosThetaTau_gen[i][j]->Sumw2();

      //daughter lepton E/Emax in tau rest frame to check if MC is correctly using the tau polarisation
      hlepPlusxTau_gen[i][j] = new TH1D(Form("%s_hlepPlusxTauGen_%s",prefix,suffix),Form("%s_lepPlusxTauGen_%s",prefix,suffix),80,0,1);
      hlepPlusxTau_gen[i][j]->GetXaxis()->SetTitle("x");
      hlepPlusxTau_gen[i][j]->Sumw2();
      
      hlepMinusxTau_gen[i][j] = new TH1D(Form("%s_hlepMinusxTauGen_%s",prefix,suffix),Form("%s_lepMinusxTauGen_%s",prefix,suffix),80,0,1);
      hlepMinusxTau_gen[i][j]->GetXaxis()->SetTitle("x");
      hlepMinusxTau_gen[i][j]->Sumw2();

/*            
      //DYEst histos
      hdilMassWithMetDYEst[i][j] = new TH1D(Form("%s_hdilMassWithMetDYEst_%s",  prefix,suffix), "Di-lepton mass with MET for DY Estimation", 40, 0., 200.);
      hdilMassWithMetDYEst[i][j]->GetXaxis()->SetTitle("M_{ll}(GeV/c^{2}), with MET > 30 GeV cut");
      hdilMassWithMetDYEst[i][j]->Sumw2();

      hdilMassNoMetDYEst[i][j] = new TH1D(Form("%s_hdilMassNoMetDYEst_%s",  prefix,suffix), "Di-lepton mass without MET for DY Estimation", 40, 0., 200.);
      hdilMassNoMetDYEst[i][j]->GetXaxis()->SetTitle("M_{ll}(GeV/c^{2}), no MET cut");
      hdilMassNoMetDYEst[i][j]->Sumw2();

      hmetInDYEst[i][j] = new TH1D(Form("%s_hmetInDYEst_%s",  prefix,suffix), "MET in Z mass for DY Estimation", 40, 0., 200.);
      hmetInDYEst[i][j]->GetXaxis()->SetTitle("MET (GeV), inside Z mass window");
      hmetInDYEst[i][j]->Sumw2();

      hmetOutDYEst[i][j] = new TH1D(Form("%s_hmetOutDYEst_%s",  prefix,suffix), "MET outside Z mass for DY Estimation", 40, 0., 200.);
      hmetOutDYEst[i][j]->GetXaxis()->SetTitle("MET (GeV), outside Z mass window");
      hmetOutDYEst[i][j]->Sumw2();
      */
    }
    
  }

 

}
