void doAll_nocuts(TString outputDir = "results")
{
    //gSystem->Load("/home/users/yanjuntu/MiniFWlib/libMiniFWLite_v5.28.00.so");
    // gSystem->Load("/home/users/yanjuntu/MiniFWlib/libMiniFWLite_5.27.06b-cms10.so");
    gSystem->Load("/home/users/yanjuntu/MiniFWlib/libMiniFWLite_CMSSW_5_3_2_patch4_V05-03-13.so");
    //gSystem->Load("/nfs-3/userdata/yanjuntu/lhapdf/lib/libLHAPDF.so");
    //gSystem->Load("/nfs-6/userdata/yanjuntu/LHAPDF/lib/libLHAPDF.so");

    /*
    gSystem->AddIncludePath(" -w -I../CORE/topmass -I/nfs-6/userdata/yanjuntu/LHAPDF/include");
    gROOT->ProcessLine(".L ../CORE/topmass/ttdilepsolve.cpp+");
    gROOT->ProcessLine(".L ../CORE/utilities.cc+");
    gROOT->ProcessLine(".L ../CORE/trackSelections.cc+");
    gROOT->ProcessLine(".L ../CORE/eventSelections.cc+");
    gROOT->ProcessLine(".L ../CORE/MITConversionUtilities.cc+");
    gROOT->ProcessLine(".L ../CORE/muonSelections.cc+");
    gROOT->ProcessLine(".L ../CORE/electronSelectionsParameters.cc+");
    gROOT->ProcessLine(".L ../CORE/electronSelections.cc+");
    gROOT->ProcessLine(".L ../CORE/metSelections.cc+");
    gROOT->ProcessLine(".L ../CORE/SimpleFakeRate.cc+");
    gROOT->ProcessLine(".L ../CORE/MT2/MT2.cc+");
    gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");
    gROOT->ProcessLine(".L ../CORE/susySelections.cc+");
    //gROOT->ProcessLine(".L ../CORE/mcSUSYkfactor.cc+");
    gROOT->ProcessLine(".L ../CORE/triggerSuperModel.cc+");
    gROOT->ProcessLine(".L ../CORE/jetSelections.cc+");
    gROOT->ProcessLine(".L ../CORE/ttbarSelections.cc+");
    gROOT->ProcessLine(".L ../CORE/triggerUtils.cc+");
    */

    gROOT->ProcessLine(".L ../../../CORE/CMS2.cc+");
    gROOT->ProcessLine(".L ../../../CORE/mcSelections.cc+");

    gSystem->CompileMacro("histtools.C", "++k", "libhisttools");
    gSystem->CompileMacro("topAFB_looper.C", "++k", "libtopAFB_looper");


    float lumiToNormalizeTo   = 4980.0e-3;

    topAFB_looper *baby = new topAFB_looper();

    vector<TString> v_baseCuts;
    v_baseCuts.push_back("applyTopPtWeighting");
    v_baseCuts.push_back("weighttaudecay");

    vector<TString> v_Cuts = v_baseCuts;


    if (true)
    {
        TChain  *ch_ttbar = new TChain("Events");

        //cout << "Doing the MadGraph ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

        //cout << "Doing the powheg ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

        //cout << "Doing the powheg tauola ttbar sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/TT_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

        cout << "Doing the Fall11 MC@NLO ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29_fix_dilepton/merged*root");


        //for samples with no negative weights
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",lumiToNormalizeTo*154./157.5);
        //for mc@NLO
        baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * (154. / 157.5) * (190.41256 / 147.4)); //must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights

        hist::color("ttdil", kGreen);
        //delete ch_ttbar;
    }
    if (false)
    {
        TChain  *ch_ttbar = new TChain("Events");
        //  cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTJets_TuneZ2_mass178_5_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v3/V04-02-29/merged*root");
        //    baby->ScanChain(ch_ttbar, v_Cuts, "ttdil_m178",lumiToNormalizeTo*154./157.5);
        //    hist::color("ttdil_m178", kGreen);
        //  cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTJets_TuneZ2_mass166_5_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v3/V04-02-29/merged*root");
        //    baby->ScanChain(ch_ttbar, v_Cuts, "ttdil_m166",lumiToNormalizeTo*154./157.5);
        //    hist::color("ttdil_m166", kGreen);
        //   cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTjets_TuneZ2_matchingup_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
        // cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-4/userdata/cms2/TTjets_TuneZ2_matchingdown_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");
        //cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2011/Fall11MC/TTjets_TuneZ2_scaleup_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v1/V04-02-29/merged*root");
        //cout << "Doing the MadGraph ttbar sys sample" << endl; ch_ttbar->Add("/nfs-6/userdata/cms2/TTjets_TuneZ2_scaledown_7TeV-madgraph-tauola_Fall11-PU_S6_START42_V14B-v2/V04-02-29/merged*root");

        cout << "Doing the scale down mc@nlo ttbar sys sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/FastSim_mcatnlo-private-SC-fac-0.5-ren-0.5_LHE2EDM_v1_RWTH_0614_fhohle-TT_mcatnlo_private_fac-0.5-ren-0.5_withSC_FastSim_coherent_Summer11FullSim_FastSimJetFixTest_AOD_0616-98893edf89c918d9b6d63453d739e0c0/VB04-02-29_Fastsim/merged*root");
        //cout << "Doing the scale up mc@nlo ttbar sys sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/FastSim_mcatnlo-private-SC-fac-2-ren-2_LHE2EDM_v1_RWTH_0614_fhohle-TT_mcatnlo_private_fac-2-ren-2_withSC_FastSim_coherent_Summer11FullSim_FastSimJetFixTest_AOD_0616-98893edf89c918d9b6d63453d739e0c0/VB04-02-29_Fastsim/merged*root");
        baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 154. / 157.5);
        hist::color("ttdil", kGreen);

    }


    TString cutstring = "";
    cutstring  = outputDir + "/hist_noCuts.root";

    cout << "Saving histograms to: " << cutstring << endl;

    hist::saveHist(cutstring.Data());
    cout << "done saving" << endl;

    hist::deleteHistos();



}

