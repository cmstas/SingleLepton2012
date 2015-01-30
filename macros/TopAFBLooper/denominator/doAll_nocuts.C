void doAll_nocuts(TString outputDir = "results")
{
    gSystem->Load("../../../Tools/MiniFWLite/libMiniFWLite.so");

    gROOT->ProcessLine(".L ../../../CORE/CMS2.cc+");
    gROOT->ProcessLine(".L ../../../CORE/mcSelections.cc+");

    gSystem->CompileMacro("histtools.C", "++k", "libhisttools");
    gSystem->CompileMacro("topAFB_looper.C", "++k", "libtopAFB_looper");


    float lumiToNormalizeTo   = 19.5;

    topAFB_looper *baby = new topAFB_looper();

    vector<TString> v_baseCuts;
    //v_baseCuts.push_back("applyTopPtWeighting");
    //v_baseCuts.push_back("weighttaudecay");

    vector<TString> v_Cuts = v_baseCuts;

    //need to put 8 TeV cms2 ttbar ntuples here
    if (true)
    {
        TChain  *ch_ttbar = new TChain("Events");

        //cout << "Doing the MadGraph ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TTJets_TuneZ2_7TeV-madgraph-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

        //cout << "Doing the powheg ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TTTo2L2Nu2B_7TeV-powheg-pythia6_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

        //cout << "Doing the powheg tauola ttbar sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2011/Summer11MC/TT_TuneZ2_7TeV-powheg-tauola_Summer11-PU_S4_START42_V11-v1/V04-02-29/merged*root");

        //cout << "Doing the Fall11 MC@NLO ttbar sample" << endl; ch_ttbar->Add("/nfs-7/userdata/cms2/TT_TuneZ2_7TeV-mcatnlo_Fall11-PU_S6_START42_V14B-v1/V04-02-29_fix_dilepton/merged*root");
        cout << "Doing the 8TeV MC@NLO ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postprocessing_TT_mcatnlo_4/preprocessing/ntuple*.root");

        //cout << "Doing the 8TeV powheg v1 ttbar sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v1/V05-03-25/merged_ntuple*.root");
        //cout << "Doing the 8TeV powheg v2 ttbar sample" << endl; ch_ttbar->Add("/hadoop/cms/store/group/snt/papers2012/Summer12_53X_MC/TT_CT10_TuneZ2star_8TeV-powheg-tauola_Summer12_DR53X-PU_S10_START53_V7A-v2/V05-03-25/merged_ntuple*.root");

        //cout << "Doing the 8TeV MC@NLO mass down ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_TT_FullLept_mass169_5_8TeV-mcatnlo/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV MC@NLO mass up ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_TT_FullLept_mass175_5_8TeV-mcatnlo/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV MC@NLO scale down ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_TT_FullLept_scaledown_8TeV-mcatnlo/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV MC@NLO scale up ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_TT_FullLept_scaleup_8TeV-mcatnlo/preprocessing/ntuple*.root");

        //cout << "Doing the 8TeV ZprimeTtoTTU_M-220_TuneZ2star_8TeV-madgraph-tauola ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_ZprimeTtoTTU_M-220_TuneZ2star_8TeV-madgraph-tauola/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV HeavyGluonToTT_axial_M-2000_TuneZ2star_8TeV-madgraph-tauola ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_HeavyGluonToTT_axial_M-2000_TuneZ2star_8TeV-madgraph-tauola/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV HeavyGluonToTT_axial_M-200_TuneZ2star_8TeV-madgraph-tauola ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_HeavyGluonToTT_axial_M-200_TuneZ2star_8TeV-madgraph-tauola/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV HeavyGluonToTT_left_M-200_TuneZ2star_8TeV-madgraph-tauola ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_HeavyGluonToTT_left_M-200_TuneZ2star_8TeV-madgraph-tauola/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV HeavyGluonToTT_right_M-2000_TuneZ2star_8TeV-madgraph-tauola ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_HeavyGluonToTT_right_M-2000_TuneZ2star_8TeV-madgraph-tauola/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV HeavyGluonToTT_right_M-200_TuneZ2star_8TeV-madgraph-tauola ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_HeavyGluonToTT_right_M-200_TuneZ2star_8TeV-madgraph-tauola/preprocessing/ntuple*.root");
        //cout << "Doing the 8TeV TT_noCorr_8TeV-mcatnlo ttbar sample" << endl; ch_ttbar->Add("/home/users/linacre/CMSSW_5_3_2_patch5_V05-03-32/crabnew/postProcessing_TT_noCorr_8TeV-mcatnlo/preprocessing/ntuple*.root");


        //for samples with no negative weights
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil",lumiToNormalizeTo*154./157.5);
        //for mc@NLO
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * (154. / 157.5) * (190.41256 / 147.4)); //must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights
        baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000.*274.00756/211.1/32852589. ); //Apply the evt_scale1fb normalisation missing from the mcatnlo ntuples. Must mutliply by ratio of per-event xsec and PREP xsec to account for negative weights.

        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo );
        
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000.*274.00756/211.1/3110611. ); //mass down
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000.*274.00756/211.1/3217835. ); //mass up
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000.*274.00756/211.1/3138155. ); //scale down
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000.*274.00756/211.1/3128986. ); //scale up

        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000./956164. );    //ZprimeTtoTTU_M_220_TuneZ2star_8TeV_madgraph_tauola
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000./971034. );    //HeavyGluonToTT_axial_M_2000_TuneZ2star_8TeV_madgraph_tauola
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000./985684. );    //HeavyGluonToTT_axial_M_200_TuneZ2star_8TeV_madgraph_tauola
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000./992939. );    //HeavyGluonToTT_left_M_200_TuneZ2star_8TeV_madgraph_tauola
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000./993386. );    //HeavyGluonToTT_right_M_2000_TuneZ2star_8TeV_madgraph_tauola
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000./996380. );    //HeavyGluonToTT_right_M_200_TuneZ2star_8TeV_madgraph_tauola
        //baby->ScanChain(ch_ttbar, v_Cuts, "ttdil", lumiToNormalizeTo * 0.945 * 234000.*274.00756/211.1/24766326. );    //TT_noCorr_8TeV_mcatnlo


        hist::color("ttdil", kGreen);
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

