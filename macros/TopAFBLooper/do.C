void do(char* path = "/nfs-6/userdata/stop/output_V00-02-32_2012", char* sample = "DY1to4Jtautau"){
  //------------------------------ 
  // load stuff
  //------------------------------ 

  gSystem->Load("libTree.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libEG.so");
  gSystem->Load("libMathCore.so");

  gSystem->Load("LHAPDF/lib/libLHAPDF.so");

  gSystem->Load("../../Tools/MiniFWLite/libMiniFWLite.so");

  //  gROOT->ProcessLine(".L ../../CORE/libCMS2NtupleMacrosCORE.so");
  gROOT->ProcessLine(".L libStopTreeLooper.so");

  StopTreeLooper *looper = new StopTreeLooper();

  //------------------------------ 
  // process sample
  //------------------------------ 

  TString sampletemp = sample;
  char* samplename = ( sampletemp.Contains("DY1to4J") ? "DY1to4Jtot" : sample );

  TChain *ch = new TChain("t");
  ch->Add(Form("%s/%s*.root", path, samplename));
  looper->setOutFileName(Form("output/%s_histos.root", sample));
  looper->loop(ch, sample);

  delete looper;
  delete ch;
}
