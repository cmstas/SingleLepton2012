{

  gSystem->Load("libGraf");
  gSystem->Load("libGpad");
  gSystem->Load("libTree");

  gSystem->SetAclicMode(TSystem::kDebug );
  gROOT->ProcessLine(".L tdrstyle_SUSY.C");
  setTDRStyle();
  gROOT->ForceStyle();

  string msg("* Welcome to ROOT v");
  msg += gROOT->GetVersion();
  msg += " *";
  string ast(msg.size(), '*');
  cout << ast << endl << msg << endl << ast << endl;

}
