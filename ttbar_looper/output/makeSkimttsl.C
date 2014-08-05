#include <fstream>
#include <sstream>
#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>

using namespace std;

void makeSkimttsl(string path = ".", string tag = "tt_mcatnlo_smallTree") {
  
  //--------------------------------------------------
  // cut for output files
  //--------------------------------------------------
  
  char* sel = "ntops == 2 && nleps == 1";
  cout << "Skimming with selection : "<<sel<<endl;

  //--------------------------------------------------
  // input and output file
  //--------------------------------------------------
  
  char* infilename = Form("%s/%s.root",path.c_str(),tag.c_str());
  char* outfilename = Form("%s/%s_ttsl.root",path.c_str(),tag.c_str());
  
  //--------------------------------------------------
  // cout stuff
  //--------------------------------------------------
  
  cout << "Reading in : " << infilename << endl;
  cout << "Writing to : " << outfilename << endl;
  cout << "Selection : " << sel << endl;
  
  //--------------------------------------------------
  // read input file, write to output files
  //--------------------------------------------------
  
  long long max_tree_size = 20000000000000000LL;
  TTree::SetMaxTreeSize(max_tree_size);
  
  TChain *chain = new TChain("t");
  chain->Add(infilename);
  
  cout << "Input tree has entries: " << chain->GetEntries() << endl;
  
  //-------------------
  // skim
  //-------------------
  
  TFile *outfile = TFile::Open(outfilename, "RECREATE");
  assert( outfile != 0 );
  TTree* outtree = chain->CopyTree( sel );
  cout << "Output tree has entries: " << outtree->GetEntries() << endl;
  outtree->Write();
  outfile->Close();

}
