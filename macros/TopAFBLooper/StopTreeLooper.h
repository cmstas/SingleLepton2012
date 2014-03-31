#ifndef STOPTREELOOPER_H
#define STOPTREELOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"

#include <iostream>
#include "Math/LorentzVector.h"
 
#include <cmath>
#include <map>

using namespace std;

class StopTreeLooper {

    public:
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float> > LorentzVector;

        StopTreeLooper();
        ~StopTreeLooper();
 
        void setOutFileName(string filename); 
        void loop(TChain *chain, TString name);

	//plotting
	void makeSIGPlots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string flav_tag ); 
	void makeCR1Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string flav_tag ); 
	void makeCR2Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string flav_tag_dl );
	void makeCR3Plots(float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string flav_tag_dl );
	void makeNJPlots( float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string flav_tag ); 
	void makeZPlots(  float evtweight, std::map<std::string, TH1F*> &h_1d, 
			   string tag_selection, string flav_tag );

    private:

	string m_outfilename_;
	// njets requirement
	int min_njets;
	//for phi corrected met
	float t1metphicorr;
	float t1metphicorrphi;

	//jets information
	int n_jets;
	int n_bjets;
	int n_ljets;
	vector<LorentzVector> jets;
	vector<float> btag;
	vector<float> sigma_jets;
	vector<int> mc;

	float pfcalo_metratio;
	float pfcalo_metdphi;


};

#endif
