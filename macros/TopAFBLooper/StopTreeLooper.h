#ifndef STOPTREELOOPER_H
#define STOPTREELOOPER_H

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"

#include <iostream>
#include "Math/LorentzVector.h"

//#include "../CORE/topmass/ttdilepsolve.cpp" 
//#include "../CORE/topmass/getTopMassEstimate.icc" 
 
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
	void makettPlots(float evtweight, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d,
			   string tag_selection, string flav_tag ); 
	void makeAccPlots(float evtweight, std::map<std::string, TH1F*> &h_1d, std::map<std::string, TH2F*> &h_2d,
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

	void solvettbar();

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
	vector<LorentzVector> bjets;
	vector<LorentzVector> nonbjets;
	vector<LorentzVector> bcandidates;
	vector<float> btag;
	vector<float> sigma_jets;
	vector<int> mc;
	



	double lep_charge_asymmetry;
	double lep_azimuthal_asymmetry;
	double lep_azimuthal_asymmetry2;
	double top_rapiditydiff_cms;
	double top_pseudorapiditydiff_cms;
	double top_rapiditydiff_Marco;
	double top_costheta_cms;
	double lepPlus_costheta_cms;
	double lepMinus_costheta_cms;
	double top_spin_correlation;
	double lep_cos_opening_angle;
	double tt_mass;
	double ttRapidity2;
	double tt_pT;
	double top1_pt;
	double top2_pt;
	double top1_p_CM;
	double top2_p_CM;
	double top_rapiditydiffsigned_cms;

	double lep_charge_asymmetry_gen;
	double lep_azimuthal_asymmetry_gen;
	double lep_azimuthal_asymmetry2_gen;
	double top_rapiditydiff_cms_gen;
	double top_pseudorapiditydiff_cms_gen;
	double top_rapiditydiff_Marco_gen;
	double top_costheta_cms_gen;
	double lepPlus_costheta_cms_gen;
	double lepMinus_costheta_cms_gen;
	double top_spin_correlation_gen;
	double lep_cos_opening_angle_gen;
	double tt_mass_gen;
	double ttRapidity2_gen;
	double tt_pT_gen;
	double top1_pt_gen;
	double top2_pt_gen;



	TLorentzVector lepPlus;
	TLorentzVector lepMinus;
	TLorentzVector jet1;
	TLorentzVector jet2;
	TLorentzVector top1_p4;
	TLorentzVector top2_p4;
	TLorentzVector nusum;
	TLorentzVector cms;

	TLorentzVector lepPlus_gen;
	TLorentzVector lepMinus_gen;
	TLorentzVector topplus_genp_p4;
	TLorentzVector topminus_genp_p4;
	TLorentzVector cms_gen;
	TLorentzVector nuPlus_gen;
	TLorentzVector nuMinus_gen;


    vector <TLorentzVector> nu1_vecs;
    vector <TLorentzVector> nu2_vecs;
    vector <TLorentzVector> top1_vecs;
    vector <TLorentzVector> top2_vecs;
    vector <double> AMWT_weights;
    double m_top;
    double m_top_B;
    double closestDeltaMET_maxwcombo;
    double closestDeltaMET_othercombo;
    double closestDeltaMET_bestcombo;

    int imaxweight;
    bool closestApproach;

	float pfcalo_metratio;
	float pfcalo_metdphi;

	//ttdilepsolve *d_llsol;


};

#endif
