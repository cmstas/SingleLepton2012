#ifndef stopUtils_h
#define stopUtils_h

#include <algorithm>
#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <list>
#include <sstream>

#include "STOPT.h"

#include "Math/LorentzVector.h"
#include "Math/Point3D.h"
#include "TMath.h"
#include "TBranch.h"
#include "TTree.h"
#include "TH1F.h"
#include "TFile.h"
#include "TBits.h"

#include "../../Tools/BTagReshaping/BTagReshaping.h"
#include "../Core/MT2.h"
#include "../Core/mt2bl_bisect.h"
#include "../Core/mt2w_bisect.h"

#include "../../CORE/utilities.h"

using namespace std;

class Candidate : public TObject {
    public:
        float chi2, mt2w, mt2bl, mt2b;
        int j1, j2, bi, oi;
        float k1, k2;
        bool match;

        ClassDef(Candidate, 2)
};

typedef vector<Candidate> CANDIDATES;

float getThetaStar(LorentzVector& vec1, LorentzVector& vec2 );

float passLRM(double pt);

int leadingJetIndex(vector<LorentzVector> jets, int iskip1, int iskip2);
unsigned int getNJets(const float etacut = 2.4);

Float_t getCSVNonb(vector<LorentzVector> jets, vector<float> csv, float discr, double ptTH, double etaTH, bool doTAU);
vector<int> getBJetIndex(double discr, int iskip1, int iskip2, vector<LorentzVector> jets, vector<float> csv, vector<float> lrm, double ptTH, double etaTH, bool doLRM, bool doTAU);

pair<int, int> getIndexPair(vector<int> listBJetIndex, vector<LorentzVector> jets, bool docleaning,int iskip1, int iskip2);


float getdltrigweight_pteta(int id_l1, float l1_pt, float l1_eta, int id_l2, float l2_pt, float l2_eta);
float getleadingtrigweight(int id_leading, float leading_pt, float leading_eta);
float gettrailingtrigweight(int id_trailing, float trailing_pt, float trailing_eta);

float getdltrigweight(int id1, int id2);
float getsltrigweight(int id1, float pt, float eta);
float getideffweight(int id1, float pt, float eta);
float getisoeffweight(int id1, float pt, float eta);

bool passEvtSelection(TString name, bool dometdphi=true);
bool passOneLeptonSelection(bool isData);
bool passTwoLeptonSelection(bool isData);
bool passMVAJetId(double corjetpt, double jeteta, double mvavalue, unsigned int tightness);
bool passTauVeto();
bool passIsoTrkVeto();
bool passIsoTrkVeto_v2();
bool passIsoTrkVeto_v3();
bool passIsoTrkVeto_v4();
bool passSingleLeptonSelection(bool isData);
bool passDileptonSelection(bool isData);
bool passDileptonSSSelection(bool isData);
bool passDileptonSelectionWithEndcapEls(bool isData);
bool passDileptonSSSelectionWithEndcapEls(bool isData);
bool passLepPlusTauSelection(bool isData);
bool passLepPlusTauSelection_v2(bool isData);
bool passLepPlusIsoTrkSelection(bool isData);
bool passLepPlusIsoTrkSelectionWHMet(bool isData);
bool passLepPlusIsoTrkSelection_noEMu(bool isData, bool isMu);
bool passLepPlusSSPionTrkSelection(bool isData);

bool pass_T2tt_LM(bool isData, TString name);
bool pass_T2tt_HM(bool isData, TString name);
bool pass_T2bw_HM(bool isData, TString name);

pair<float,float> getPhiCorrMET( float met, float metphi, int nvtx, bool ismc);

float getDataMCRatio(float eta);
float getDataMCRatioOld(float eta);
float getDataMCRatioSystDown(float eta);
float getDataMCRatioSystUp(float eta);
float getDataMCRatioFix(float eta);

float vtxweight_n( const int nvertices, TH1F *hist, bool isData );
float getdphi( float phi1 , float phi2 );
float dRbetweenVectors(LorentzVector& vec1,LorentzVector& vec2 );
float getMinDphi(float metPhi, LorentzVector& vec1, LorentzVector& vec2 );
float getMT( float pt1 , float phi1 , float pt2 , float phi2 );

struct DorkyEventIdentifier {
    // this is a workaround for not having unique event id's in MC
    unsigned long int run, event,lumi;
    bool operator < (const DorkyEventIdentifier &) const;
    bool operator == (const DorkyEventIdentifier &) const;
};

bool is_duplicate (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &already_seen);
int load_badlaserevents  (char* filename, std::set<DorkyEventIdentifier> &events_lasercalib);
bool is_badLaserEvent (const DorkyEventIdentifier &id, std::set<DorkyEventIdentifier> &events_lasercalib);

enum MT2Type { MT2b, MT2bl, MT2w };
double calculateMT2w(vector<LorentzVector>& jets, vector<float>& btag, LorentzVector& lep, float met, float metphi, MT2Type mt2type = MT2w);
double mt2wWrapper(LorentzVector& lep, LorentzVector& jet_o, LorentzVector& jet_b, float met, float metphi, MT2Type mt2type = MT2w);
double calculateChi2(vector<LorentzVector>& jets, vector<float>& sigma_jets);
double calculateChi2SNT(vector<LorentzVector>& jets, vector<float>& sigma_jets, vector<float>& btag);
double fc2 (double c1, double m12, double m22, double m02, bool verbose);
double fchi2 (double c1, double pt1, double sigma1, double pt2, double sigma2, double m12, double m22, double m02);
void minuitFunction(int&, double* , double &result, double par[], int);
double getChi2(LorentzVector& jets_b, LorentzVector& jets_j1, LorentzVector& jets_j2, float sigma_b, float sigma_j1, float sigma_j2);

int getRegionNumber(float mstep, float mlsp);
double TopPtWeight(double topPt);

float getBtagSF(float pt, float eta, int id, bool fastsim = false);

static const float BTAG_MED = 0.679;
static const float PDG_TOP_MASS = 173.5;
static const float PDG_W_MASS = 80.385;



#endif
