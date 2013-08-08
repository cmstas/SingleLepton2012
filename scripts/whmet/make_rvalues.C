#include <vector>
#include "TFile.h"
#include "TH2F.h"

void fill_hists(int xbin, int ybin, float* nsig, std::vector<TH2F*>& hists, int force_sr = -1) {

  // results from LandS:
  // (obs, -2s,-1s,median,1s,2s)

  // ------ bg syst == 40%
  // - MET > 100 SR:
  // 7.8119, 4.7331, 5.64857, 7.46563, 10.7313, 14.8455
  // - MET > 150 SR:
  // 5.19088, 3.58843, 3.94099, 5.28198, 7.94436, 12.708
  // - MET > 175 SR:
  // 4.61047, 1.52635, 3.4525, 4.45819, 6.21637, 9.60939

  // float limits[3][6] = { {7.8119, 4.7331, 5.64857, 7.46563, 10.7313, 14.8455},
  // 			 {5.19088, 3.58843, 3.94099, 5.28198, 7.94436, 12.708},
  // 			 {4.61047, 1.52635, 3.4525, 4.45819, 6.21637, 9.60939} };

  // ------ bg syst == 50%
  // - MET > 100 SR:
  // 8.2622, 2.75011, 6.06945, 7.94679, 11.1489, 17.0695
  // - MET > 150 SR:
  // 5.32859, 3.67591, 4.07899, 5.39979, 8.08415, 13.7163
  // - MET > 175 SR:
  // 4.6892, 2.70454, 3.50945, 4.52161, 6.12764, 9.08404

  // float limits[3][6] = { {8.2622, 2.75011, 6.06945, 7.94679, 11.1489, 17.0695},
  // 			 {5.32859, 3.67591, 4.07899, 5.39979, 8.08415, 13.7163},
  // 			 {4.6892, 2.70454, 3.50945, 4.52161, 6.12764, 9.08404} };

  // ------ yields Jul4, realistic bg systs, 13.5% signal syst
  // - MET > 100 SR:
  // 8.32223, 4.45433, 5.60809, 7.73236, 11.044, 16.4581
  // - MET > 150 SR:
  // 6.1081, 3.24599, 4.0928, 5.56034, 8.1045, 11.3873
  // - MET > 175 SR:
  // 4.64393, 2.83425, 3.29131, 4.44521, 6.44999, 9.21175

  // float limits[3][6] = { {8.32223, 4.45433, 5.60809, 7.73236, 11.044, 16.4581},
  // 			 {6.1081, 3.24599, 4.0928, 5.56034, 8.1045, 11.3873},
  // 			 {4.64393, 2.83425, 3.29131, 4.44521, 6.44999, 9.21175} };

  // ------ yields Jul8, realistic bg systs, 13.5% signal syst
  // - MET > 100 SR:
  // 7.50019, 4.10332, 5.33914, 7.39089, 10.7542, 15.0601
  // - MET > 125 SR:
  // 6.51427, 3.59698, 4.49405, 6.15729, 8.98277, 12.6002
  // - MET > 150 SR:
  // 6.21536, 3.24626, 4.17659, 5.56758, 8.07044, 11.3485
  // - MET > 175 SR:
  // 4.65729, 2.9127, 3.34101, 4.3766, 6.33029, 9.04181

  // float limits[4][6] = { {7.50019, 4.10332, 5.33914, 7.39089, 10.7542, 15.0601},
  // 			 {6.51427, 3.59698, 4.49405, 6.15729, 8.98277, 12.6002},
  // 			 {6.21536, 3.24626, 4.17659, 5.56758, 8.07044, 11.3485},
  // 			 {4.65729, 2.9127, 3.34101, 4.3766, 6.33029, 9.04181} };

  // ------ yields Jul9, realistic bg systs, unblinded, 13.5% signal syst
  // - MET > 100 SR:
  // 7.76184, 4.29746, 5.45724, 7.60852, 11.0095, 15.5778
  // - MET > 125 SR:
  // 7.94583, 3.85777, 4.77677, 6.53362, 9.36938, 13.3155
  // - MET > 150 SR:
  // 5.15232, 3.22892, 4.03365, 5.48868, 8.00462, 12.327
  // - MET > 175 SR:
  // 5.97494, 2.99493, 3.45368, 4.60862, 6.64367, 9.43689

  // float limits[4][6] = { {7.76184, 4.29746, 5.45724, 7.60852, 11.0095, 15.5778},
  // 			 {7.94583, 3.85777, 4.77677, 6.53362, 9.36938, 13.3155},
  // 			 {5.15232, 3.22892, 4.03365, 5.48868, 8.00462, 12.327},
  // 			 {5.97494, 2.99493, 3.45368, 4.60862, 6.64367, 9.43689} };

  // ------ yields Jul9, realistic bg systs, unblinded, 13.5% signal syst
  //   --- expected limits from Asimov data
  //   --- 5x as many toys as previous limits
  // - MET > 100 SR:
  // 7.66363, 4.39113, 5.54092, 7.62098, 11.0678, 16.0925
  // - MET > 125 SR:
  // 7.89839, 3.76562, 4.6066, 6.39732, 9.26266, 13.3343
  // - MET > 150 SR:
  // 5.18007, 3.33889, 4.15183, 5.60542, 8.20232, 12.0207
  // - MET > 175 SR:
  // 6.0058, 2.79165, 3.30581, 4.45221, 6.48366, 9.43796

  // float limits[4][6] = { {7.66363, 4.39113, 5.54092, 7.62098, 11.0678, 16.0925},
  // 			 {7.89839, 3.76562, 4.6066, 6.39732, 9.26266, 13.3343},
  // 			 {5.18007, 3.33889, 4.15183, 5.60542, 8.20232, 12.0207},
  // 			 {6.0058, 2.79165, 3.30581, 4.45221, 6.48366, 9.43796} };

  // ------ yields Jul31, adding WZ->lvbb background, realistic signal systs (10\% stat)
  //   --- expected limits from Asimov data, obs from normal (smaller toys)
  //   --- 5x as many toys as previous limits
  // - MET > 100 SR:
  // 7.44818, 4.72358, 5.76072, 7.92515, 11.3552, 18.8855
  // - MET > 125 SR:
  // 7.44357, 3.93421, 4.8732, 6.70977, 9.85079, 14.9105
  // - MET > 150 SR:
  // 4.77414, 3.82953, 4.64403, 6.37997, 9.16173, 14.5715
  // - MET > 175 SR:
  // 5.79868, 2.94812, 3.50832, 4.72142, 6.82149, 9.84271

  // float limits[4][6] = { {7.44818, 4.72358, 5.76072, 7.92515, 11.3552, 18.8855},
  // 			 {7.44357, 3.93421, 4.8732, 6.70977, 9.85079, 14.9105},
  // 			 {4.77414, 3.82953, 4.64403, 6.37997, 9.16173, 14.5715},
  // 			 {5.79868, 2.94812, 3.50832, 4.72142, 6.82149, 9.84271} };

  // ------ yields Aug7, using Wbb NLO xsec, realistic signal systs (10\% stat)
  //   --- expected limits from Asimov data, obs from normal (smaller toys)
  //   --- 5x as many toys as previous limits
  // - MET > 100 SR:
  // 7.082, 4.5777, 5.70832, 7.79504, 11.3513, 16.0816
  // - MET > 125 SR:
  // 7.42825, 3.79622, 4.67491, 6.47775, 9.47928, 13.7308
  // - MET > 150 SR:
  // 5.07166, 3.39727, 4.14893, 5.63084, 8.20909, 11.8783
  // - MET > 175 SR:
  // 5.71997, 2.96702, 3.52906, 4.76454, 6.91856, 10.0066

  float limits[4][6] = { {7.082, 4.5777, 5.70832, 7.79504, 11.3513, 16.0816},
  			 {7.42825, 3.79622, 4.67491, 6.47775, 9.47928, 13.7308},
  			 {5.07166, 3.39727, 4.14893, 5.63084, 8.20909, 11.8783},
  			 {5.71997, 2.96702, 3.52906, 4.76454, 6.91856, 10.0066} };

  int best_sr = -1;
  float best_rval = 100.;
  for (int i = 0; i < 4; ++i) {
    float rval = limits[i][3]/nsig[i];
    std::cout << rval << std::endl;
    if (rval < best_rval) {
      best_rval = rval;
      best_sr = i;
    }
  }

  std::cout << "selected sr: " << best_sr << ", nsig: " << nsig[best_sr] << std::endl;

  if (force_sr > -1) {
    std::cout << "forcing sr: " << force_sr << ", nsig: " << nsig[force_sr] << std::endl;
    best_sr = force_sr;
  }

  for (int j = 0; j < 6; ++j) {
    hists[j]->SetBinContent(xbin,ybin,limits[best_sr][j]/nsig[best_sr]);
  }

  return;
}

void make_rvalues() {

  int force_sr = -1;

  TFile* f = new TFile("r-values.root","RECREATE");

  const int nbinsx = 8;
  const double xbins[nbinsx+1] = {120,140,160,190,210,290,310,390,410};
  //  const int nbinsx = 6;
  // const float xmin = 125.;
  // const float xmax = 375.;
  const int nbinsy = 1;
  const double ymin = 0.;
  const double ymax = 50.;

  TH2F* hObs = new TH2F("hObs","", nbinsx, xbins, nbinsy, ymin, ymax);
  TH2F* hExp2m = new TH2F("hExp2m","", nbinsx, xbins, nbinsy, ymin, ymax);
  TH2F* hExp1m = new TH2F("hExp1m","", nbinsx, xbins, nbinsy, ymin, ymax);
  TH2F* hExp = new TH2F("hExp","", nbinsx, xbins, nbinsy, ymin, ymax);
  TH2F* hExp1p = new TH2F("hExp1p","", nbinsx, xbins, nbinsy, ymin, ymax);
  TH2F* hExp2p = new TH2F("hExp2p","", nbinsx, xbins, nbinsy, ymin, ymax);

  std::vector<TH2F*> hists;
  hists.push_back(hObs);
  hists.push_back(hExp2m);
  hists.push_back(hExp1m);
  hists.push_back(hExp);
  hists.push_back(hExp1p);
  hists.push_back(hExp2p);

  // mass 130
  std::cout << "-- mass 130,1" << std::endl;
  float vals_130[] = {7.2,6.3,5.1,3.8}; 
  fill_hists(1,1,vals_130,hists);

  // mass 150
  std::cout << "-- mass 150,1" << std::endl;
  float vals_150[] = {7.0,5.9,4.6,3.3}; 
  fill_hists(2,1,vals_150,hists,force_sr);

  // mass 175
  std::cout << "-- mass 175,1" << std::endl;
  float vals_175[] = {7.2,5.5,4.4,3.5};
  fill_hists(3,1,vals_175,hists,force_sr);

  // mass 200
  std::cout << "-- mass 200,1" << std::endl;
  float vals_200[] = {8.2,7.0,5.3,3.9};
  fill_hists(4,1,vals_200,hists);

  // mass 250
  std::cout << "-- mass 250,1" << std::endl;
  float vals_250[] = {7.5,6.7,5.8,4.9};
  fill_hists(5,1,vals_250,hists);

  // mass 300
  std::cout << "-- mass 300,1" << std::endl;
  float vals_300[] = {6.7,6.3,5.5,4.6};
  fill_hists(6,1,vals_300,hists);

  // mass 350
  std::cout << "-- mass 350,1" << std::endl;
  float vals_350[] = {4.6,4.3,3.9,3.4};
  fill_hists(7,1,vals_350,hists);

  // mass 400
  std::cout << "-- mass 400,1" << std::endl;
  float vals_400[] = {3.3,3.2,3.1,2.7};
  fill_hists(8,1,vals_400,hists);

  f->Write();
  f->Close();

}
