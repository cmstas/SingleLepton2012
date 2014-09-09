Double_t pi = TMath::Pi();

Double_t bins_lepChargeAsym[] =  { -2., -68./60., -48./60., -32./60., -20./60., -8./60., 0., 8./60., 20./60., 32./60., 48./60., 68./60., 2.};
//Double_t bins_lepChargeAsym[] =  { -2., -1.4, -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.4, 2.};
//Double_t bins_lepChargeAsym[] =  { -2., -0.8, -0.4, 0., 0.4, 0.8, 2.};
Double_t bins_lepAzimAsym2[] = {0., 5.*pi/60., 10.*pi/60., 15.*pi/60., 20.*pi/60., 25.*pi/60., 30.*pi/60., 35.*pi/60., 40.*pi/60., 45.*pi/60., 50.*pi/60., 55.*pi/60., pi};
//Double_t bins_lepAzimAsym2[] = {0., 2.*pi/20., 4.*pi/20., 5.5*pi/20., 7.*pi/20., 8.5*pi/20., 10.*pi/20., 11.5*pi/20., 13.*pi/20., 14.5*pi/20., 16.*pi/20., 18*pi/20., pi}; 
//Double_t bins_lepAzimAsym2[] = {0., 4.*pi/20., 7.*pi/20., 10.*pi/20., 13.*pi/20., 16.*pi/20., pi};
Double_t bins_lepAzimAsym[] = {-pi, -50.*pi/60., -40.*pi/60., -30.*pi/60., -20.*pi/60., -10.*pi/60.,  0., 10.*pi/60., 20.*pi/60., 30.*pi/60., 40.*pi/60., 50.*pi/60., pi};
//Double_t bins_lepAzimAsym[] = {-pi, -16.*pi/20., -13.*pi/20., -10.*pi/20., -7.*pi/20., -4.*pi/20.,  0., 4.*pi/20., 7.*pi/20., 10.*pi/20., 13.*pi/20., 16.*pi/20., pi};
//Double_t bins_lepAzimAsym[] = {-1., -0.8, -0.4, 0., 0.4, 0.8, 1.};
Double_t bins_topCosTheta[] = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.}; 
//Double_t bins_topCosTheta[] = {-1., -0.7, -0.4, 0., 0.4, 0.7, 1.}; 
Double_t bins_pseudorapiditydiff[] =  { -2., -1.0, -28./60., 0., 28./60., 1.0, 2.}; 
//Double_t bins_pseudorapiditydiff[] =  { -2., -1.0, -0.5, 0., 0.5, 1.0, 2.}; 
Double_t bins_rapiditydiff[] =  { -2., -1.0, -20./60., 0., 20./60., 1.0, 2.}; 
//Double_t bins_rapiditydiff[] =  { -2., -0.8, -0.3, 0., 0.3, 0.8, 2.}; 
Double_t bins_rapiditydiffMarco[] =  { -2., -44./60., -20./60., 0., 20./60., 44./60., 2.}; 
//Double_t bins_rapiditydiffMarco[] =  { -2., -0.7, -0.3, 0., 0.3, 0.7, 2.}; 
Double_t bins_lepCosTheta[] = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.};
//Double_t bins_lepCosTheta[] = {-1., -0.6, -0.3, 0., 0.3, 0.6, 1.}; 
Double_t bins_topSpinCorr[] = {-1., -0.4, -10./60., 0., 10./60., 0.4, 1.}; 
//Double_t bins_topSpinCorr[] = {-1., -0.5, -0.2, 0., 0.2, 0.5, 1.}; 
Double_t bins_lepCosOpeningAngle[] = {-1., -2./3., -1./3., 0., 1./3., 2./3., 1.}; 
//Double_t bins_lepCosOpeningAngle[] = {-1., -0.6, -0.3, 0., 0.3, 0.6, 1.}; 


Double_t ybinsmtt[] = {0., 430., 530., 1200.};
Double_t ybinsttpt[] = {0., 41., 92., 300.};
Double_t ybinsttrapidity2[] = {0., 0.34, 0.75, 1.5};

//minimum bin width should not be smaller than guide from resolution plots (TopAFBLooper/resolutionplots.sh)
//approximate resolution for lepChargeAsym: ~0
//approximate resolution for lepAzimAsym2: ~0
//approximate resolution for lepAzimAsym: ~0
//approximate resolution for topCosTheta: 0.14
//approximate resolution for pseudorapiditydiff: 0.41 
//approximate resolution for rapiditydiff: 0.35
//approximate resolution for rapiditydiffMarco: 0.34
//approximate resolution for lepCosTheta: 0.17
//approximate resolution for topSpinCorr: 0.15
//approximate resolution for lepCosOpeningAngle: 0.18
