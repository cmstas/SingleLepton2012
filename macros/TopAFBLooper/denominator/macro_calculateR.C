void macro_calculateR(){
gROOT->ProcessLine(".L calculateR.C+");
calculateR("lepAzimAsym2");
calculateR("lepAzimAsym");
calculateR("lepChargeAsym");
calculateR("topCosTheta");
calculateR("lepPlusCosTheta");
calculateR("lepMinusCosTheta");
calculateR("topSpinCorr");
calculateR("lepCosOpeningAngle");
calculateR("rapiditydiff");
calculateR("pseudorapiditydiff");
calculateR("rapiditydiffMarco");
calculateR("lepCosTheta");
}