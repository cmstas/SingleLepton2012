void macro_compare_channels_acceptance(){
gROOT->ProcessLine(".L compare_channels_acceptance.C+");
compare_channels_acceptance("lepAzimAsym2");
compare_channels_acceptance("lepAzimAsym");
compare_channels_acceptance("lepChargeAsym");
compare_channels_acceptance("topCosTheta");
compare_channels_acceptance("lepPlusCosTheta");
compare_channels_acceptance("lepMinusCosTheta");
compare_channels_acceptance("topSpinCorr");
compare_channels_acceptance("lepCosOpeningAngle");
compare_channels_acceptance("rapiditydiff");
compare_channels_acceptance("pseudorapiditydiff");
compare_channels_acceptance("rapiditydiffMarco");
compare_channels_acceptance("lepCosTheta");
}