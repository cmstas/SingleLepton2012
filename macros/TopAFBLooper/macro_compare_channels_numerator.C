void macro_compare_channels_numerator(){
gROOT->ProcessLine(".L compare_channels_numerator.C+");
compare_channels_numerator("lepAzimAsym2");
compare_channels_numerator("lepAzimAsym");
compare_channels_numerator("lepChargeAsym");
compare_channels_numerator("topCosTheta");
compare_channels_numerator("lepPlusCosTheta");
compare_channels_numerator("lepMinusCosTheta");
compare_channels_numerator("topSpinCorr");
compare_channels_numerator("lepCosOpeningAngle");
compare_channels_numerator("rapiditydiff");
compare_channels_numerator("pseudorapiditydiff");
compare_channels_numerator("rapiditydiffMarco");
compare_channels_numerator("lepCosTheta");
}