void macro_compare_channels_denominator(){
gROOT->ProcessLine(".L compare_channels_denominator.C+");
compare_channels_denominator("lepAzimAsym2");
compare_channels_denominator("lepAzimAsym");
compare_channels_denominator("lepChargeAsym");
compare_channels_denominator("topCosTheta");
compare_channels_denominator("lepPlusCosTheta");
compare_channels_denominator("lepMinusCosTheta");
compare_channels_denominator("topSpinCorr");
compare_channels_denominator("lepCosOpeningAngle");
compare_channels_denominator("rapiditydiff");
compare_channels_denominator("pseudorapiditydiff");
compare_channels_denominator("rapiditydiffMarco");
compare_channels_denominator("lepCosTheta");
}