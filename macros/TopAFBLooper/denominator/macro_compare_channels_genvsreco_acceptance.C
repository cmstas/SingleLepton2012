void macro_compare_channels_genvsreco_acceptance(){
gROOT->ProcessLine(".L compare_channels_genvsreco_acceptance.C+");
compare_channels_genvsreco_acceptance("lepAzimAsym2");
compare_channels_genvsreco_acceptance("lepAzimAsym");
compare_channels_genvsreco_acceptance("lepChargeAsym");
compare_channels_genvsreco_acceptance("topCosTheta");
compare_channels_genvsreco_acceptance("lepPlusCosTheta");
compare_channels_genvsreco_acceptance("lepMinusCosTheta");
compare_channels_genvsreco_acceptance("topSpinCorr");
compare_channels_genvsreco_acceptance("lepCosOpeningAngle");
compare_channels_genvsreco_acceptance("rapiditydiff");
compare_channels_genvsreco_acceptance("pseudorapiditydiff");
compare_channels_genvsreco_acceptance("rapiditydiffMarco");
compare_channels_genvsreco_acceptance("lepCosTheta");
}