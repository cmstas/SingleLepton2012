void macro_comparescaleandtptvariations(){
gROOT->ProcessLine(".L comparescaleandtptvariations.C+");
comparescaleandtptvariations("ttpTGen");
comparescaleandtptvariations("ttMassGen");
comparescaleandtptvariations("lepAzimAsym2Gen");
comparescaleandtptvariations("lepAzimAsymGen");
comparescaleandtptvariations("lepChargeAsymGen");
comparescaleandtptvariations("topCosThetaGen");
comparescaleandtptvariations("lepPlusCosThetaGen");
comparescaleandtptvariations("lepMinusCosThetaGen");
comparescaleandtptvariations("topSpinCorrGen");
comparescaleandtptvariations("lepCosOpeningAngleGen");
comparescaleandtptvariations("rapiditydiffGen");
comparescaleandtptvariations("pseudorapiditydiffGen");
comparescaleandtptvariations("rapiditydiffMarcoGen");
comparescaleandtptvariations("lepCosThetaGen");
}