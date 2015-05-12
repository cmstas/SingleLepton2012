void macro_denominatorPDFvars(){
gROOT->ProcessLine(".L denominatorPDFvars.C+");
denominatorPDFvars("lepAzimAsym2");
denominatorPDFvars("lepAzimAsym");
denominatorPDFvars("lepChargeAsym");
denominatorPDFvars("topCosTheta");
denominatorPDFvars("lepPlusCosTheta");
denominatorPDFvars("lepMinusCosTheta");
denominatorPDFvars("topSpinCorr");
denominatorPDFvars("lepCosOpeningAngle");
denominatorPDFvars("rapiditydiff");
denominatorPDFvars("pseudorapiditydiff");
denominatorPDFvars("rapiditydiffMarco");
denominatorPDFvars("lepCosTheta");
}