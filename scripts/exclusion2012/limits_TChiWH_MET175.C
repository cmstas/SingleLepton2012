float getUpperLimit_TChiWH_MET175( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.5;
     if(seff >= 0.05 && seff < 0.10) ul = 5.6;
     if(seff >= 0.10 && seff < 0.15) ul = 5.6;
     if(seff >= 0.15 && seff < 0.20) ul = 5.8;
     if(seff >= 0.20 && seff < 0.25) ul = 6.0;
     if(seff >= 0.25 && seff < 0.30) ul = 6.2;
     return ul;
}
float getExpectedUpperLimit_TChiWH_MET175( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.5;
     if(seff >= 0.05 && seff < 0.10) ul = 4.5;
     if(seff >= 0.10 && seff < 0.15) ul = 4.6;
     if(seff >= 0.15 && seff < 0.20) ul = 4.7;
     if(seff >= 0.20 && seff < 0.25) ul = 4.8;
     if(seff >= 0.25 && seff < 0.30) ul = 5.0;
     return ul;
}
float getExpectedP1UpperLimit_TChiWH_MET175( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 6.4;
     if(seff >= 0.05 && seff < 0.10) ul = 6.5;
     if(seff >= 0.10 && seff < 0.15) ul = 6.7;
     if(seff >= 0.15 && seff < 0.20) ul = 7.0;
     if(seff >= 0.20 && seff < 0.25) ul = 7.2;
     if(seff >= 0.25 && seff < 0.30) ul = 7.6;
     return ul;
}
float getExpectedM1UpperLimit_TChiWH_MET175( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 3.4;
     if(seff >= 0.05 && seff < 0.10) ul = 3.4;
     if(seff >= 0.10 && seff < 0.15) ul = 3.5;
     if(seff >= 0.15 && seff < 0.20) ul = 3.5;
     if(seff >= 0.20 && seff < 0.25) ul = 3.6;
     if(seff >= 0.25 && seff < 0.30) ul = 3.6;
     return ul;
}
