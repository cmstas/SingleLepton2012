float getUpperLimit_TChiWH_MET100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.2;
     if(seff >= 0.05 && seff < 0.10) ul = 7.2;
     if(seff >= 0.10 && seff < 0.15) ul = 7.3;
     if(seff >= 0.15 && seff < 0.20) ul = 7.5;
     if(seff >= 0.20 && seff < 0.25) ul = 7.6;
     if(seff >= 0.25 && seff < 0.30) ul = 7.9;
     return ul;
}
float getExpectedUpperLimit_TChiWH_MET100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.6;
     if(seff >= 0.05 && seff < 0.10) ul = 7.6;
     if(seff >= 0.10 && seff < 0.15) ul = 7.8;
     if(seff >= 0.15 && seff < 0.20) ul = 8.0;
     if(seff >= 0.20 && seff < 0.25) ul = 8.2;
     if(seff >= 0.25 && seff < 0.30) ul = 8.5;
     return ul;
}
float getExpectedP1UpperLimit_TChiWH_MET100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 10.6;
     if(seff >= 0.05 && seff < 0.10) ul = 10.9;
     if(seff >= 0.10 && seff < 0.15) ul = 11.2;
     if(seff >= 0.15 && seff < 0.20) ul = 11.7;
     if(seff >= 0.20 && seff < 0.25) ul = 12.2;
     if(seff >= 0.25 && seff < 0.30) ul = 12.9;
     return ul;
}
float getExpectedM1UpperLimit_TChiWH_MET100( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.6;
     if(seff >= 0.05 && seff < 0.10) ul = 5.6;
     if(seff >= 0.10 && seff < 0.15) ul = 5.7;
     if(seff >= 0.15 && seff < 0.20) ul = 5.8;
     if(seff >= 0.20 && seff < 0.25) ul = 5.9;
     if(seff >= 0.25 && seff < 0.30) ul = 6.0;
     return ul;
}
