float getUpperLimit_TChiWH_MET150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.9;
     if(seff >= 0.05 && seff < 0.10) ul = 5.0;
     if(seff >= 0.10 && seff < 0.15) ul = 5.0;
     if(seff >= 0.15 && seff < 0.20) ul = 5.1;
     if(seff >= 0.20 && seff < 0.25) ul = 5.3;
     if(seff >= 0.25 && seff < 0.30) ul = 5.3;
     return ul;
}
float getExpectedUpperLimit_TChiWH_MET150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 5.5;
     if(seff >= 0.05 && seff < 0.10) ul = 5.5;
     if(seff >= 0.10 && seff < 0.15) ul = 5.6;
     if(seff >= 0.15 && seff < 0.20) ul = 5.8;
     if(seff >= 0.20 && seff < 0.25) ul = 5.9;
     if(seff >= 0.25 && seff < 0.30) ul = 6.0;
     return ul;
}
float getExpectedP1UpperLimit_TChiWH_MET150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.7;
     if(seff >= 0.05 && seff < 0.10) ul = 7.9;
     if(seff >= 0.10 && seff < 0.15) ul = 8.1;
     if(seff >= 0.15 && seff < 0.20) ul = 8.4;
     if(seff >= 0.20 && seff < 0.25) ul = 8.8;
     if(seff >= 0.25 && seff < 0.30) ul = 9.3;
     return ul;
}
float getExpectedM1UpperLimit_TChiWH_MET150( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.0;
     if(seff >= 0.05 && seff < 0.10) ul = 4.1;
     if(seff >= 0.10 && seff < 0.15) ul = 4.1;
     if(seff >= 0.15 && seff < 0.20) ul = 4.2;
     if(seff >= 0.20 && seff < 0.25) ul = 4.2;
     if(seff >= 0.25 && seff < 0.30) ul = 4.3;
     return ul;
}
