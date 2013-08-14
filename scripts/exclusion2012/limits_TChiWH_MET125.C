float getUpperLimit_TChiWH_MET125( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 7.2;
     if(seff >= 0.05 && seff < 0.10) ul = 7.2;
     if(seff >= 0.10 && seff < 0.15) ul = 7.4;
     if(seff >= 0.15 && seff < 0.20) ul = 7.6;
     if(seff >= 0.20 && seff < 0.25) ul = 7.8;
     if(seff >= 0.25 && seff < 0.30) ul = 8.1;
     return ul;
}
float getExpectedUpperLimit_TChiWH_MET125( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 6.4;
     if(seff >= 0.05 && seff < 0.10) ul = 6.4;
     if(seff >= 0.10 && seff < 0.15) ul = 6.6;
     if(seff >= 0.15 && seff < 0.20) ul = 6.7;
     if(seff >= 0.20 && seff < 0.25) ul = 6.9;
     if(seff >= 0.25 && seff < 0.30) ul = 7.1;
     return ul;
}
float getExpectedP1UpperLimit_TChiWH_MET125( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 9.0;
     if(seff >= 0.05 && seff < 0.10) ul = 9.2;
     if(seff >= 0.10 && seff < 0.15) ul = 9.5;
     if(seff >= 0.15 && seff < 0.20) ul = 9.9;
     if(seff >= 0.20 && seff < 0.25) ul = 10.3;
     if(seff >= 0.25 && seff < 0.30) ul = 10.9;
     return ul;
}
float getExpectedM1UpperLimit_TChiWH_MET125( float seff ){
     float ul = 9999.;
     if(seff >= 0.00 && seff < 0.05) ul = 4.7;
     if(seff >= 0.05 && seff < 0.10) ul = 4.7;
     if(seff >= 0.10 && seff < 0.15) ul = 4.8;
     if(seff >= 0.15 && seff < 0.20) ul = 4.8;
     if(seff >= 0.20 && seff < 0.25) ul = 4.9;
     if(seff >= 0.25 && seff < 0.30) ul = 5.0;
     return ul;
}
