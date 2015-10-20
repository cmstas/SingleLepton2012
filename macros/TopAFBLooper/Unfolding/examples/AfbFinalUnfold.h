//C++ includes
#include <algorithm>
#include <fstream>

//ROOT includes
#include "TH1D.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TMath.h"
#include "TUnfoldSys.h"
#include "TMinuit.h"

//Analysis-specific includes
#include "../../denominator/acceptanceplots.h" //import binning from acceptanceplots

double xsection = 154.0*0.06451;

int lineWidth=3;

TString observablename;
TString acceptanceName;
TString xaxislabel;
TString yaxislabel;
TString asymlabel;

float xmin1D=-1.0;
float xmax1D= 1.0;
float xmin=-1.0;
float xmax= 1.0;
float ymin=-1.0;
float ymax= 1.0;

const int nbins1D=6;
const int nbinsx2D=6; //Normal configuration: use 6 bins in x, before and after unfolding.
const int nbinsx2Dalt = 12; //Alternate configuration: use 12 bins in x before unfolding, and 6 bins after.
const int nbinsy2D=3;
//const int nbinsunwrapped = nbinsx2D*nbinsy2D;
//const int nbinsunwrappedalt = nbinsx2Dalt*nbinsy2D;

double xbins1D[nbins1D+1];
double xbins2D[nbinsx2D+1];
double xbins2Dalt[nbinsx2Dalt+1];
double ybins2D[nbinsy2D+1];


Double_t stat_corr  [nbinsx2Dalt]; //errors include syst error in the unfolding
Double_t stat_uncorr[nbinsx2Dalt]; //errors do not include syst error in the unfolding
Double_t syst_corr  [nbinsx2Dalt];
Double_t hardcodeddata  [nbinsx2Dalt];
Double_t theory_default[nbinsx2Dalt];
Double_t theory_scaleup[nbinsx2Dalt];
Double_t theory_scaledown[nbinsx2Dalt];
Double_t uncorr_default[nbinsx2Dalt];
Double_t uncorr_scaleup[nbinsx2Dalt];
Double_t uncorr_scaledown[nbinsx2Dalt];

Float_t sign(Float_t t) 
{
  if( t >= 0.0 )
    return 1.0;
  else
    return -1.0;
}


void GetAfb(TH1D* h, Float_t &afb, Float_t  &afberr){

  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  //event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"");
  event_total = event_plus + event_minus;

  //std::cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<std::endl;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

}

void GetAfb(TH2D* h, Float_t &afb, Float_t  &afberr){

  Int_t nbinsx = h->GetNbinsX();
  Int_t nbinsy = h->GetNbinsY();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  event_minus  = h->IntegralAndError(0, nbinsx/2, 0, nbinsy+1, event_minus_err,"");
  event_plus   = h->IntegralAndError(nbinsx/2+1, nbinsx+1, 0, nbinsy+1, event_plus_err,"");
  event_total = event_plus + event_minus;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

}


void GetAfb_integratewidth(TH1D* h, Float_t &afb, Float_t  &afberr){

  Int_t nbins = h->GetNbinsX();
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"width");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_plus_err,"width");
  event_total = event_plus + event_minus;

  afb = (event_plus-event_minus)/(event_plus+event_minus);
  afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

}


//void GetAfbBinByBin(TH1D* h, Float_t &afbbin, Float_t  &afberrbin){
void GetAfbBinByBin(TH1D* h){

  Int_t nbins = h->GetNbinsX();
  const int nbins2 = nbins/2 +1;
  Float_t event_minus[nbins2];
  Float_t event_plus[nbins2];
  Float_t event_total[nbins2];
  Double_t event_plus_err[nbins2];
  Double_t event_minus_err[nbins2];
  Double_t afbbin[nbins2];
  Double_t afberrbin[nbins2];

  Double_t event_plus_total = 0.;
  Double_t event_minus_total = 0.;
  Double_t event_total_total = 0.;

  for(int i=0;i<nbins2;i++){
    //event_minus[i]  = h-> IntegralAndError(i, i, event_plus_err[i],"");
    event_minus[i]  = h-> IntegralAndError(i, i, event_minus_err[i],"");
    event_minus_total += event_minus[i];
    //event_plus[i]   = h-> IntegralAndError(nbins+1-i, nbins+1-i, event_minus_err[i],"");
    event_plus[i]   = h-> IntegralAndError(nbins+1-i, nbins+1-i, event_plus_err[i],"");
    event_plus_total += event_plus[i];
    event_total[i] = event_plus[i] + event_minus[i];
    event_total_total += event_total[i];

    //std::cout<<event_minus[i]<<" "<<event_minus_err[i]<<" "<<event_plus[i]<<" "<<event_plus_err[i]<<" "<<event_total[i]<<std::endl;

    afbbin[i] = (event_plus[i]-event_minus[i])/(event_plus[i]+event_minus[i]);
    afberrbin[i]   = sqrt(4*(event_plus[i]*event_plus[i]*event_minus_err[i]*event_minus_err[i] 
      + event_minus[i]*event_minus[i]*event_plus_err[i]*event_plus_err[i])/
      (event_total[i]*event_total[i]*event_total[i]*event_total[i]));
    std::cout<<i<<" AFB = "<<afbbin[i]<<" +/- "<<afberrbin[i]<<std::endl;
  }
}


void GetAvsY(TH1* histogram, TMatrixD &covarianceM, std::vector<double> &myafb, std::vector<double> &myerr, ofstream& second_output_file){

  myafb.clear();
  myerr.clear();

    //Get info from histogram
  int nbins = histogram->GetNbinsX();
  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Output
  double afb[8], err[8];

    //Setup Some Needed Vectors
  double alpha[16], beta[16], dfdn[16];

    //Get Asymmetry for each Y bin
  for(int j=0;j<nbins/2;j++) {

    int forBin = nbins/2 + j;
    int bacBin = nbins/2 - j - 1;

    for(int i=0;i<nbins;i++){
      if( i == forBin ){ 
        alpha[i] = 1; 
        beta[i] = 1;
      }else if( i == bacBin ){ 
        alpha[i] = -1; 
        beta[i] = 1;
      }else{
        alpha[i] = 0;
        beta[i] = 0;
      }
    }

    double sum = 0. , diff = 0.;
    for(int i=0;i<nbins;i++){
      sum += beta[i] * n[i];
      diff += alpha[i] * n[i];
    }

      //Calculate Everything
    if(sum > 0){ 

  //Error Calculation
      for(int i=0;i<nbins;i++){
        dfdn[i] = ( alpha[i] * sum - beta[i] * diff ) / pow(sum,2);
      }

      double afberr = 0.;
      for(int i=0;i<nbins;i++){
        for(int k=0;k<nbins;k++){
          afberr += covarianceM(i,k) * dfdn[i] * dfdn[k];
      //if(i==k) std::cout<<"DAH: "<<n[i]<<" "<<k<<" "<<covarianceM(i,k)<<std::endl;
        }
      }
      afberr = sqrt(afberr);

      err[j] = afberr;
      afb[j] = diff / sum; 

    }else{ 

      afb[j] = 0.; 
      err[j] = 0.;

    }
    myafb.push_back(afb[j]);
    myerr.push_back(err[j]);

    std::cout<<j<<" AFB = "<<afb[j]<<" +/- "<<err[j]<<std::endl;
    second_output_file << acceptanceName << " " << observablename << " AFB" << j << ": " << afb[j] << " +/- " << err[j] << std::endl;    	
  }
}

void GetAvsY2d (TH2* h2, std::vector<double> &myafb, std::vector<double> &myerr, ofstream& second_output_file){

  myafb.clear();
  myerr.clear();

  Float_t rowafb;
  Float_t rowerr;

  int ny = h2->GetNbinsY();

  TH1D* rowhist = new TH1D();

  for(int i=1; i<=ny; i++) {

	rowhist = h2->ProjectionX("rowprojection", i, i, "");  // I verified that this works correctly even without option "e".

	GetAfb(rowhist, rowafb, rowerr);

	myafb.push_back(rowafb);
	myerr.push_back(rowerr);

	// std::cout << i << " AFB = " << rowafb << " +/- " << rowerr << std::endl;
	// second_output_file << acceptanceName << " " << observablename << " AFB" << i << ": " << rowafb << " +/- " << rowerr << std::endl;
  }

  delete rowhist;
}

void GetCorrectedAfb(TH1D* histogram, TMatrixD &covarianceM, Float_t &afb, Float_t  &afberr){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    //Get histogram info
  int nbins = histogram->GetNbinsX();
  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Setup Alpha Vector
  double alpha[16];
  for(int i=0;i<nbins;i++) {
	if(i < nbins/2 ){ alpha[i] = -1;
	}
	else{ alpha[i] = 1;}
  }

    //Components of the error calculation
  double sum_n = 0.;
  double sum_alpha_n = 0.;
  for(int i=0;i<nbins;i++){
    sum_n += n[i];
    sum_alpha_n += alpha[i] * n[i];
  }

  double dfdn[16];
  for(int i=0;i<nbins;i++){
    dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2);
  }

    //Error Calculation
  afberr = 0.;
  for(int i=0;i<nbins;i++){
    for(int j=0;j<nbins;j++){
      afberr += covarianceM(i,j) * dfdn[i] * dfdn[j];
    }
  }
  afberr = sqrt(afberr);

    //Calculate Afb
  afb = sum_alpha_n / sum_n;

    //    std::cout<<"AFB = "<<afb<<" "<<afberr<<std::endl;
}



void GetCorrectedAfb2d(TH2D* histogram, TMatrixD &covarianceM, std::vector<double> &myafb, std::vector<double> &myerr, TMatrixD &AFBcovarianceM, ofstream& second_output_file){

    //calculate AFB, using covariance matrix, m_correctE(j,i), for the uncertainty

  myafb.clear();
  myerr.clear();

    //Get histogram info
  Int_t nbinsx = histogram->GetNbinsX();
  Int_t nbinsy = histogram->GetNbinsY();

  const Int_t numbinsx = nbinsx;
  const Int_t numbinsy = nbinsy;
  const Int_t numbinsxy = numbinsx*numbinsy;
  const Int_t nbins2 = numbinsx/2;

  double afb[numbinsy+1];
  double afberr[numbinsy+1];
  //double afbcov[numbinsy][numbinsy];

  double afbdoublediff[nbins2][numbinsy];
  double afbdoubledifferr[nbins2][numbinsy];

  memset( afb, 0, sizeof(afb) );  //Initialize these arrays to zero
  memset( afberr, 0, sizeof(afberr) );
  //memset( afbcov, 0, sizeof(afbcov) );

  memset( afbdoublediff, 0, sizeof(afbdoublediff) );  //Initialize these arrays to zero
  memset( afbdoubledifferr, 0, sizeof(afbdoubledifferr) );

  double n[numbinsx][numbinsy];
  for(int i=0;i<numbinsx;i++){
    for(int j=0;j<numbinsy;j++){
      n[i][j] = histogram->GetBinContent(i+1,j+1);
    }
  }

    //Setup Alpha Vector
  double alpha[numbinsx];
  for(int i=0;i<numbinsx;i++) if(i < numbinsx/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
  double sum_n_doublediff[nbins2][numbinsy]; //= {0.};
  double sum_alpha_n_doublediff[nbins2][numbinsy]; //= {0.};

  double sum_n[numbinsy]; //= {0.};
  double sum_alpha_n[numbinsy]; //= {0.};
  double sum_n_Inclusive = 0.;
  double sum_alpha_n_Inclusive = 0.;

  memset( sum_n, 0, sizeof(sum_n) );
  memset( sum_alpha_n, 0, sizeof(sum_alpha_n) );
  memset( sum_n_doublediff, 0, sizeof(sum_n_doublediff) );
  memset( sum_alpha_n_doublediff, 0, sizeof(sum_alpha_n_doublediff) );

  for(int i=0;i<numbinsx;i++){
    for(int j=0;j<numbinsy;j++){
      sum_n[j] += n[i][j];
      sum_alpha_n[j] += alpha[i] * n[i][j];
      sum_n_Inclusive += n[i][j];
      sum_alpha_n_Inclusive += alpha[i] * n[i][j];
      if(i<nbins2) {
        sum_n_doublediff[i][j] = n[i][j] + n[numbinsx-1-i][j] ;
        sum_alpha_n_doublediff[i][j] = alpha[i] * n[i][j] + alpha[numbinsx-1-i] * n[numbinsx-1-i][j] ;        
      }
    }
  }
    //AFB = SUM( n[i][j]*alpha[i] ) / SUM( n[i][j] )
    //dAFB/dn_ij = [ SUM( n[i][j] ) * alpha[i] -  SUM( n[i][j]*alpha[i] ) * ( 1 ) ]  / SUM( n[i][j] )^2
  double dfdn[numbinsx][numbinsy];
  double dfdnInclusive[numbinsx][numbinsy];
  double dfdn_doublediff[numbinsx][numbinsy];

  memset( dfdn, 0, sizeof(dfdn) );
  memset( dfdnInclusive, 0, sizeof(dfdnInclusive) );
  memset( dfdn_doublediff, 0, sizeof(dfdn_doublediff) );


  for(int i=0;i<numbinsx;i++){
    for(int j=0;j<numbinsy;j++){
      dfdn[i][j] = ( alpha[i] * sum_n[j] - sum_alpha_n[j] ) / pow(sum_n[j],2);
      dfdnInclusive[i][j] = ( alpha[i] * sum_n_Inclusive - sum_alpha_n_Inclusive ) / pow(sum_n_Inclusive,2);
    }
  }

  for(int j=0;j<numbinsy;j++){
    for(int i=0;i<numbinsx;i++){
      int k = -999;
      if (i < nbins2) k = i;
      else k = numbinsx-1-i;
      dfdn_doublediff[i][j] = ( alpha[i] * sum_n_doublediff[k][j] - sum_alpha_n_doublediff[k][j] ) / pow(sum_n_doublediff[k][j],2);
    }
  }


  for(int i=0;i<numbinsy;i++){
    for(int j=0;j<numbinsy;j++){
      AFBcovarianceM(i,j) = 0;
    }
  }

    //Error Calculation
  for(int i=0;i<numbinsxy;i++){
    for(int j=0;j<numbinsxy;j++){
      int i_2di = i % numbinsx;
      int i_2dj = i / numbinsx;
      int j_2di = j % numbinsx;
      int j_2dj = j / numbinsx;
      //std::cout<<"i: "<<i<<" "<<i_2di<<" "<<i_2dj<<" j: "<<j<<" "<<j_2di<<" "<<j_2dj<<std::endl;
      afberr[0] += covarianceM(i,j) * dfdnInclusive[i_2di][i_2dj] * dfdnInclusive[j_2di][j_2dj];
      if(i_2dj==j_2dj) afberr[i_2dj+1] += covarianceM(i,j) * dfdn[i_2di][i_2dj] * dfdn[j_2di][j_2dj]; 
      //if(i_2dj==j_2dj) cout<<"bin: "<<i_2dj<<": "<<covarianceM(i,j) <<" "<< dfdn[i_2di][i_2dj] <<" "<< dfdn[j_2di][j_2dj]<<" "<<covarianceM(i,j) * dfdn[i_2di][i_2dj] * dfdn[j_2di][j_2dj]<<" "<<afberr[i_2dj+1]<<endl;
      AFBcovarianceM(i_2dj,j_2dj) += covarianceM(i,j) * dfdn[i_2di][i_2dj] * dfdn[j_2di][j_2dj];
    }
  }

  for(int k=0;k<nbins2;k++){
    for(int i=0;i<numbinsxy;i++){
      for(int j=0;j<numbinsxy;j++){
        int i_2di = i % numbinsx;
        int i_2dj = i / numbinsx;
        int j_2di = j % numbinsx;
        int j_2dj = j / numbinsx;

        if( (i_2dj==j_2dj) &&  (i_2di==k || i_2di==numbinsx-1-k ) && (j_2di==k || j_2di==numbinsx-1-k ) ) {
          afbdoubledifferr[k][i_2dj] += covarianceM(i,j) * dfdn_doublediff[i_2di][i_2dj] * dfdn_doublediff[j_2di][j_2dj];
          //cout<<"DDbin: "<<i_2dj<<"x"<<k<<": "<<covarianceM(i,j) <<" "<< dfdn_doublediff[i_2di][i_2dj] <<" "<< dfdn_doublediff[j_2di][j_2dj]<<" "<<covarianceM(i,j) * dfdn_doublediff[i_2di][i_2dj] * dfdn_doublediff[j_2di][j_2dj]<<" "<<afbdoubledifferr[k][i_2dj]<<endl;
        }
      }
    }
  }

    //Calculate Afb

  afb[0] = sum_alpha_n_Inclusive / sum_n_Inclusive;

  for (int j = 0; j < numbinsy+1; ++j)
  {
    afberr[j] = sqrt(afberr[j]);
    if(j>0) afb[j] = sum_alpha_n[j-1] / sum_n[j-1];
    myafb.push_back(afb[j]);
    myerr.push_back(afberr[j]);

    // std::cout<<j<<" AFB = "<<afb[j]<<" +/- "<<afberr[j]<<std::endl;
    // second_output_file << acceptanceName << " " << observablename << " AFB" << j << ": " << afb[j] << " +/- " << afberr[j] << std::endl; 

  }

  for(int j=0;j<numbinsy;j++){
    for(int k=0;k<nbins2;k++){
      afbdoubledifferr[k][j] = sqrt(afbdoubledifferr[k][j]);
      afbdoublediff[k][j] = sum_alpha_n_doublediff[k][j] / sum_n_doublediff[k][j];
      //std::cout<<k<<" AFB = "<<afbdoubledifferr[k][j]<<" +/- "<<afbdoublediff[k][j]<<std::endl;
      myerr.push_back(afbdoubledifferr[k][j]);
      myafb.push_back(afbdoublediff[k][j]);
    }
  }

}





void GetCorrectedAfb_integratewidth(TH1D* histogram, TMatrixD &covarianceM, Float_t &afb, Float_t  &afberr){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    //Get histogram info
  int nbins = histogram->GetNbinsX();
  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1) * histogram->GetBinWidth(i+1);
  }

    //Setup Alpha Vector
  double alpha[16];
  for(int i=0;i<nbins;i++) if(i < nbins/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
  double sum_n = 0.;
  double sum_alpha_n = 0.;
  for(int i=0;i<nbins;i++){
    sum_n += n[i];
    sum_alpha_n += alpha[i] * n[i];
  }

  double dfdn[16];
  for(int i=0;i<nbins;i++){
    dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2);
  }

    //Error Calculation
  afberr = 0.;
  for(int i=0;i<nbins;i++){
    for(int j=0;j<nbins;j++){
      afberr += covarianceM(i,j) * dfdn[i] * dfdn[j];
    }
  }
  afberr = sqrt(afberr);

    //Calculate Afb
  afb = sum_alpha_n / sum_n;

    //    std::cout<<"AFB = "<<afb<<" "<<afberr<<std::endl;
}


void GetCorrectedAfb_integratewidth_V(TH1D* histogram, TMatrixD &covarianceM, Float_t &afb, Float_t  &afberr){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

    //Get histogram info
  int nbins = histogram->GetNbinsX();
  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1) * histogram->GetBinWidth(i+1);
  }

    //Setup Alpha Vector
  double alpha[16];
  for(int i=0;i<nbins;i++) if(i < nbins/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
  double sum_n = 0.;
  double sum_alpha_n = 0.;
  for(int i=0;i<nbins;i++){
    sum_n += n[i];
    sum_alpha_n += alpha[i] * n[i];
  }

  double dfdn[16];
  for(int i=0;i<nbins;i++){
    dfdn[i] = ( alpha[i] * sum_n - sum_alpha_n ) / pow(sum_n,2);
  }

    //Error Calculation
  afberr = 0.;
  for(int i=0;i<nbins;i++){
    for(int j=0;j<nbins;j++){
      afberr += covarianceM(i,j) * dfdn[i] * dfdn[j] * histogram->GetBinWidth(i+1) * histogram->GetBinWidth(j+1);
    }
  }
  afberr = sqrt(afberr);

    //Calculate Afb
  afb = sum_alpha_n / sum_n;

    //    std::cout<<"AFB = "<<afb<<" "<<afberr<<std::endl;
}


void GetCorrectedAfbBinByBin(TH1D* histogram, TMatrixD &covarianceM, std::vector<double> &myafb, std::vector<double> &myerr, ofstream& second_output_file){

    //Need to calculate AFB and Error for the fully corrected distribution, m_correctE(j,i)

  myafb.clear();
  myerr.clear();

    //Get histogram info
  int nbins = histogram->GetNbinsX();
  const int nbins2 = nbins/2;

  Double_t afbbin[nbins2];
  Double_t afberrbin[nbins2];

  memset( afbbin, 0, sizeof(afbbin) );  //Initialize these arrays to zero
  memset( afberrbin, 0, sizeof(afberrbin) );

  double n[16];
  for(int i=0;i<nbins;i++){
    n[i] = histogram->GetBinContent(i+1);
  }

    //Setup Alpha Vector
  double alpha[16];
  for(int i=0;i<nbins;i++) if(i < nbins/2 ){ alpha[i] = -1;}else{ alpha[i] = 1;}

    //Components of the error calculation
  double sum_n[nbins2];
  double sum_alpha_n[nbins2];
  double sum_n_total = 0.;
  double sum_alpha_n_total = 0.;

  memset( sum_n, 0, sizeof(sum_n) );
  memset( sum_alpha_n, 0, sizeof(sum_alpha_n) );

  for(int i=0;i<nbins2;i++){
    sum_n[i] = n[i] + n[nbins-1-i];
    sum_alpha_n[i] = alpha[i] * n[i] + alpha[nbins-1-i] * n[nbins-1-i];
    sum_n_total += sum_n[i];
    sum_alpha_n_total += sum_alpha_n[i];
  }

  double dfdn[16];
  for(int i=0;i<nbins;i++){
    int k = -999;
    if (i < nbins2) k = i;
    else k = nbins-1-i;
    dfdn[i] = ( alpha[i] * sum_n[k] - sum_alpha_n[k] ) / pow(sum_n[k],2);
  }

    //Error Calculation

  for(int k=0;k<nbins2;k++){
    afberrbin[k] = 0.;
    for(int i=0;i<nbins;i++){
      for(int j=0;j<nbins;j++){
        if( (i==k || i==nbins-1-k ) && (j==k || j==nbins-1-k ) ) {
          afberrbin[k] += covarianceM(i,j) * dfdn[i] * dfdn[j];
              //std::cout<<covarianceM(i,j)<<" "<<dfdn[i]<<" "<<dfdn[j]<<" "<<std::endl;
        }
      }
    }
    afberrbin[k] = sqrt(afberrbin[k]);
    afbbin[k] = sum_alpha_n[k] / sum_n[k];
    std::cout<<k<<" AFB = "<<afbbin[k]<<" +/- "<<afberrbin[k]<<std::endl;
    second_output_file << acceptanceName << " " << observablename << " AFB" << k << ": " << afbbin[k] << " +/- " << afberrbin[k] << std::endl;

    myafb.push_back(afbbin[k]);
    myerr.push_back(afberrbin[k]);
  }
  //double afb = sum_alpha_n_total / sum_n_total;
  //std::cout<<"AFB = "<<afb<<std::endl;
}







void Initialize1DBinning(int iVar){


  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      //xaxislabel="|#eta_{l+}|-|#eta_{l-}|";
      xaxislabel="#Delta|#eta_{#font[12]{l}}|";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{C}^{lep}";
      xbins1D[0]=-2.0; xbins1D[1]=-0.8; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.8; xbins1D[6]=2.0;
      //xbins2Dalt[0]=-2.0; xbins2Dalt[2]=-0.8; xbins2Dalt[4]=-0.4; xbins2Dalt[6]=0.0; xbins2Dalt[8]=0.4; xbins2Dalt[10]=0.8; xbins2Dalt[12]=2.0;
      //xbins2Dalt[1]=-1.4; xbins2Dalt[3]=-0.6; xbins2Dalt[5]=-0.2; xbins2Dalt[7]=0.2; xbins2Dalt[9]=0.6; xbins2Dalt[11]=1.4;
      std::copy(bins_lepChargeAsym, bins_lepChargeAsym + nbinsx2Dalt + 1, xbins2Dalt);
stat_corr[0] = 0.002646;
stat_corr[1] = 0.004095;
stat_corr[2] = 0.005030;
stat_corr[3] = 0.006198;
stat_corr[4] = 0.006825;
stat_corr[5] = 0.008881;
stat_corr[6] = 0.008686;
stat_corr[7] = 0.006642;
stat_corr[8] = 0.006402;
stat_corr[9] = 0.005072;
stat_corr[10] = 0.004155;
stat_corr[11] = 0.002671;
syst_corr[0] = 0.002749;
syst_corr[1] = 0.003118;
syst_corr[2] = 0.004575;
syst_corr[3] = 0.005071;
syst_corr[4] = 0.006178;
syst_corr[5] = 0.005837;
syst_corr[6] = 0.006953;
syst_corr[7] = 0.005181;
syst_corr[8] = 0.004333;
syst_corr[9] = 0.003753;
syst_corr[10] = 0.003053;
syst_corr[11] = 0.002622;
      xmin1D=-2.0;
      xmax1D= 2.0;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="|#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}|";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      //Double_t pi = 3.141592653589793;
      xbins1D[0]=0.0; xbins1D[1]=4.*pi/20.; xbins1D[2]=7.*pi/20.; xbins1D[3]=10.*pi/20.; xbins1D[4]=13.*pi/20.; xbins1D[5]=16.*pi/20.; xbins1D[6]=pi;
      //xbins2Dalt[0]=0.0; xbins2Dalt[2]=4.*pi/20.; xbins2Dalt[4]=7.*pi/20.; xbins2Dalt[6]=10.*pi/20.; xbins2Dalt[8]=13.*pi/20.; xbins2Dalt[10]=16.*pi/20.; xbins2Dalt[12]=pi;
      //xbins2Dalt[1]=2.*pi/20.; xbins2Dalt[3]=5.5*pi/20.; xbins2Dalt[5]=8.5*pi/20.; xbins2Dalt[7]=11.5*pi/20.; xbins2Dalt[9]=14.5*pi/20.; xbins2Dalt[11]=18.*pi/20.;
      std::copy(bins_lepAzimAsym2, bins_lepAzimAsym2 + nbinsx2Dalt + 1, xbins2Dalt);
stat_corr[0] = 0.006411;
stat_corr[1] = 0.006039;
stat_corr[2] = 0.005614;
stat_corr[3] = 0.005493;
stat_corr[4] = 0.005597;
stat_corr[5] = 0.005617;
stat_corr[6] = 0.005613;
stat_corr[7] = 0.005743;
stat_corr[8] = 0.005673;
stat_corr[9] = 0.005683;
stat_corr[10] = 0.005819;
stat_corr[11] = 0.005818;
syst_corr[0] = 0.006339;
syst_corr[1] = 0.005889;
syst_corr[2] = 0.005597;
syst_corr[3] = 0.006652;
syst_corr[4] = 0.004198;
syst_corr[5] = 0.003349;
syst_corr[6] = 0.003054;
syst_corr[7] = 0.003544;
syst_corr[8] = 0.005065;
syst_corr[9] = 0.006177;
syst_corr[10] = 0.007218;
syst_corr[11] = 0.006920;
      xmin1D=0.0;
      xmax1D=pi;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosTheta, bins_lepCosTheta + nbins1D + 1, xbins1D);
stat_corr[0] = 0.014481;
stat_corr[1] = 0.014018;
stat_corr[2] = 0.014256;
stat_corr[3] = 0.014459;
stat_corr[4] = 0.015634;
stat_corr[5] = 0.018516;
syst_corr[0] = 0.038579;
syst_corr[1] = 0.016556;
syst_corr[2] = 0.011227;
syst_corr[3] = 0.011516;
syst_corr[4] = 0.022296;
syst_corr[5] = 0.027031;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosTheta, bins_lepCosTheta + nbins1D + 1, xbins1D);
stat_corr[0] = 0.014476;
stat_corr[1] = 0.013998;
stat_corr[2] = 0.014394;
stat_corr[3] = 0.014567;
stat_corr[4] = 0.015772;
stat_corr[5] = 0.018573;
syst_corr[0] = 0.022806;
syst_corr[1] = 0.015086;
syst_corr[2] = 0.028456;
syst_corr[3] = 0.010661;
syst_corr[4] = 0.018041;
syst_corr[5] = 0.027773;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}#kern[-0.7]{*}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosTheta, bins_lepCosTheta + nbins1D + 1, xbins1D);
stat_corr[0] = 0.010241;
stat_corr[1] = 0.009913;
stat_corr[2] = 0.010139;
stat_corr[3] = 0.010265;
stat_corr[4] = 0.011115;
stat_corr[5] = 0.013120;
syst_corr[0] = 0.025667;
syst_corr[1] = 0.014574;
syst_corr[2] = 0.014193;
syst_corr[3] = 0.008997;
syst_corr[4] = 0.019432;
syst_corr[5] = 0.025580;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Polarization combining positively and negatively charged leptons (CPV, hard-coded distributions)
    case 12:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}^{CPV}";
      acceptanceName="lepCosThetaCPV";
      asymlabel="A_{P}^{CPV}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosTheta, bins_lepCosTheta + nbins1D + 1, xbins1D);
hardcodeddata[0] = 0.512634;
hardcodeddata[1] = 0.486733;
hardcodeddata[2] = 0.500935;
hardcodeddata[3] = 0.516042;
hardcodeddata[4] = 0.491587;
hardcodeddata[5] = 0.492067;
stat_corr[0] = 0.011828;
stat_corr[1] = 0.010230;
stat_corr[2] = 0.009966;
stat_corr[3] = 0.007616;
stat_corr[4] = 0.008056;
stat_corr[5] = 0.009275;
syst_corr[0] = 0.013094;
syst_corr[1] = 0.008539;
syst_corr[2] = 0.012121;
syst_corr[3] = 0.013870;
syst_corr[4] = 0.008778;
syst_corr[5] = 0.012532;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*} #kern[+0.3]{c}os#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.5; xbins1D[2]=-0.2; xbins1D[3]=0.0; xbins1D[4]=0.2; xbins1D[5]=0.5; xbins1D[6]=1.0;
      std::copy(bins_topSpinCorr, bins_topSpinCorr + nbins1D + 1, xbins1D);
stat_corr[0] = 0.009609;
stat_corr[1] = 0.027576;
stat_corr[2] = 0.033787;
stat_corr[3] = 0.033890;
stat_corr[4] = 0.029060;
stat_corr[5] = 0.010769;
syst_corr[0] = 0.009116;
syst_corr[1] = 0.027285;
syst_corr[2] = 0.027023;
syst_corr[3] = 0.031841;
syst_corr[4] = 0.027917;
syst_corr[5] = 0.011565;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapiditydiff_Marco";
      //xaxislabel="|y_{top}|-|y_{tbar}|";
      xaxislabel="#Delta|y_{t}|";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      //xbins1D[0]=-2.0; xbins1D[1]=-0.7; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.7; xbins1D[6]=2.0;
      std::copy(bins_rapiditydiffMarco, bins_rapiditydiffMarco + nbins1D + 1, xbins1D);
stat_corr[0] = 0.003615;
stat_corr[1] = 0.012856;
stat_corr[2] = 0.013327;
stat_corr[3] = 0.013256;
stat_corr[4] = 0.012843;
stat_corr[5] = 0.003611;
syst_corr[0] = 0.003637;
syst_corr[1] = 0.012526;
syst_corr[2] = 0.010109;
syst_corr[3] = 0.015468;
syst_corr[4] = 0.008035;
syst_corr[5] = 0.002589;
      xmin1D=-2.0;
      xmax1D= 2.0;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.7; xbins1D[2]=-0.4; xbins1D[3]=0.0; xbins1D[4]=0.4; xbins1D[5]=0.7; xbins1D[6]=1.0;
      std::copy(bins_topCosTheta, bins_topCosTheta + nbins1D + 1, xbins1D);
stat_corr[0] = 0.015279;
stat_corr[1] = 0.013734;
stat_corr[2] = 0.012428;
stat_corr[3] = 0.012622;
stat_corr[4] = 0.013600;
stat_corr[5] = 0.015356;
syst_corr[0] = 0.019106;
syst_corr[1] = 0.013487;
syst_corr[2] = 0.009754;
syst_corr[3] = 0.012733;
syst_corr[4] = 0.011470;
syst_corr[5] = 0.017153;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos#kern[+0.3]{#varphi}";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos#kern[+0.3]{#varphi}}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosOpeningAngle, bins_lepCosOpeningAngle + nbins1D + 1, xbins1D);
stat_corr[0] = 0.011557;
stat_corr[1] = 0.015963;
stat_corr[2] = 0.014612;
stat_corr[3] = 0.014758;
stat_corr[4] = 0.016970;
stat_corr[5] = 0.016595;
syst_corr[0] = 0.009921;
syst_corr[1] = 0.012321;
syst_corr[2] = 0.009965;
syst_corr[3] = 0.024970;
syst_corr[4] = 0.015402;
syst_corr[5] = 0.037653;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 9:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins1D[0]=-pi; xbins1D[1]=-0.67*pi; xbins1D[2]=-0.33*pi; xbins1D[3]=0.0; xbins1D[4]=0.33*pi; xbins1D[5]=0.67*pi; xbins1D[6]=pi;
      std::copy(bins_lepAzimAsym, bins_lepAzimAsym + nbinsx2Dalt + 1, xbins2Dalt);
stat_corr[0] = 0.002927;
stat_corr[1] = 0.002814;
stat_corr[2] = 0.002817;
stat_corr[3] = 0.002780;
stat_corr[4] = 0.002738;
stat_corr[5] = 0.003100;
stat_corr[6] = 0.003023;
stat_corr[7] = 0.002759;
stat_corr[8] = 0.002784;
stat_corr[9] = 0.002827;
stat_corr[10] = 0.002855;
stat_corr[11] = 0.002925;
syst_corr[0] = 0.003453;
syst_corr[1] = 0.003387;
syst_corr[2] = 0.001736;
syst_corr[3] = 0.002245;
syst_corr[4] = 0.002749;
syst_corr[5] = 0.003002;
syst_corr[6] = 0.002776;
syst_corr[7] = 0.002959;
syst_corr[8] = 0.002058;
syst_corr[9] = 0.001397;
syst_corr[10] = 0.003312;
syst_corr[11] = 0.003486;
      xmin1D=-pi;
      xmax1D= pi;
      break;
    }
    //   Top Asy I
    case 10:
    {
      observablename="top_pseudorapiditydiff_cms";
      xaxislabel="|#eta_{t}|-|#eta_{#bar{t}}|";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      //xbins1D[0]=-2.0; xbins1D[1]=-1.0; xbins1D[2]=-0.5; xbins1D[3]=0.0; xbins1D[4]=0.5; xbins1D[5]=1.0; xbins1D[6]=2.0;
      std::copy(bins_pseudorapiditydiff, bins_pseudorapiditydiff + nbins1D + 1, xbins1D);
stat_corr[0] = 0.005248;
stat_corr[1] = 0.009735;
stat_corr[2] = 0.010001;
stat_corr[3] = 0.010019;
stat_corr[4] = 0.009636;
stat_corr[5] = 0.005234;
syst_corr[0] = 0.008124;
syst_corr[1] = 0.015570;
syst_corr[2] = 0.007862;
syst_corr[3] = 0.008307;
syst_corr[4] = 0.008855;
syst_corr[5] = 0.004897;
      xmin1D=-2.0;
      xmax1D= 2.0;
      break;
    }
    //   Top Asy II
    case 11:
    {
      observablename="top_rapiditydiff_cms";
      xaxislabel="(y_{top}-y_{#bar{t}})(y_{top}+y_{#bar{t}})";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      //xbins1D[0]=-2.0; xbins1D[1]=-0.8; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.8; xbins1D[6]=2.0;
      std::copy(bins_rapiditydiff, bins_rapiditydiff + nbins1D + 1, xbins1D);
stat_corr[0] = 0.004327;
stat_corr[1] = 0.006495;
stat_corr[2] = 0.011174;
stat_corr[3] = 0.011144;
stat_corr[4] = 0.006467;
stat_corr[5] = 0.004310;
syst_corr[0] = 0.004502;
syst_corr[1] = 0.005017;
syst_corr[2] = 0.009560;
syst_corr[3] = 0.015023;
syst_corr[4] = 0.004748;
syst_corr[5] = 0.004847;
      xmin1D=-2.0;
      xmax1D= 2.0;
      break;
    }
    default:
    {
      std::cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinning(int iVar){

  std::copy( ybinsmtt, ybinsmtt+nbinsy2D+1, ybins2D );
  //ybins2D[0]=0.0; ybins2D[1]=430; ybins2D[2]=530.0; ybins2D[3]=1200.0;  // KIT binning
  //ybins2D[0]=0.0; ybins2D[1]=410; ybins2D[2]=510.0; ybins2D[3]=800.0;   // SnT binning

  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="#Delta|#eta_{#font[12]{l}}|";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{C}^{lep}";
      xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.8; xbins2D[6]=2.0;
      // xbins2Dalt[0]=-2.0; xbins2Dalt[2]=-0.8; xbins2Dalt[4]=-0.4; xbins2Dalt[6]=0.0; xbins2Dalt[8]=0.4; xbins2Dalt[10]=0.8; xbins2Dalt[12]=2.0;
      // xbins2Dalt[1]=-1.4; xbins2Dalt[3]=-0.6; xbins2Dalt[5]=-0.2; xbins2Dalt[7]=0.2; xbins2Dalt[9]=0.6; xbins2Dalt[11]=1.4;
	  std::copy( bins_lepChargeAsym, bins_lepChargeAsym+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2Dalt[0];
      xmax=xbins2Dalt[12];
stat_corr[0] = 0.013835;
stat_corr[1] = 0.015634;
stat_corr[2] = 0.015603;
syst_corr[0] = 0.012147;
syst_corr[1] = 0.022843;
syst_corr[2] = 0.015408;
theory_default[0] = 2.2854757178882583E-003;
theory_default[1] = 5.9512702066694088E-003;
theory_default[2] = 1.1139752121814698E-002;
theory_scaledown[0] = 2.4501106513824173E-003;
theory_scaledown[1] = 6.1209924089345499E-003;
theory_scaledown[2] = 1.1210203897297082E-002;
theory_scaleup[0] = 2.1609926593942857E-003;
theory_scaleup[1] = 5.8106644384500646E-003;
theory_scaleup[2] = 1.1051329824914245E-002;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="|#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}|";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      //Double_t pi = 3.141592653589793;
      xbins2D[0]=0.0; xbins2D[1]=4.*pi/20.; xbins2D[2]=7.*pi/20.; xbins2D[3]=10.*pi/20.; xbins2D[4]=13.*pi/20.; xbins2D[5]=16.*pi/20.; xbins2D[6]=pi;
      // xbins2Dalt[0]=0.0; xbins2Dalt[2]=4.*pi/20.; xbins2Dalt[4]=7.*pi/20.; xbins2Dalt[6]=10.*pi/20.; xbins2Dalt[8]=13.*pi/20.; xbins2Dalt[10]=16.*pi/20.; xbins2Dalt[12]=pi;
      // xbins2Dalt[1]=2.*pi/20.; xbins2Dalt[3]=5.5*pi/20.; xbins2Dalt[5]=8.5*pi/20.; xbins2Dalt[7]=11.5*pi/20.; xbins2Dalt[9]=14.5*pi/20.; xbins2Dalt[11]=18.*pi/20.;
	  std::copy( bins_lepAzimAsym2, bins_lepAzimAsym2+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2Dalt[0];
      xmax=xbins2Dalt[12];
stat_corr[0] = 0.013555;
stat_corr[1] = 0.016247;
stat_corr[2] = 0.017365;
syst_corr[0] = 0.011587;
syst_corr[1] = 0.013935;
syst_corr[2] = 0.026465;
theory_default[0] = -6.9732520956928390E-002;
theory_default[1] = 0.13757572640895188;
theory_default[2] = 0.32220011359666922;
theory_scaledown[0] = -7.7390640843771508E-002;
theory_scaledown[1] = 0.12885142091447371     ;
theory_scaledown[2] = 0.31053693999666665     ;
theory_scaleup[0] = -6.2696189297629046E-002;
theory_scaleup[1] = 0.14456051533545783     ;
theory_scaleup[2] = 0.33152342645249983   ;
uncorr_default[0] =  7.9717858494246896E-002;
uncorr_default[1] = 0.21583763408365253     ;
uncorr_default[2] =  0.35599189667219111    ;
uncorr_scaledown[0] =  7.3675604395192940E-002;
uncorr_scaledown[1] = 0.20849165447785850     ;
uncorr_scaledown[2] = 0.34509694460959389     ;
uncorr_scaleup[0] =  8.4394038380300018E-002;
uncorr_scaleup[1] = 0.22114286723734239     ;
uncorr_scaleup[2] = 0.36488777812911555     ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.027648;
stat_corr[1] = 0.026979;
stat_corr[2] = 0.023121;
syst_corr[0] = 0.024340;
syst_corr[1] = 0.040043;
syst_corr[2] = 0.039233;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.027791;
stat_corr[1] = 0.027502;
stat_corr[2] = 0.023537;
syst_corr[0] = 0.034023;
syst_corr[1] = 0.040907;
syst_corr[2] = 0.036260;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}#kern[-0.7]{*}";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.019611;
stat_corr[1] = 0.019257;
stat_corr[2] = 0.016504;
syst_corr[0] = 0.024603;
syst_corr[1] = 0.038421;
syst_corr[2] = 0.034595;
theory_default[0] = 7.1082369726289908E-004/2.;
theory_default[1] = 2.3378532348680807E-003/2.;
theory_default[2] = 6.3313438371528459E-003/2.;
theory_scaledown[0] = 6.6874849711424591E-005/2.;
theory_scaledown[1] = 1.4348503548028096E-003/2.;
theory_scaledown[2] = 5.1933101595061102E-003/2.;
theory_scaleup[0] = 1.5657140676943848E-003/2.;
theory_scaleup[1] = 3.5699273929340894E-003/2.;
theory_scaleup[2] = 7.9729622808530475E-003/2.;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons (CPV, hard-coded distributions)
    case 12:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}^{CPV}";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepCosThetaCPV";
      asymlabel="A_{P}^{CPV}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
    std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
hardcodeddata[0] = -0.028583;
hardcodeddata[1] = 0.048836;
hardcodeddata[2] = -0.017654;
stat_corr[0] = 0.016042;
stat_corr[1] = 0.016699;
stat_corr[2] = 0.011597;
syst_corr[0] = 0.013420;
syst_corr[1] = 0.012704;
syst_corr[2] = 0.015862;
theory_default[0] = 0.;
theory_default[1] = 0.;
theory_default[2] = 0.;
theory_scaledown[0] = 0.00001;
theory_scaledown[1] = 0.00001;
theory_scaledown[2] = 0.00001;
theory_scaleup[0] = 0.;
theory_scaleup[1] = 0.;
theory_scaleup[2] = 0.;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*} #kern[+0.3]{c}os#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.5; xbins2D[2]=-0.2; xbins2D[3]=0.0; xbins2D[4]=0.2; xbins2D[5]=0.5; xbins2D[6]=1.0;
	  std::copy( bins_topSpinCorr, bins_topSpinCorr+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.031606;
stat_corr[1] = 0.028903;
stat_corr[2] = 0.026475;
syst_corr[0] = 0.030814;
syst_corr[1] = 0.031671;
syst_corr[2] = 0.029691;
theory_default[0] = 0.46750205191882088     /-4.;
theory_default[1] = 0.31317519190246929     /-4.;
theory_default[2] = 0.12261795699405617     /-4.;
theory_scaledown[0] = 0.46319347417725298     /-4.;
theory_scaledown[1] = 0.31165874096446244     /-4.;
theory_scaledown[2] = 0.12883792536745478     /-4.;
theory_scaleup[0] = 0.47011975407220463     /-4.;
theory_scaleup[1] = 0.31449303134754297     /-4.;
theory_scaleup[2] = 0.11864642708770839     /-4.;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapiditydiff_Marco";
      //xaxislabel="|y_{top}|-|y_{tbar}|";
      xaxislabel="#Delta|y_{t}|";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      // xbins2D[0]=-2.0; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiffMarco, bins_rapiditydiffMarco+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.030726;
stat_corr[1] = 0.028637;
stat_corr[2] = 0.021358;
syst_corr[0] = 0.020164;
syst_corr[1] = 0.016851;
syst_corr[2] = 0.012830;
theory_default[0] = 0.8191/100.;
theory_default[1] = 1.2336/100.;
theory_default[2] = 1.4604/100.;
theory_scaledown[0] = 0.8514/100.;
theory_scaledown[1] = 1.2675/100.;
theory_scaledown[2] = 1.4769/100.;
theory_scaleup[0] = 0.7782/100.;
theory_scaleup[1] = 1.2035/100.;
theory_scaleup[2] = 1.4331/100.;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.7; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.7; xbins2D[6]=1.0;
	  std::copy( bins_topCosTheta, bins_topCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.033664;
stat_corr[1] = 0.027623;
stat_corr[2] = 0.021358;
syst_corr[0] = 0.026631;
syst_corr[1] = 0.020886;
syst_corr[2] = 0.020771;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos#kern[+0.3]{#varphi}";
	  yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos#kern[+0.3]{#varphi}}";
      std::copy( bins_lepCosOpeningAngle, bins_lepCosOpeningAngle+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.023536;
stat_corr[1] = 0.030750;
stat_corr[2] = 0.026147;
syst_corr[0] = 0.023502;
syst_corr[1] = 0.021852;
syst_corr[2] = 0.020929;
theory_default[0] = -0.37967779291362413     /-2.;
theory_default[1] = -0.20308561122601917     /-2.;
theory_default[2] =  -7.9999796989889432E-002/-2.;
theory_scaledown[0] = -0.37980026787542420     /-2.;
theory_scaledown[1] = -0.20577418488172378     /-2.;
theory_scaledown[2] =  -8.4616320034969911E-002/-2.;
theory_scaleup[0] = -0.37951909898526098     /-2.;
theory_scaleup[1] = -0.20081880635556171     /-2.;
theory_scaleup[2] =  -7.6374174747447288E-002/-2.;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 9:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins2D[0]=-pi; xbins2D[1]=-0.67*pi; xbins2D[2]=-0.33*pi; xbins2D[3]=0.0; xbins2D[4]=0.33*pi; xbins2D[5]=0.67*pi; xbins2D[6]=pi;
	  std::copy( bins_lepAzimAsym, bins_lepAzimAsym+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2Dalt[0];
      xmax=xbins2Dalt[12];
stat_corr[0] = 0.013350;
stat_corr[1] = 0.015831;
stat_corr[2] = 0.016306;
syst_corr[0] = 0.012829;
syst_corr[1] = 0.019967;
syst_corr[2] = 0.011459;
      break;
    }
    //   Top Asy I
    case 10:
    {
      observablename="top_pseudorapiditydiff_cms";
      xaxislabel="|#eta_{t}|-|#eta_{#bar{t}}|";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      // xbins2D[0]=-2.0; xbins2D[1]=-1.0; xbins2D[2]=-0.5; xbins2D[3]=0.0; xbins2D[4]=0.5; xbins2D[5]=1.0; xbins2D[6]=2.0;
	  std::copy( bins_pseudorapiditydiff, bins_pseudorapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.034780;
stat_corr[1] = 0.029588;
stat_corr[2] = 0.022018;
syst_corr[0] = 0.035417;
syst_corr[1] = 0.017418;
syst_corr[2] = 0.015974;
      break;
    }
    //   Top Asy II
    case 11:
    {
      observablename="top_rapiditydiff_cms";
      xaxislabel="(y_{top}-y_{#bar{t}})(y_{top}+y_{#bar{t}})";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      // xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.8; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiff, bins_rapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.031085;
stat_corr[1] = 0.026499;
stat_corr[2] = 0.021183;
syst_corr[0] = 0.022211;
syst_corr[1] = 0.017211;
syst_corr[2] = 0.014910;
      break;
    }
    default:
    {
      std::cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinningttpt(int iVar){

  std::copy( ybinsttpt, ybinsttpt+nbinsy2D+1, ybins2D );
  //ybins2D[0]=0.0; ybins2D[1]=41.0; ybins2D[2]=92.0; ybins2D[3]=300.0;  //KIT
  //ybins2D[0]=0.0; ybins2D[1]=24; ybins2D[2]=52.0; ybins2D[3]=100.0;  //SnT

  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="#Delta|#eta_{#font[12]{l}}|";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{C}^{lep}";
      xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.8; xbins2D[6]=2.0;
      std::copy(bins_lepChargeAsym, bins_lepChargeAsym + nbinsx2Dalt + 1, xbins2Dalt);
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.009680;
stat_corr[1] = 0.015919;
stat_corr[2] = 0.018381;
syst_corr[0] = 0.008818;
syst_corr[1] = 0.018239;
syst_corr[2] = 0.010411;
theory_default[0] =  7.1593856807750456E-003;
theory_default[1] = -2.8365229918471788E-003;
theory_default[2] = -2.4853842372138179E-004;
theory_scaledown[0] =  7.4219358633036405E-003;
theory_scaledown[1] = -2.6480402174896131E-003;
theory_scaledown[2] = -1.8538677370730610E-004;
theory_scaleup[0] =  6.9446338508168184E-003;
theory_scaleup[1] = -3.0510257170354711E-003;
theory_scaleup[2] = -3.4258795512555877E-004;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="|#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}|";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      //Double_t pi = 3.141592653589793;
      xbins2D[0]=0.0; xbins2D[1]=4.*pi/20.; xbins2D[2]=7.*pi/20.; xbins2D[3]=10.*pi/20.; xbins2D[4]=13.*pi/20.; xbins2D[5]=16.*pi/20.; xbins2D[6]=pi;
	  std::copy( bins_lepAzimAsym2, bins_lepAzimAsym2+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.009665;
stat_corr[1] = 0.016347;
stat_corr[2] = 0.017955;
syst_corr[0] = 0.015918;
syst_corr[1] = 0.016934;
syst_corr[2] = 0.019506;
theory_default[0] = 0.12572265215619549;
theory_default[1] = 0.17730999154723914;
theory_default[2] = 0.10262414675791232;
theory_scaledown[0] = 0.11750437273931535     ;
theory_scaledown[1] = 0.18068890055195069     ;
theory_scaledown[2] = 0.10458438626234656     ;
theory_scaleup[0] = 0.13202651204583515     ;
theory_scaleup[1] = 0.17451905187496153     ;
theory_scaleup[2] =  0.10097884989366024     ;
uncorr_default[0] = 0.22220976181669910     ;
uncorr_default[1] = 0.26857894246502284     ;
uncorr_default[2] = 0.16656665689188777     ;
uncorr_scaledown[0] = 0.21643512364946638     ;
uncorr_scaledown[1] = 0.27174929871286002     ;
uncorr_scaledown[2] = 0.16833912907072063     ;
uncorr_scaleup[0] = 0.22631675529417361     ;
uncorr_scaleup[1] = 0.26588952844158220     ;
uncorr_scaleup[2] = 0.16516067482158125     ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.017159;
stat_corr[1] = 0.031246;
stat_corr[2] = 0.033723;
syst_corr[0] = 0.030499;
syst_corr[1] = 0.029990;
syst_corr[2] = 0.032914;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.017274;
stat_corr[1] = 0.032483;
stat_corr[2] = 0.034299;
syst_corr[0] = 0.022402;
syst_corr[1] = 0.043310;
syst_corr[2] = 0.059827;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}#kern[-0.7]{*}";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.012192;
stat_corr[1] = 0.022546;
stat_corr[2] = 0.024068;
syst_corr[0] = 0.025264;
syst_corr[1] = 0.035629;
syst_corr[2] = 0.040167;
theory_default[0] =  4.0906016806477302E-003/2.;
theory_default[1] = -2.5301124136594020E-003/2.;
theory_default[2] = -1.5738490152040310E-003/2.;
theory_scaledown[0] =  3.1011380834262948E-003/2.;
theory_scaledown[1] = -1.8882542671156361E-003/2.;
theory_scaledown[2] = -1.2015870115943490E-003/2.;
theory_scaleup[0] =  5.4119065998264569E-003/2.;
theory_scaleup[1] = -3.2501230390085197E-003/2.;
theory_scaleup[2] = -1.9913120462740344E-003/2.;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons (CPV, hard-coded distributions)
    case 12:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}^{CPV}";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepCosThetaCPV";
      asymlabel="A_{P}^{CPV}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
    std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
hardcodeddata[0] = 0.013834;
hardcodeddata[1] = -0.016522;
hardcodeddata[2] = -0.011906;
stat_corr[0] = 0.010337;
stat_corr[1] = 0.019123;
stat_corr[2] = 0.016800;
syst_corr[0] = 0.008671;
syst_corr[1] = 0.012947;
syst_corr[2] = 0.026431;
theory_default[0] = 0.;
theory_default[1] = 0.;
theory_default[2] = 0.;
theory_scaledown[0] = 0.00001;
theory_scaledown[1] = 0.00001;
theory_scaledown[2] = 0.00001;
theory_scaleup[0] = 0.;
theory_scaleup[1] = 0.;
theory_scaleup[2] = 0.;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*} #kern[+0.3]{c}os#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.5; xbins2D[2]=-0.2; xbins2D[3]=0.0; xbins2D[4]=0.2; xbins2D[5]=0.5; xbins2D[6]=1.0;
	  std::copy( bins_topSpinCorr, bins_topSpinCorr+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.019954;
stat_corr[1] = 0.038313;
stat_corr[2] = 0.043456;
syst_corr[0] = 0.018771;
syst_corr[1] = 0.039446;
syst_corr[2] = 0.037553;
theory_default[0] = 0.32575069797174833     /-4.;
theory_default[1] = 0.30172042152922113     /-4.;
theory_default[2] = 0.16011095908856965     /-4.;
theory_scaledown[0] = 0.33275360908640378     /-4.;
theory_scaledown[1] = 0.30486960295536952     /-4.;
theory_scaledown[2] = 0.16173782686029431     /-4.;
theory_scaleup[0] = 0.32107426318799120     /-4.;
theory_scaleup[1] = 0.29798436513324439     /-4.;
theory_scaleup[2] = 0.15791658414606186     /-4.;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapiditydiff_Marco";
      //xaxislabel="|y_{top}|-|y_{tbar}|";
      xaxislabel="#Delta|y_{t}|";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiffMarco, bins_rapiditydiffMarco+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.017959;
stat_corr[1] = 0.033168;
stat_corr[2] = 0.035523;
syst_corr[0] = 0.013415;
syst_corr[1] = 0.028209;
syst_corr[2] = 0.032569;
theory_default[0] = 1.2727/100.;
theory_default[1] = -0.4664/100.;
theory_default[2] = -0.1398/100.;
theory_scaledown[0] = 1.3310/100.;
theory_scaledown[1] = -0.4397/100.;
theory_scaledown[2] = -0.1378/100.;
theory_scaleup[0] = 1.2126/100.;
theory_scaleup[1] = -0.4905/100.;
theory_scaleup[2] = -0.1579/100.;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.7; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.7; xbins2D[6]=1.0;
	  std::copy( bins_topCosTheta, bins_topCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.017956;
stat_corr[1] = 0.032204;
stat_corr[2] = 0.033988;
syst_corr[0] = 0.012867;
syst_corr[1] = 0.026835;
syst_corr[2] = 0.029921;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos#kern[+0.3]{#varphi}";
	  yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos#kern[+0.3]{#varphi}}";
      std::copy( bins_lepCosOpeningAngle, bins_lepCosOpeningAngle+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.016582;
stat_corr[1] = 0.030646;
stat_corr[2] = 0.035281;
syst_corr[0] = 0.017349;
syst_corr[1] = 0.019668;
syst_corr[2] = 0.037118;
theory_default[0] = -0.23035649296477598     /-2.;
theory_default[1] = -0.21600252525269378     /-2.;
theory_default[2] = -0.15394842051237889     /-2.;
theory_scaledown[0] = -0.23572256019027907     /-2.;
theory_scaledown[1] = -0.21753199627949082     /-2.;
theory_scaledown[2] = -0.15506885779939464     /-2.;
theory_scaleup[0] = -0.22657554420214174     /-2.;
theory_scaleup[1] = -0.21402905918130047     /-2.;
theory_scaleup[2] = -0.15254738928307032     /-2.;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 9:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins2D[0]=-pi; xbins2D[1]=-0.67*pi; xbins2D[2]=-0.33*pi; xbins2D[3]=0.0; xbins2D[4]=0.33*pi; xbins2D[5]=0.67*pi; xbins2D[6]=pi;
	  std::copy( bins_lepAzimAsym, bins_lepAzimAsym+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2Dalt[0];
      xmax=xbins2Dalt[12];
stat_corr[0] = 0.009507;
stat_corr[1] = 0.016015;
stat_corr[2] = 0.017831;
syst_corr[0] = 0.008425;
syst_corr[1] = 0.018556;
syst_corr[2] = 0.012965;
      break;
    }
    //   Top Asy I
    case 10:
    {
      observablename="top_pseudorapiditydiff_cms";
      xaxislabel="|#eta_{t}|-|#eta_{#bar{t}}|";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-1.0; xbins2D[2]=-0.5; xbins2D[3]=0.0; xbins2D[4]=0.5; xbins2D[5]=1.0; xbins2D[6]=2.0;
	  std::copy( bins_pseudorapiditydiff, bins_pseudorapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.018227;
stat_corr[1] = 0.034196;
stat_corr[2] = 0.036488;
syst_corr[0] = 0.013098;
syst_corr[1] = 0.027217;
syst_corr[2] = 0.040665;
      break;
    }
    //   Top Asy II
    case 11:
    {
      observablename="top_rapiditydiff_cms";
      xaxislabel="(y_{top}-y_{#bar{t}})(y_{top}+y_{#bar{t}})";
      yaxislabel="p_{T}^{t#bar{t}}";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.8; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiff, bins_rapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.017068;
stat_corr[1] = 0.030759;
stat_corr[2] = 0.033878;
syst_corr[0] = 0.012184;
syst_corr[1] = 0.024887;
syst_corr[2] = 0.032782;
      break;
    }
    default:
    {
      std::cout<<"Set the variable switch";
    }
  }
}

void Initialize2DBinningttrapidity2(int iVar){

  std::copy( ybinsttrapidity2, ybinsttrapidity2+nbinsy2D+1, ybins2D );
  //ybins2D[0]=0.0; ybins2D[1]=0.34; ybins2D[2]=0.75; ybins2D[3]=1.5;  //KIT
  //ybins2D[0]=0.0; ybins2D[1]=0.3; ybins2D[2]=0.7; ybins2D[3]=1.5;  //SnT

  switch (iVar)
  {
    //   Lepton Charge Asymmetry
    case 0:
    {
      observablename="lep_charge_asymmetry";
      xaxislabel="#Delta|#eta_{#font[12]{l}}|";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{C}^{lep}";
      xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.8; xbins2D[6]=2.0;
      std::copy(bins_lepChargeAsym, bins_lepChargeAsym + nbinsx2Dalt + 1, xbins2Dalt);
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.015112;
stat_corr[1] = 0.015118;
stat_corr[2] = 0.012636;
syst_corr[0] = 0.010533;
syst_corr[1] = 0.011198;
syst_corr[2] = 0.007772;
theory_default[0] = 1.8121201623237826E-003;
theory_default[1] = 3.8289879454593279E-003;
theory_default[2] = 1.0581939026885376E-002;
theory_scaledown[0] = 1.8802737680046660E-003;
theory_scaledown[1] = 3.9471877447466294E-003;
theory_scaledown[2] = 1.0942112950465519E-002;
theory_scaleup[0] = 1.7149721049203807E-003;
theory_scaleup[1] = 3.7381661697158614E-003;
theory_scaleup[2] = 1.0282819076233869E-002;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="|#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}|";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      //Double_t pi = 3.141592653589793;
      xbins2D[0]=0.0; xbins2D[1]=4.*pi/20.; xbins2D[2]=7.*pi/20.; xbins2D[3]=10.*pi/20.; xbins2D[4]=13.*pi/20.; xbins2D[5]=16.*pi/20.; xbins2D[6]=pi;
	  std::copy( bins_lepAzimAsym2, bins_lepAzimAsym2+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.014424;
stat_corr[1] = 0.015507;
stat_corr[2] = 0.013076;
syst_corr[0] = 0.011797;
syst_corr[1] = 0.015084;
syst_corr[2] = 0.015116;
theory_default[0] = 0.12143852540849276;
theory_default[1] = 0.11773003778632196;
theory_default[2] = 9.7956675608496635E-002;
theory_scaledown[0] = 0.11071632577339391     ;
theory_scaledown[1] = 0.10574476757453949     ;
theory_scaledown[2] =  8.8071804415588872E-002;
theory_scaleup[0] = 0.13021384877415476     ;
theory_scaleup[1] = 0.12634949704727166     ;
theory_scaleup[2] = 0.10551368688593787     ;
uncorr_default[0] = 0.21429142505201348     ;
uncorr_default[1] = 0.20930284910896546     ;
uncorr_default[2] =   0.19167062791048511   ;
uncorr_scaledown[0] = 0.20514123507198134     ;
uncorr_scaledown[1] = 0.19911579639863697     ;
uncorr_scaledown[2] = 0.18327255887038590     ;
uncorr_scaleup[0] = 0.22113149572428242     ;
uncorr_scaleup[1] = 0.21632459729578701     ;
uncorr_scaleup[2] = 0.19804471163300374     ;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.024086;
stat_corr[1] = 0.025472;
stat_corr[2] = 0.022283;
syst_corr[0] = 0.035638;
syst_corr[1] = 0.039748;
syst_corr[2] = 0.033703;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.023729;
stat_corr[1] = 0.025769;
stat_corr[2] = 0.022797;
syst_corr[0] = 0.040623;
syst_corr[1] = 0.025322;
syst_corr[2] = 0.035251;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}#kern[-0.7]{*}";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.016920;
stat_corr[1] = 0.018117;
stat_corr[2] = 0.015955;
syst_corr[0] = 0.034345;
syst_corr[1] = 0.029422;
syst_corr[2] = 0.030605;
theory_default[0] = 3.8492364883239345E-003/2.;
theory_default[1] = 3.5952321772769460E-003/2.;
theory_default[2] = 2.1670253761306621E-003/2.;
theory_scaledown[0] = 2.8415728465513968E-003/2.;
theory_scaledown[1] = 2.6445453653837435E-003/2.;
theory_scaledown[2] = 1.5069758359409560E-003/2.;
theory_scaleup[0] = 5.2646488199438613E-003/2.;
theory_scaleup[1] = 4.9261986930265447E-003/2.;
theory_scaleup[2] = 3.0961161700422133E-003/2.;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons (CPV, hard-coded distributions)
    case 12:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l}}^{CPV}";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepCosThetaCPV";
      asymlabel="A_{P}^{CPV}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
    std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
hardcodeddata[0] = 0.006452;
hardcodeddata[1] = 0.011011;
hardcodeddata[2] = -0.016149;
stat_corr[0] = 0.020541;
stat_corr[1] = 0.015429;
stat_corr[2] = 0.011178;
syst_corr[0] = 0.015563;
syst_corr[1] = 0.016608;
syst_corr[2] = 0.010332;
theory_default[0] = 0.;
theory_default[1] = 0.;
theory_default[2] = 0.;
theory_scaledown[0] = 0.00001;
theory_scaledown[1] = 0.00001;
theory_scaledown[2] = 0.00001;
theory_scaleup[0] = 0.;
theory_scaleup[1] = 0.;
theory_scaleup[2] = 0.;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*} #kern[+0.3]{c}os#kern[+0.3]{#theta}_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*}";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.5; xbins2D[2]=-0.2; xbins2D[3]=0.0; xbins2D[4]=0.2; xbins2D[5]=0.5; xbins2D[6]=1.0;
	  std::copy( bins_topSpinCorr, bins_topSpinCorr+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.027759;
stat_corr[1] = 0.027664;
stat_corr[2] = 0.026640;
syst_corr[0] = 0.026837;
syst_corr[1] = 0.023370;
syst_corr[2] = 0.028555;
theory_default[0] = 0.33658157068751166     /-4.;
theory_default[1] = 0.33554839556169519     /-4.;
theory_default[2] = 0.29448759496760302     /-4.;
theory_scaledown[0] = 0.33745719476246222     /-4.;
theory_scaledown[1] = 0.33664705475432627     /-4.;
theory_scaledown[2] = 0.29864735454600722     /-4.;
theory_scaleup[0] = 0.33749451185900098     /-4.;
theory_scaleup[1] = 0.33580189737431837     /-4.;
theory_scaleup[2] = 0.29105379660190372     /-4.;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapiditydiff_Marco";
      //xaxislabel="|y_{top}|-|y_{tbar}|";
      xaxislabel="#Delta|y_{t}|";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiffMarco, bins_rapiditydiffMarco+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.031062;
stat_corr[1] = 0.027001;
stat_corr[2] = 0.022091;
syst_corr[0] = 0.012456;
syst_corr[1] = 0.015659;
syst_corr[2] = 0.024084;
theory_default[0] = 0.3022/100.;
theory_default[1] = 0.7980/100.;
theory_default[2] = 1.9331/100.;
theory_scaledown[0] = 0.3248/100.;
theory_scaledown[1] = 0.8098/100.;
theory_scaledown[2] = 1.9759/100.;
theory_scaleup[0] = 0.2848/100.;
theory_scaleup[1] = 0.7728/100.;
theory_scaleup[2] = 1.8865/100.;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.7; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.7; xbins2D[6]=1.0;
	  std::copy( bins_topCosTheta, bins_topCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.034300;
stat_corr[1] = 0.023391;
stat_corr[2] = 0.020632;
syst_corr[0] = 0.029624;
syst_corr[1] = 0.017308;
syst_corr[2] = 0.022987;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos#kern[+0.3]{#varphi}";
	  yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos#kern[+0.3]{#varphi}}";
      std::copy( bins_lepCosOpeningAngle, bins_lepCosOpeningAngle+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.023341;
stat_corr[1] = 0.023710;
stat_corr[2] = 0.021224;
syst_corr[0] = 0.021092;
syst_corr[1] = 0.019302;
syst_corr[2] = 0.022364;
theory_default[0] = -0.24487937806196763     /-2.;
theory_default[1] = -0.24549998678247889     /-2.;
theory_default[2] = -0.22152212699184887     /-2.;
theory_scaledown[0] = -0.24775792274269123     /-2.;
theory_scaledown[1] = -0.24848923053444055     /-2.;
theory_scaledown[2] = -0.22678986082154506     /-2.;
theory_scaleup[0] = -0.24350240482743912     /-2.;
theory_scaleup[1] = -0.24370571689956300     /-2.;
theory_scaleup[2] = -0.21735053981565339     /-2.;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 9:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins2D[0]=-pi; xbins2D[1]=-0.67*pi; xbins2D[2]=-0.33*pi; xbins2D[3]=0.0; xbins2D[4]=0.33*pi; xbins2D[5]=0.67*pi; xbins2D[6]=pi;
	  std::copy( bins_lepAzimAsym, bins_lepAzimAsym+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2Dalt[0];
      xmax=xbins2Dalt[12];
stat_corr[0] = 0.014136;
stat_corr[1] = 0.015129;
stat_corr[2] = 0.012711;
syst_corr[0] = 0.007671;
syst_corr[1] = 0.010124;
syst_corr[2] = 0.008784;
      break;
    }
    //   Top Asy I
    case 10:
    {
      observablename="top_pseudorapiditydiff_cms";
      xaxislabel="|#eta_{t}|-|#eta_{#bar{t}}|";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-1.0; xbins2D[2]=-0.5; xbins2D[3]=0.0; xbins2D[4]=0.5; xbins2D[5]=1.0; xbins2D[6]=2.0;
	  std::copy( bins_pseudorapiditydiff, bins_pseudorapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.032537;
stat_corr[1] = 0.026447;
stat_corr[2] = 0.021911;
syst_corr[0] = 0.032061;
syst_corr[1] = 0.023732;
syst_corr[2] = 0.016591;
      break;
    }
    //   Top Asy II
    case 11:
    {
      observablename="top_rapiditydiff_cms";
      xaxislabel="(y_{top}-y_{#bar{t}})(y_{top}+y_{#bar{t}})";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.8; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiff, bins_rapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
stat_corr[0] = 0.028827;
stat_corr[1] = 0.025497;
stat_corr[2] = 0.025246;
syst_corr[0] = 0.018089;
syst_corr[1] = 0.014869;
syst_corr[2] = 0.030576;
      break;
    }
    default:
    {
      std::cout<<"Set the variable switch";
    }
  }
}

void fillUnderOverFlow(TH1D *h1, float value, double weight, double Nsolns)
{
  double min = h1->GetXaxis()->GetXmin();
  double max = h1->GetXaxis()->GetXmax();

  if (value >= max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value <= min) value = h1->GetBinCenter(1);

  int bin_number = h1->FindBin(value);
  double orig_content = h1->GetBinContent(bin_number);
  double orig_error = h1->GetBinError(bin_number);

  //h1->Fill(value, weight);
  h1->SetBinContent( bin_number, orig_content+weight );
  h1->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*Nsolns ) );
}

//--------------------------------------------------------------------

void fillUnderOverFlow(TH2 *h2, float xvalue, float yvalue, double weight, double Nsolns)
{
  double maxx = h2->GetXaxis()->GetXmax();
  double minx = h2->GetXaxis()->GetXmin();
  double maxy = h2->GetYaxis()->GetXmax();
  double miny = h2->GetYaxis()->GetXmin();

  if (xvalue >= maxx) xvalue = h2->GetXaxis()->GetBinCenter(h2->GetNbinsX());
  if (xvalue <= minx) xvalue = h2->GetXaxis()->GetBinCenter(1);
  if (yvalue >= maxy) yvalue = h2->GetYaxis()->GetBinCenter(h2->GetNbinsY());
  if (yvalue <= miny) yvalue = h2->GetYaxis()->GetBinCenter(1);

  int bin_number = h2->FindBin(xvalue,yvalue);
  double orig_content = h2->GetBinContent(bin_number);
  double orig_error = h2->GetBinError(bin_number);

  //h2->Fill(xvalue, yvalue, weight);
  h2->SetBinContent( bin_number, orig_content+weight );
  h2->SetBinError( bin_number, sqrt( orig_error*orig_error + weight*weight*Nsolns ) );
}




//Given a TH2D and the x,y coordinates of a point, this function tells you what
//unwrapped bin that point would go into.
int getUnwrappedBin (TH2 *h2, double xval, double yval)
{
  int nx = h2->GetNbinsX();
  int ny = h2->GetNbinsY();
  int nxroot = nx+2;
  int nyroot = ny+2;

  int rootbin = h2->FindBin(xval, yval);

  //Deal with y-under/overflow, and corner cases
  if(rootbin < nxroot)  rootbin += nxroot;
  else if (rootbin >= (nyroot-1)*nxroot)  rootbin -= nxroot;

  //Deal with x-under/overflow
  if(rootbin%nxroot == 0) rootbin+=1;
  else if((rootbin+1)%nxroot == 0) rootbin -= 1;

  int unwrappedbin = rootbin - nx - 2*(rootbin/nxroot);

  if (unwrappedbin > nx*ny || unwrappedbin < 1) std::cout << "ERROR: unwrapped bin " << unwrappedbin << " doesn't exist! (ROOT bin " << rootbin << ")." << std::endl;
  return unwrappedbin;
}






// This is copied from the KIT group, with one important modification: to make a TH2 into a TH1,
// they concatenate columns, while we concatenate rows!
void unwrap2dhisto(TH2* h2, TH1* h1)
{
  int n=0;
  for(int y=1; y<=h2->GetNbinsY(); y++)
  {
	for(int x=1; x<=h2->GetNbinsX(); x++)
    {
	  n++;
	  h1->SetBinContent(n, h2->GetBinContent(x,y));
	  h1->SetBinError(n, h2->GetBinError(x,y));
	}
  }
  h1->Sumw2();
}


void unwrap2dhisto_3ch(TH2* h2, TH1* h1)
{
  int n=0;
  int bins_per_channel = h2->GetNbinsX()/3;
  int x_3ch = -99;

  for(int channel=0; channel<3; channel++)
	{
	  for(int y=1; y<=h2->GetNbinsY(); y++)
		{
		  for(int x=1; x<=bins_per_channel; x++)
			{
			  n++;
			  x_3ch = x + (channel*bins_per_channel);
			  h1->SetBinContent(n, h2->GetBinContent(x_3ch,y));
			  h1->SetBinError(n, h2->GetBinError(x_3ch,y));
			}
		}
	}

}


void rewrap1dhisto(TH1* h1, TH2* h2)
{
  int n=0;
  for(int y=1; y<=h2->GetNbinsY(); y++)
  {
	for(int x=1; x<=h2->GetNbinsX(); x++)
    {
	  n++;
	  h2->SetBinContent(x, y, h1->GetBinContent(n));
	  h2->SetBinError(x, y, h1->GetBinError(n));
	}
  }
  h1->Sumw2();

}

void rewrap1dhisto_3ch(TH1* h1, TH2* h2)
{
  int n=0;
  int bins_per_channel = h2->GetNbinsX()/3;
  int x_3ch = -99;

  for(int channel=0; channel<3; channel++)
	{
	  for(int y=1; y<=h2->GetNbinsY(); y++)
		{
		  for(int x=1; x<=bins_per_channel; x++)
			{
			  n++;
			  x_3ch = x + (channel*bins_per_channel);
			  h2->SetBinContent(x_3ch, y, h1->GetBinContent(n));
			  h2->SetBinError(x_3ch, y, h1->GetBinError(n));
			}
		}
	}

}

// The following is also copied from the KIT group code
double scaleBias = -99.9;
TUnfoldSys* myUnfold_TUnfoldGlobalPointerForTMinuit;
TH1D* myUnfold_hdataGlobalPointerForTMinuit;

void myUnfold_globalFunctionForMinuit(int &npar, double *gin, double &f, double *par, int iflag)
{
  const double logtau = par[0];
  const double myScaleBias = par[1];
  myUnfold_TUnfoldGlobalPointerForTMinuit->DoUnfold(pow(10, logtau), myUnfold_hdataGlobalPointerForTMinuit, myScaleBias);
  
  f = fabs(myUnfold_TUnfoldGlobalPointerForTMinuit->GetRhoAvg());
}



void minimizeRhoAverage(TUnfoldSys* unfold, TH1D* hdata, double log10min, double log10max)
{
  myUnfold_TUnfoldGlobalPointerForTMinuit = unfold;
  myUnfold_hdataGlobalPointerForTMinuit = hdata;
  
  // Instantiate Minuit for 2 parameters
  TMinuit minuit(2);
  minuit.SetFCN(myUnfold_globalFunctionForMinuit);
  minuit.SetPrintLevel(-1); // -1 no output, 1 output
  
  minuit.DefineParameter(0, "logtau", (log10min+log10max)/2, 1, log10min, log10max);
  minuit.DefineParameter(1, "scaleBias", scaleBias, 0, scaleBias, scaleBias);
  minuit.FixParameter(1);
  
  minuit.SetMaxIterations(100);
  minuit.Migrad();
  
  double bestlogtau = -1000;
  double bestlogtau_err = -1000; // error is meaningless because we don't have a likelihood, but method expects it
  minuit.GetParameter(0, bestlogtau, bestlogtau_err);
  unfold->DoUnfold(pow(10, bestlogtau), hdata, scaleBias); 
  
}
