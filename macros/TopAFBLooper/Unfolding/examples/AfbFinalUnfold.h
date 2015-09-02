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


Double_t stat_corr  [nbins1D]; //errors include syst error in the unfolding
Double_t stat_uncorr[nbins1D]; //errors do not include syst error in the unfolding
Double_t syst_corr  [nbinsx2Dalt];
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
syst_corr[0] = 0.002787;
syst_corr[1] = 0.002705;
syst_corr[2] = 0.003916;
syst_corr[3] = 0.004461;
syst_corr[4] = 0.005688;
syst_corr[5] = 0.004609;
syst_corr[6] = 0.006036;
syst_corr[7] = 0.004693;
syst_corr[8] = 0.003632;
syst_corr[9] = 0.003175;
syst_corr[10] = 0.002658;
syst_corr[11] = 0.002571;
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
syst_corr[0] = 0.005652;
syst_corr[1] = 0.005259;
syst_corr[2] = 0.005539;
syst_corr[3] = 0.006371;
syst_corr[4] = 0.003555;
syst_corr[5] = 0.002393;
syst_corr[6] = 0.002174;
syst_corr[7] = 0.002732;
syst_corr[8] = 0.004657;
syst_corr[9] = 0.005967;
syst_corr[10] = 0.006943;
syst_corr[11] = 0.006558;
      xmin1D=0.0;
      xmax1D=pi;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosTheta, bins_lepCosTheta + nbins1D + 1, xbins1D);
syst_corr[0] = 0.037501;
syst_corr[1] = 0.022551;
syst_corr[2] = 0.011134;
syst_corr[3] = 0.012406;
syst_corr[4] = 0.019139;
syst_corr[5] = 0.026255;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosTheta, bins_lepCosTheta + nbins1D + 1, xbins1D);
syst_corr[0] = 0.029792;
syst_corr[1] = 0.014941;
syst_corr[2] = 0.030211;
syst_corr[3] = 0.011674;
syst_corr[4] = 0.024992;
syst_corr[5] = 0.025815;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l}}#kern[-0.60]{*})";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosTheta, bins_lepCosTheta + nbins1D + 1, xbins1D);
syst_corr[0] = 0.028416;
syst_corr[1] = 0.012288;
syst_corr[2] = 0.013677;
syst_corr[3] = 0.007939;
syst_corr[4] = 0.020046;
syst_corr[5] = 0.024818;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)} cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.5; xbins1D[2]=-0.2; xbins1D[3]=0.0; xbins1D[4]=0.2; xbins1D[5]=0.5; xbins1D[6]=1.0;
      std::copy(bins_topSpinCorr, bins_topSpinCorr + nbins1D + 1, xbins1D);
syst_corr[0] = 0.007457;
syst_corr[1] = 0.024063;
syst_corr[2] = 0.035687;
syst_corr[3] = 0.028366;
syst_corr[4] = 0.034363;
syst_corr[5] = 0.011739;
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
syst_corr[0] = 0.002156;
syst_corr[1] = 0.011404;
syst_corr[2] = 0.017593;
syst_corr[3] = 0.014166;
syst_corr[4] = 0.010871;
syst_corr[5] = 0.002012;
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
syst_corr[0] = 0.017568;
syst_corr[1] = 0.022009;
syst_corr[2] = 0.009245;
syst_corr[3] = 0.020965;
syst_corr[4] = 0.009781;
syst_corr[5] = 0.016227;
      xmin1D=-1.0;
      xmax1D= 1.0;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos(#varphi)";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos(#varphi)}";
      //xbins1D[0]=-1.0; xbins1D[1]=-0.6; xbins1D[2]=-0.3; xbins1D[3]=0.0; xbins1D[4]=0.3; xbins1D[5]=0.6; xbins1D[6]=1.0;
      std::copy(bins_lepCosOpeningAngle, bins_lepCosOpeningAngle + nbins1D + 1, xbins1D);
syst_corr[0] = 0.009324;
syst_corr[1] = 0.008636;
syst_corr[2] = 0.009730;
syst_corr[3] = 0.023219;
syst_corr[4] = 0.015106;
syst_corr[5] = 0.040866;
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
syst_corr[0] = 0.003164;
syst_corr[1] = 0.003269;
syst_corr[2] = 0.001312;
syst_corr[3] = 0.001881;
syst_corr[4] = 0.002315;
syst_corr[5] = 0.002580;
syst_corr[6] = 0.002471;
syst_corr[7] = 0.002612;
syst_corr[8] = 0.001706;
syst_corr[9] = 0.000888;
syst_corr[10] = 0.003022;
syst_corr[11] = 0.003132;
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
syst_corr[0] = 0.008950;
syst_corr[1] = 0.014146;
syst_corr[2] = 0.008326;
syst_corr[3] = 0.007341;
syst_corr[4] = 0.010421;
syst_corr[5] = 0.005649;
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
syst_corr[0] = 0.003717;
syst_corr[1] = 0.004244;
syst_corr[2] = 0.015422;
syst_corr[3] = 0.011554;
syst_corr[4] = 0.003817;
syst_corr[5] = 0.003607;
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
syst_corr[0] = 0.013214;
syst_corr[1] = 0.024206;
syst_corr[2] = 0.013814;
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
syst_corr[0] = 0.009138;
syst_corr[1] = 0.014262;
syst_corr[2] = 0.024671;
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
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)}";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.022466;
syst_corr[1] = 0.031345;
syst_corr[2] = 0.037566;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.029994;
syst_corr[1] = 0.046919;
syst_corr[2] = 0.032904;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l}}#kern[-0.60]{*})";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.022363;
syst_corr[1] = 0.038842;
syst_corr[2] = 0.033216;
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
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)} cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      yaxislabel="M_{t#bar{t}}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      // xbins2D[0]=-1.0; xbins2D[1]=-0.5; xbins2D[2]=-0.2; xbins2D[3]=0.0; xbins2D[4]=0.2; xbins2D[5]=0.5; xbins2D[6]=1.0;
	  std::copy( bins_topSpinCorr, bins_topSpinCorr+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.028123;
syst_corr[1] = 0.031495;
syst_corr[2] = 0.039671;
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
syst_corr[0] = 0.022308;
syst_corr[1] = 0.021053;
syst_corr[2] = 0.008190;
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
syst_corr[0] = 0.031191;
syst_corr[1] = 0.016606;
syst_corr[2] = 0.020135;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos(#varphi)";
	  yaxislabel="M_{t#bar{t}}";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos(#varphi)}";
      std::copy( bins_lepCosOpeningAngle, bins_lepCosOpeningAngle+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.020869;
syst_corr[1] = 0.022364;
syst_corr[2] = 0.031806;
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
syst_corr[0] = 0.011692;
syst_corr[1] = 0.020405;
syst_corr[2] = 0.014942;
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
syst_corr[0] = 0.019642;
syst_corr[1] = 0.018600;
syst_corr[2] = 0.016936;
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
syst_corr[0] = 0.026669;
syst_corr[1] = 0.022259;
syst_corr[2] = 0.011321;
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
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepChargeAsym";
      asymlabel="A_{C}^{lep}";
      xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.8; xbins2D[6]=2.0;
      std::copy(bins_lepChargeAsym, bins_lepChargeAsym + nbinsx2Dalt + 1, xbins2Dalt);
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.008513;
syst_corr[1] = 0.016099;
syst_corr[2] = 0.009374;
      break;
    }
    //   Lepton Azimuthal Asymmetry 2
    case 1:
    {
      observablename="lep_azimuthal_asymmetry2";
      xaxislabel="|#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}|";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepAzimAsym2";
      asymlabel="A_{#Delta#phi}";
      //Double_t pi = 3.141592653589793;
      xbins2D[0]=0.0; xbins2D[1]=4.*pi/20.; xbins2D[2]=7.*pi/20.; xbins2D[3]=10.*pi/20.; xbins2D[4]=13.*pi/20.; xbins2D[5]=16.*pi/20.; xbins2D[6]=pi;
	  std::copy( bins_lepAzimAsym2, bins_lepAzimAsym2+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.014755;
syst_corr[1] = 0.014183;
syst_corr[2] = 0.017534;
theory_default[0] = 0.12572265215619549;
theory_default[1] = 0.17730999154723914;
theory_default[2] = 0.10262414675791232;
      break;
    }
    //   Top Polarization
    case 2:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)}";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.028864;
syst_corr[1] = 0.029459;
syst_corr[2] = 0.034268;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.026633;
syst_corr[1] = 0.041555;
syst_corr[2] = 0.058129;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l}}#kern[-0.60]{*})";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.027057;
syst_corr[1] = 0.032500;
syst_corr[2] = 0.039668;
      break;
    }
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)} cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.5; xbins2D[2]=-0.2; xbins2D[3]=0.0; xbins2D[4]=0.2; xbins2D[5]=0.5; xbins2D[6]=1.0;
	  std::copy( bins_topSpinCorr, bins_topSpinCorr+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.015347;
syst_corr[1] = 0.034877;
syst_corr[2] = 0.075586;
      break;
    }
    //   Top Asy III
    case 6:
    {
      observablename="top_rapiditydiff_Marco";
      //xaxislabel="|y_{top}|-|y_{tbar}|";
      xaxislabel="#Delta|y_{t}|";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="rapiditydiffMarco";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-0.7; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.7; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiffMarco, bins_rapiditydiffMarco+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.010535;
syst_corr[1] = 0.028918;
syst_corr[2] = 0.043742;
      break;
    }
    //   Top Charge Asymmetry
    case 7:
    {
      observablename="top_costheta_cms";
      xaxislabel="cos(#theta_{top})";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="topCosTheta";
      asymlabel="A_{FB}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.7; xbins2D[2]=-0.4; xbins2D[3]=0.0; xbins2D[4]=0.4; xbins2D[5]=0.7; xbins2D[6]=1.0;
	  std::copy( bins_topCosTheta, bins_topCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.011759;
syst_corr[1] = 0.039841;
syst_corr[2] = 0.029903;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos(#varphi)";
	  yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos(#varphi)}";
      std::copy( bins_lepCosOpeningAngle, bins_lepCosOpeningAngle+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.021592;
syst_corr[1] = 0.052984;
syst_corr[2] = 0.038834;
      break;
    }
    //   Lepton Azimuthal Asymmetry
    case 9:
    {
      observablename="lep_azimuthal_asymmetry";
      xaxislabel="#Delta#phi_{#font[12]{l#lower[-0.4]{+}l#lower[-0.48]{-}}}";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="lepAzimAsym";
      asymlabel="A_{#Delta#phi}";
      xbins2D[0]=-pi; xbins2D[1]=-0.67*pi; xbins2D[2]=-0.33*pi; xbins2D[3]=0.0; xbins2D[4]=0.33*pi; xbins2D[5]=0.67*pi; xbins2D[6]=pi;
	  std::copy( bins_lepAzimAsym, bins_lepAzimAsym+nbinsx2Dalt+1, xbins2Dalt );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2Dalt[0];
      xmax=xbins2Dalt[12];
syst_corr[0] = 0.007224;
syst_corr[1] = 0.016631;
syst_corr[2] = 0.010907;
      break;
    }
    //   Top Asy I
    case 10:
    {
      observablename="top_pseudorapiditydiff_cms";
      xaxislabel="|#eta_{t}|-|#eta_{#bar{t}}|";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="pseudorapiditydiff";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-1.0; xbins2D[2]=-0.5; xbins2D[3]=0.0; xbins2D[4]=0.5; xbins2D[5]=1.0; xbins2D[6]=2.0;
	  std::copy( bins_pseudorapiditydiff, bins_pseudorapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.009142;
syst_corr[1] = 0.022402;
syst_corr[2] = 0.033419;
      break;
    }
    //   Top Asy II
    case 11:
    {
      observablename="top_rapiditydiff_cms";
      xaxislabel="(y_{top}-y_{#bar{t}})(y_{top}+y_{#bar{t}})";
      yaxislabel="p_{T,t#bar{t}}";
      acceptanceName="rapiditydiff";
      asymlabel="A_{C}";
      //xbins2D[0]=-2.0; xbins2D[1]=-0.8; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.8; xbins2D[6]=2.0;
	  std::copy( bins_rapiditydiff, bins_rapiditydiff+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.008857;
syst_corr[1] = 0.025205;
syst_corr[2] = 0.048779;
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
syst_corr[0] = 0.006804;
syst_corr[1] = 0.009500;
syst_corr[2] = 0.005959;
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
syst_corr[0] = 0.010557;
syst_corr[1] = 0.011978;
syst_corr[2] = 0.016107;
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
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)}";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepPlusCosTheta";
      asymlabel="A_{P+}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.036935;
syst_corr[1] = 0.025054;
syst_corr[2] = 0.040610;
      break;
    }
      //   Top Polarization using negatively charged leptons
    case 3:
    {
      observablename="lepMinus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepMinusCosTheta";
      asymlabel="A_{P-}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.040557;
syst_corr[1] = 0.026845;
syst_corr[2] = 0.033840;
      break;
    }
      //   Top Polarization combining positively and negatively charged leptons
    case 4:
    {
      observablename="lepPlus_costheta_cms";
      xaxislabel="cos(^{}#theta_{#font[12]{l}}#kern[-0.60]{*})";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepCosTheta";
      asymlabel="A_{P}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.6; xbins2D[2]=-0.3; xbins2D[3]=0.0; xbins2D[4]=0.3; xbins2D[5]=0.6; xbins2D[6]=1.0;
	  std::copy( bins_lepCosTheta, bins_lepCosTheta+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.036169;
syst_corr[1] = 0.025865;
syst_corr[2] = 0.033269;
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
    //   Top Spin Correlation
    case 5:
    {
      observablename="top_spin_correlation";
      xaxislabel="cos(^{}#theta_{#font[12]{l#lower[-0.4]{+}}}#kern[-1.70]{*}#kern[+0.40]{)} cos(^{}#theta_{#font[12]{l#lower[-0.48]{-}}}#kern[-1.15]{*})";
      yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="topSpinCorr";
      asymlabel="A_{c1c2}";
      //xbins2D[0]=-1.0; xbins2D[1]=-0.5; xbins2D[2]=-0.2; xbins2D[3]=0.0; xbins2D[4]=0.2; xbins2D[5]=0.5; xbins2D[6]=1.0;
	  std::copy( bins_topSpinCorr, bins_topSpinCorr+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.024469;
syst_corr[1] = 0.019281;
syst_corr[2] = 0.027787;
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
syst_corr[0] = 0.026455;
syst_corr[1] = 0.020100;
syst_corr[2] = 0.021573;
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
syst_corr[0] = 0.019665;
syst_corr[1] = 0.025502;
syst_corr[2] = 0.022769;
      break;
    }
    //   Lepton opening angle
    case 8:
    {
      observablename="lep_cos_opening_angle";
      xaxislabel="cos(#varphi)";
	  yaxislabel="|y_{t#bar{t}}|";
      acceptanceName="lepCosOpeningAngle";
      asymlabel="A_{cos(#varphi)}";
      std::copy( bins_lepCosOpeningAngle, bins_lepCosOpeningAngle+nbinsx2D+1, xbins2D );
      ymin=ybins2D[0];
      ymax=ybins2D[3];
      xmin=xbins2D[0];
      xmax=xbins2D[6];
syst_corr[0] = 0.050125;
syst_corr[1] = 0.015460;
syst_corr[2] = 0.020665;
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
syst_corr[0] = 0.005535;
syst_corr[1] = 0.003738;
syst_corr[2] = 0.006707;
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
syst_corr[0] = 0.024531;
syst_corr[1] = 0.019521;
syst_corr[2] = 0.014407;
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
syst_corr[0] = 0.024684;
syst_corr[1] = 0.015457;
syst_corr[2] = 0.027697;
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
