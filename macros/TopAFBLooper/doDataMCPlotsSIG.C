//////////////////////////////////////////////////////////////////////////////
// Assumptions:
//
// - files are in output directory, and their names have following pattern:
//     "wmn"
//     "J0m"
//     "J1m"
//     "J2m"
//     "J3m"
//     "J4m" + user-provided string
//     "J5m"
//     "bb"
//     "cc"
//     "wtn"
//     "zmm"
//     "ztt"
//     "tt"
// - Canvases will be saved in plots directory
// - Possible interventions:
//   - change stack order: search for "stackingOrder" vector
//   - data vs lumi normalization: "Renormalization of data vs MC"
//   - stack Y axis max: "Special requests for Y axis range!"
//   - definition of Poisson errors (!): GetPoissonizedGraph
//
// Usage:
//
// - compile the macro in root:
//  [] .L doDataMCPlots.C++
// - launch functions for 1D and 2D plotting:
//  [] doDataMCPlots("_staco","myDataFile.root",310.,0.611,1,0)
//
//////////////////////////////////////////////////////////////////////////////
#include <TH1F.h>
#include <TH2F.h>
#include <THStack.h>
#include <TLegend.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include <TLine.h>
#include <TList.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <sstream>
#include <iomanip>

using namespace std;

// This is meant to be passed as the third argument, the predicate, of the standard library sort algorithm
inline bool sortByIntegral(TH1F *h1, TH1F *h2 )
{
    return h1->Integral() < h2->Integral();
}

double error_poisson_up(double data);
double error_poisson_down(double data);

double compatibilityTest(TH1F *hist, THStack *hstack);

double straightline ( double *x, double *parm)
{
    return parm[0] + parm[1] * x[0];
}

float addSqr( float e1 , float e2 )
{
    float err2 = pow(e1, 2) + pow(e2, 2);
    return sqrt(err2);
}

void zeroHist( TH1F *&h )
{

    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
        h->SetBinContent(i, 0.);
        h->SetBinError(i, 0.);
    }
}

TGraphAsymmErrors *GetPoissonizedGraph(TH1F *histo);

void myText(Double_t x, Double_t y, Color_t color, char *text)
{

    TLatex l;
    l.SetNDC();
    l.SetTextColor(color);
    l.SetTextSize(0.05);
    l.DrawLatex(x, y, text);
}

void myOtherText(Double_t x, Double_t y, double tsize, char *text)
{
    TLatex l;
    l.SetTextSize(tsize);
    l.SetNDC();
    l.DrawLatex(x, y, text);
}

void myCMSLabel(Double_t x, Double_t y, int type)
{

    // CMS Label types:
    // 0: no label
    // 1: CMS

    if (type == 0)
        return;

    TLatex l;
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextSize(0.04);
    l.DrawLatex(x, y, "CMS");
    l.SetTextFont(42);
    l.SetTextSize(0.04);
    l.DrawLatex(x + 0.13, y, "Preliminary");

}

void printYields( vector<TH1F *> mc , const char *labels[] , TH1F *chdata , bool latex = false );
void initSymbols(bool);
void printLine(bool);
void printHeader();
void GetAfberr(TH1F* h, Float_t &afb, Float_t  &afberr);
float GetAfb(TH1F* h);
float GetAfbLHalf(TH1F* h);
float GetAfbRHalf(TH1F* h);
void print( TH1F *h , string label , bool correlatedError = false );
const char *pm;
const char *percent;
const char *delim;
const char *delimstart;
const char *delimend;
const char *ee;
const char *mm;
const char *em;
const char *e;
const char *m;
int   width1;
int   width2;
int   extrawidth;
int   linelength;

float histError( TH1F *h , int minbin , int maxbin )
{

    float err2 = 0;
    for ( int i = minbin ; i <= maxbin ; ++i )
    {
        err2 += pow( h->GetBinError(i) , 2 );
    }

    return sqrt(err2);
}

void scaleHistError( TH1F *&h, float sf, float esf )
{

    float ierr = 0.;
    for (int i = 1; i <= h->GetNbinsX(); ++i)
    {
        ierr = h->GetBinContent(i) > 0. ? h->GetBinError(i) / h->GetBinContent(i) : 0.;
        ierr = addSqr( ierr, esf / sf ) * h->GetBinContent(i);
        h->SetBinError(i, ierr);
    }

}

pair<float, float> datamcsf( TH1F *h_dt, vector<TH1F *> h_mc, int minbin, int maxbin )
{

    float n_dt = h_dt->Integral(minbin, maxbin + 1);
    float e_dt = histError( h_dt, minbin, maxbin + 1);
    //Add up MC samples
    if (h_mc.size() < 1) return make_pair(0., 0.);

    int nbins = h_mc.at(0)->GetNbinsX();
    float xmin = h_mc.at(0)->GetBinLowEdge(1);
    float xmax = h_mc.at(0)->GetBinLowEdge(nbins);
    cout << "Going from " << minbin << " to " << maxbin
         << " nbin range " << h_mc.at(0)->GetXaxis()->GetBinLowEdge(minbin)
         << " - " << h_mc.at(0)->GetXaxis()->GetBinUpEdge(maxbin) << endl;
    TH1F *hmctot = new TH1F("hmctot", "hmctot", nbins, xmin, xmax);
    hmctot->Sumw2();

    for (unsigned int imc = 0 ; imc < h_mc.size() ; imc++)
    {
        h_mc.at(imc)->Sumw2();
        if ( imc == 0 ) hmctot = (TH1F *) h_mc.at(imc)->Clone();
        else           hmctot->Add(h_mc.at(imc));
    }

    float n_mc = hmctot->Integral(minbin, maxbin);
    float e_mc = histError( hmctot, minbin, maxbin );

    float sf_dtmc = n_dt / n_mc;
    float e_sf_dtmc = sqrt(pow( (e_dt / n_dt), 2) + pow( (e_mc / n_mc), 2 )) * sf_dtmc;

    return make_pair(sf_dtmc, e_sf_dtmc);
}


//////////////////////////////////////////////////////////////////////////////
void doDataMCPlotsSIG(const char *region = "SIG", const char *ttbar_tag = "mcatnlo")
{

    TString ttbar_string = ttbar_tag;

    //derive scale factors

    //list of samples
    //const int MCID = 10;
    const int MCID = 10;

    const int nRegions = 15;

    enum regions {SIG=0,
        CR0,
        CR1,
        CR1v,
        CR2,
        CR2v,
        CR3,
        CR3v,
        CR4,
        CR4v,
        CR40,
        CR5,
        CR5v,
        CR6,
        CR6v
    };


    const char *histtag[nRegions] = 
    {
        "h_sig_",
        "h_cr0_",
        "h_cr1_",
        "h_cr1v_",
        "h_cr2_",
        "h_cr2v_",
        "h_cr3_",
        "h_cr3v_",
        "h_cr4_",
        "h_cr4v_",
        "h_cr40_",
        "h_cr5_",
        "h_cr5v_",
        "h_cr6_",
        "h_cr6v_"
    };
    int histtag_number = SIG;
    if(strncmp(region,"CR0",1000) == 0) histtag_number = CR0;
    if(strncmp(region,"CR1",1000) == 0) histtag_number = CR1;
    if(strncmp(region,"CR1v",1000) == 0) histtag_number = CR1v;
    if(strncmp(region,"CR2",1000) == 0) histtag_number = CR2;
    if(strncmp(region,"CR2v",1000) == 0) histtag_number = CR2v;
    if(strncmp(region,"CR3",1000) == 0) histtag_number = CR3;
    if(strncmp(region,"CR3v",1000) == 0) histtag_number = CR3v;
    if(strncmp(region,"CR4",1000) == 0) histtag_number = CR4;
    if(strncmp(region,"CR4v",1000) == 0) histtag_number = CR4v;
    if(strncmp(region,"CR40",1000) == 0) histtag_number = CR40;
    if(strncmp(region,"CR5",1000) == 0) histtag_number = CR5;
    if(strncmp(region,"CR5v",1000) == 0) histtag_number = CR5v;
    if(strncmp(region,"CR6",1000) == 0) histtag_number = CR6;
    if(strncmp(region,"CR6v",1000) == 0) histtag_number = CR6v;


    const char *dirtag[nRegions] = 
    {
        "SIG",
        "CR0",
        "CR1",
        "CR1v",
        "CR2",
        "CR2v",
        "CR3",
        "CR3v",
        "CR4",
        "CR4v",
        "CR40",
        "CR5",
        "CR5v",
        "CR6",
        "CR6v"
    };

    //Be careful changing the order of these (some hard coding may still exist). At least always keep ttdl first.
    const char *mcsample[MCID] =
    {
        "ttdl",
        "ttsl",
        //"ttfake",
        "w1to4jets",
        "tW_lepsl",
        "tW_lepdl",
        "diboson",
        //"DY1to4Jtot",
        "DY1to4Jeemm",
        "DY1to4Jtautau",
        "ttV",
        "triboson",
        //"ttfake_mcatnlo",
        //"tttt_mcatnlo"
    };

    enum sample {TTDL = 0,
                 TTSL,
                 //TTFA,
                 WJETS,
                 TWSL,
                 TWDL,
                 DIBO,
                 DY,
                 DYTT,
                 TTV,
                 VVV
                };


    const int CRtargetsamplenumber[nRegions] = 
    {
        TTDL,
        DY,
        DY,
        DY,
        DY,
        DY,
        DY,
        DY,
        TWDL,
        DYTT,
        DYTT,
        TTSL,
        TTSL,
        TTSL,
        TTSL
    };
/*
    //apply SFs to other backgrounds when deriving the SFs for this CR?
    const int SFsApply[nRegions] = 
    {
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        1,
        0,  //not clear whether the fake and DY SFs are applicable in this region
        0,  //not clear whether the fake and DY SFs are applicable in this region
        1,
        1,
        0,  //no MET cut means we don't want to scale the DY
        0   //no MET cut means we don't want to scale the DY
    };
*/

    //0: only fakes
    //1: full SFs
    //2: full except DY b 
    //3: only DY b and fakes
    const int SFsApply[nRegions] = 
    {
        1, //SIG
        2, //CR0
        1, //CR1
        2, //CR1v
        3, //CR2
        0, //CR2v
        3, //CR3
        0, //CR3v
        1, //CR4
        2, //CR4v
        0, //CR40
        1, //CR5
        2, //CR5v
        3, //CR6
        0  //CR6v
    };


    const char *legend[MCID] =
    {
        "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}}",
        //"\\ttdl\\",
        "t#bar{t} #rightarrow #font[12]{l^{#pm}} + jets",
        //"\\ttsl\\",
        //"t#bar{t} #rightarrow all jets",
        "W+jets",
        // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (e/#mu)",
        // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (lost)",
        // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{had}#rightarrow1-prong)",
        // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{had}#rightarrow3-prong)",
        // "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}} (#tau_{lep})",
        "single top (s/t-chan, 1-lep)",
        //"rare",
        "single top (tW, 2-lep)",
        //    "single top",
        "WW/WZ/ZZ",
        //"DY+jets",
        "DY #rightarrow ee/#mu#mu+jets",
        "DY #rightarrow #tau#tau+jets",
        "ttW/Z/#gamma",
        "triboson",
        //"t#bar{t} #rightarrow all jets",
        //"t#bar{t}t#bar{t}"
    };
    //    "triboson"};
    //    "QCD"};

    //  const int mccolor[]={7,2,6,4,5,kOrange,9,kAzure-9, 8, kViolet,kGreen+1, 15,12,13,27};
    //  const int mccolor[]={7,2,6,4,kOrange,8,9,kAzure-9, 5, kViolet,kGreen+1, 15,12,13,27};
    //const int mccolor[] = {7, 2, kRed-2, 6, 4, kOrange, 8, kViolet + 1, kAzure - 9, kGreen - 2, kViolet, kGreen + 1, 15, 12, 13, 27};
    const int mccolor[] = {7, 2, 6, 4, kOrange, 8, kViolet + 1, kViolet - 7, kAzure - 9, kGreen - 2, kViolet, kGreen + 1, 15, 12, 13, 27};
    //  const int mccolor[]={7,2,6,4,kOrange,8,5,kAzure-9, 9,5,kViolet,kGreen+1, 15,12,13,27};

    //-------------------------------
    // SINGLE LEPTON - MT SCALING
    //-------------------------------

    const int nCh = 3;
    const char *leptag[nCh] = {"_diel", "_dimu", "_mueg"};
    const char *leplabel[nCh] = {"ee", "#mu#mu", "e#mu"};
    enum channel {DIEL = 0, DIMU, MUEG };
    const char *histoname = "h_sig_m_top";

    bool doverbose = true;
    bool scalettdiltodata = true;
    bool scalebkgtoCRs = true;
    bool scaletoyieldsafterttbarsol = false;
    bool printnumbersonplots = true;

    double lumi = 19.5;

    double KSsum = 0.;
    double chi2sum = 0.;

    double KSsum_channels[nCh] = {0.};
    double chi2sum_channels[nCh] = {0.};

    int isr = 0;
    const char *metcut[1] = {""};
    const char *channelhist[2] = {"channel", "channel_withttbarsol"};

    //calculate background SFs from CRs
    double bkgsf[4][MCID][nCh];
    double bkgsferr[4][MCID][nCh];
    float CRsf[nRegions][nCh];
    float CRsferr[nRegions][nCh];
    float CRsftot[nRegions];
    float CRsftoterr[nRegions];
    if(scalebkgtoCRs || (scalettdiltodata && histtag_number != SIG) ) {

        TFile *dt_dl_temp[nCh][nRegions];
        TFile *mc_dl_temp[MCID][nRegions];

        for (int i = 0; i < nRegions; ++i)
        {   

            dt_dl_temp[DIEL][i] = TFile::Open(Form("%soutput/data_diel_histos.root",dirtag[i]));
            dt_dl_temp[DIMU][i] = TFile::Open(Form("%soutput/data_dimu_histos.root",dirtag[i]));
            dt_dl_temp[MUEG][i] = TFile::Open(Form("%soutput/data_mueg_histos.root",dirtag[i]));

            for (int j = 0; j < MCID; ++j)
            {
                if (j < 2)
                    mc_dl_temp[j][i] = TFile::Open(Form("%soutput/%s_%s_histos.root", dirtag[i],
                                                mcsample[j], ( j != TTDL && ttbar_string.Contains("FullLept") ? "mcatnlo" : ttbar_tag )));
                else
                    mc_dl_temp[j][i] = TFile::Open(Form("%soutput/%s_histos.root", dirtag[i],
                                                mcsample[j]));
            }
        }

        for (int j = 0; j < MCID; ++j)
        {
            for (int leptype = 0; leptype < nCh; ++leptype)
            {
                for (int sftype = 0; sftype < 4; ++sftype)
                {
                    bkgsf[sftype][j][leptype] = 1.;
                    bkgsferr[sftype][j][leptype] = 0.;
                }
            }
        }

        int iteration = 0;
        //6 iterations is plenty
        while (iteration < 6) {
            iteration++;

            for (int i = 0; i < nRegions; ++i)
            {   
                if(i==CR40) continue;
                CRsftot[i] = 0;
                CRsftoterr[i] = 0;

                for (int leptype = 0; leptype < nCh; ++leptype)
                {
                    TH1F *h_dt1d;
                    TH1F *h_mc1d[MCID];
                    TH1F *h_mc1d_tot;

                    h_dt1d = (TH1F *)dt_dl_temp[leptype][i]->Get(Form("%s%s%s%s", histtag[i], channelhist[scaletoyieldsafterttbarsol], "", leptag[leptype]));
                    h_dt1d->SetName(Form("%s%s", histtag[i], channelhist[scaletoyieldsafterttbarsol]));

                    bool doinit = false;
                    for (int j = 0; j < MCID; ++j)
                    {

                        h_mc1d[j] = (TH1F *)mc_dl_temp[j][i]->Get(Form("%s%s%s%s", histtag[i], channelhist[scaletoyieldsafterttbarsol], "", leptag[leptype]));


                        if (h_mc1d[j] == 0)
                        {
                            h_mc1d[j] = (TH1F *)dt_dl_temp[leptype][i]->Get(Form("%s%s%s%s", histtag[i], channelhist[scaletoyieldsafterttbarsol], "", leptag[leptype]));
                            h_mc1d[j]->SetName(Form("%s_%s%s", mcsample[j], histtag[i], channelhist[scaletoyieldsafterttbarsol]));
                            zeroHist(h_mc1d[j]);
                        }
                        h_mc1d[j]->SetName(Form("%s_%s%s", mcsample[j], histtag[i], channelhist[scaletoyieldsafterttbarsol]));

                        if(j == TTDL && j != CRtargetsamplenumber[i]) h_mc1d[j]->Scale(bkgsf[1][j][leptype]); //use previously calculated SF for TTDL
                        if(  j != TTDL && j != CRtargetsamplenumber[i] && ( (CRtargetsamplenumber[i] != TTSL ) || (CRtargetsamplenumber[i] == TTSL && j != TTSL && j != WJETS && j != TWSL )  )  ) h_mc1d[j]->Scale(bkgsf[SFsApply[i]][j][leptype]); //use previously calculated SFs for the other backgrounds, if SFsApply[region]
                        //if(  j != TTDL && j != CRtargetsamplenumber[i] && ( (CRtargetsamplenumber[i] != TTSL ) || (CRtargetsamplenumber[i] == TTSL && j != TTSL && j != WJETS && j != TWSL && j != TTFA )  )  ) h_mc1d[j]->Scale(bkgsf[SFsApply[i]][j][leptype]); //use previously calculated SFs for the other backgrounds, if SFsApply[region]

                        if (!doinit)
                        {
                            h_mc1d_tot = (TH1F *)h_mc1d[j]->Clone(Form("mctot_%s%s", histtag[i], channelhist[scaletoyieldsafterttbarsol]));
                            doinit = true;
                        }
                        else h_mc1d_tot->Add(h_mc1d[j]);

                    }


                    double mcallerr = 0.;
                    double dtallerr = 0.;
                    double mcCRsampleerrtemp = 0.;
                    double mcCRsampleerr = 0.;

                    float mcall = h_mc1d_tot->IntegralAndError(0, h_mc1d_tot->GetNbinsX()+1, mcallerr );
                    //add the uncertainties from the other SFs:
                    double mcallerrsftotalsquared = 0.;
                    float mcallerrsf[MCID];
                    for (int j = 0; j < MCID; ++j)
                    {
                        if(j == CRtargetsamplenumber[i]) continue; //only want to add the uncertainty from the other SFs 
                        //if( j==WJETS || j==TWSL || j==TTFA ) continue; //these uncertainties are correlated and these components are always scaled together, so sum them linearly before combining 
                        if( j==WJETS || j==TWSL ) continue; //these uncertainties are correlated and these components are always scaled together, so sum them linearly before combining (below)
                        if( j==DYTT && CRtargetsamplenumber[i] != DY && CRtargetsamplenumber[i] != DYTT ) continue; //DYTT is correlated with DY, but these components are not scaled together. For non-DY CRs, continue and sum the uncertainties linearly.

                        mcallerrsf[j] = bkgsferr[SFsApply[i]][j][leptype] * h_mc1d[j]->Integral();

                        if( j==TTSL ) mcallerrsf[j] += bkgsferr[SFsApply[i]][WJETS][leptype] * h_mc1d[WJETS]->Integral();
                        if( j==TTSL ) mcallerrsf[j] += bkgsferr[SFsApply[i]][TWSL][leptype] * h_mc1d[TWSL]->Integral();
                        //if( j==TTSL ) mcallerrsf[j] += bkgsferr[SFsApply[i]][TTFA][leptype] * h_mc1d[TTFA]->Integral();

                        //for non-DY CRs, add the DYTT and DY uncertainties assuming 100% correlation (this is slightly conservative because the true correlation is <100%)
                        if( j==DY && CRtargetsamplenumber[i] != DYTT && CRtargetsamplenumber[i] != DY ) mcallerrsf[j] += bkgsferr[SFsApply[i]][DYTT][leptype] * h_mc1d[DYTT]->Integral();

                        mcallerrsftotalsquared += mcallerrsf[j]*mcallerrsf[j] ; 
                    }

                    //cout<<"mcallerr: "<<mcallerr<<" additional error from SFs: "<<sqrt(mcallerrsftotalsquared)<<endl;
                    mcallerr = sqrt(mcallerr*mcallerr+mcallerrsftotalsquared);

                    float dtall = h_dt1d->IntegralAndError(0, h_dt1d->GetNbinsX()+1, dtallerr );
                    float mcCRsample = h_mc1d[CRtargetsamplenumber[i]]->IntegralAndError(0, h_mc1d[CRtargetsamplenumber[i]]->GetNbinsX()+1, mcCRsampleerrtemp );
                    mcCRsampleerr += mcCRsampleerrtemp*mcCRsampleerrtemp;
                    if(CRtargetsamplenumber[i]==TTSL) {
                        mcCRsample += h_mc1d[WJETS]->IntegralAndError(0, h_mc1d[WJETS]->GetNbinsX()+1, mcCRsampleerrtemp );  //add w+jets to ttsl
                        mcCRsampleerr += mcCRsampleerrtemp*mcCRsampleerrtemp;
                        mcCRsample += h_mc1d[TWSL]->IntegralAndError(0, h_mc1d[TWSL]->GetNbinsX()+1, mcCRsampleerrtemp );  //add 1l single top to ttsl
                        mcCRsampleerr += mcCRsampleerrtemp*mcCRsampleerrtemp;
                        //mcCRsample += h_mc1d[TTFA]->IntegralAndError(0, h_mc1d[TTFA]->GetNbinsX()+1, mcCRsampleerrtemp );  //add ttfake to ttsl
                        //mcCRsampleerr += mcCRsampleerrtemp*mcCRsampleerrtemp;
                    }
                    mcCRsampleerr = sqrt(mcCRsampleerr);
                    double CRsftemp = 1. + (dtall-mcall)/mcCRsample;
                    double CRsftemperr = sqrt( pow(mcCRsampleerr*(mcall-mcCRsample-dtall)/mcCRsample/mcCRsample , 2) + (dtallerr*dtallerr + mcallerr*mcallerr - mcCRsampleerr*mcCRsampleerr)/mcCRsample/mcCRsample );
                    //cout<<sqrt( pow(mcCRsampleerr*(mcall-mcCRsample-dtall)/mcCRsample/mcCRsample , 2) )<<" "<<sqrt((mcallerr*mcallerr - mcCRsampleerr*mcCRsampleerr)/mcCRsample/mcCRsample)<<" "<<sqrt((dtallerr*dtallerr)/mcCRsample/mcCRsample)<<endl;
                    //cout<<"Uncertainty increase factor from full error propagation: "<< CRsftemperr / (dtallerr/mcCRsample) <<endl;
                    CRsf[i][leptype] = CRsftemp;
                    CRsferr[i][leptype] = CRsftemperr;
                    //cout<<"SF from "<<dirtag[i]<<" for "<< (CRtargetsamplenumber[i] == TTSL ? "fakes" : legend[CRtargetsamplenumber[i]]) <<" for "<<leplabel[leptype]<<": "<<CRsftemp<<" +/- "<<CRsftemperr<<endl;

                    //immediately fill ttdl SFs so they can be used for the CRs later in the loop
                    if(i == SIG) {
                        bkgsf[1][TTDL][leptype] = CRsf[i][leptype];
                        bkgsferr[1][TTDL][leptype] = CRsferr[i][leptype];
                    }

                    //calculate weighted average
                    if( (CRtargetsamplenumber[i] == DY && leptype < MUEG) || (CRtargetsamplenumber[i] == DYTT && leptype == MUEG) || (CRtargetsamplenumber[i] == TTSL && leptype == DIMU) || ( CRtargetsamplenumber[i] != DY && CRtargetsamplenumber[i] != DYTT && CRtargetsamplenumber[i] != TTSL) ) {
                        CRsftot[i] += CRsf[i][leptype]/CRsferr[i][leptype]/CRsferr[i][leptype];
                        CRsftoterr[i] += 1./CRsferr[i][leptype]/CRsferr[i][leptype];
                    }

                }

                CRsftot[i] /= CRsftoterr[i];
                CRsftoterr[i] = 1./sqrt(CRsftoterr[i]);
            }

            cout<<bkgsf[1][TTDL][0]<<" ± "<<bkgsferr[1][TTDL][0]<<" "<<bkgsf[1][TTDL][1]<<" ± "<<bkgsferr[1][TTDL][1]<<" "<<bkgsf[1][TTDL][2]<<" ± "<<bkgsferr[1][TTDL][2]<<" "<<bkgsf[1][TTSL][0]<<" ± "<<bkgsferr[1][TTSL][0]<<" "<<bkgsf[1][DY][0]<<" ± "<<bkgsferr[1][DY][0]<<" "<<bkgsf[1][DYTT][0]<<" ± "<<bkgsferr[1][DYTT][0]<<" "<<endl;


            for (int leptype = 0; leptype < nCh; ++leptype)
            {

                     bkgsf[1][TTSL][leptype] = CRsftot[CR5];
                     bkgsferr[1][TTSL][leptype] = CRsftoterr[CR5];
                     //bkgsf[1][TTFA][leptype] = CRsftot[CR5];
                     //bkgsferr[1][TTFA][leptype] = CRsftoterr[CR5];
                     bkgsf[1][WJETS][leptype] = CRsftot[CR5];
                     bkgsferr[1][WJETS][leptype] = CRsftoterr[CR5];
                     bkgsf[1][TWSL][leptype] = CRsftot[CR5];
                     bkgsferr[1][TWSL][leptype] = CRsftoterr[CR5];
                     bkgsf[1][TWDL][leptype] = 1.;
                     bkgsferr[1][TWDL][leptype] = 0.;
                     bkgsf[1][DIBO][leptype] = 1.;
                     bkgsferr[1][DIBO][leptype] = 0.;
                     bkgsf[1][DY][leptype] = CRsftot[CR1];
                     bkgsferr[1][DY][leptype] = CRsftoterr[CR1];
                     bkgsf[1][DYTT][leptype] = CRsftot[CR2];
                     bkgsferr[1][DYTT][leptype] = CRsftoterr[CR2];
                     bkgsf[1][TTV][leptype] = 1.;
                     bkgsferr[1][TTV][leptype] = 0.;
                     bkgsf[1][VVV][leptype] = 1.;
                     bkgsferr[1][VVV][leptype] = 0.;

                     bkgsf[2][TTSL][leptype] = CRsftot[CR5];
                     bkgsferr[2][TTSL][leptype] = CRsftoterr[CR5];
                     //bkgsf[2][TTFA][leptype] = CRsftot[CR5];
                     //bkgsferr[2][TTFA][leptype] = CRsftoterr[CR5];
                     bkgsf[2][WJETS][leptype] = CRsftot[CR5];
                     bkgsferr[2][WJETS][leptype] = CRsftoterr[CR5];
                     bkgsf[2][TWSL][leptype] = CRsftot[CR5];
                     bkgsferr[2][TWSL][leptype] = CRsftoterr[CR5];
                     bkgsf[2][TWDL][leptype] = 1.;
                     bkgsferr[2][TWDL][leptype] = 0.;
                     bkgsf[2][DIBO][leptype] = 1.;
                     bkgsferr[2][DIBO][leptype] = 0.;
                     bkgsf[2][DY][leptype] = CRsftot[CR0];
                     bkgsferr[2][DY][leptype] = CRsftoterr[CR0];
                     bkgsf[2][DYTT][leptype] = 1.;
                     bkgsferr[2][DYTT][leptype] = 1.;
                     bkgsf[2][TTV][leptype] = 1.;
                     bkgsferr[2][TTV][leptype] = 0.;
                     bkgsf[2][VVV][leptype] = 1.;
                     bkgsferr[2][VVV][leptype] = 0.;

                     bkgsf[3][TTSL][leptype] = CRsftot[CR5];
                     bkgsferr[3][TTSL][leptype] = CRsftoterr[CR5];
                     //bkgsf[3][TTFA][leptype] = CRsftot[CR5];
                     //bkgsferr[3][TTFA][leptype] = CRsftoterr[CR5];
                     bkgsf[3][WJETS][leptype] = CRsftot[CR5];
                     bkgsferr[3][WJETS][leptype] = CRsftoterr[CR5];
                     bkgsf[3][TWSL][leptype] = CRsftot[CR5];
                     bkgsferr[3][TWSL][leptype] = CRsftoterr[CR5];
                     bkgsf[3][TWDL][leptype] = 1.;
                     bkgsferr[3][TWDL][leptype] = 0.;
                     bkgsf[3][DIBO][leptype] = 1.;
                     bkgsferr[3][DIBO][leptype] = 0.;
                     bkgsf[3][DY][leptype] = CRsftot[CR2];
                     bkgsferr[3][DY][leptype] = CRsftoterr[CR2];
                     bkgsf[3][DYTT][leptype] = CRsftot[CR2];
                     bkgsferr[3][DYTT][leptype] = CRsftoterr[CR2];
                     bkgsf[3][TTV][leptype] = 1.;
                     bkgsferr[3][TTV][leptype] = 0.;
                     bkgsf[3][VVV][leptype] = 1.;
                     bkgsferr[3][VVV][leptype] = 0.;

                     bkgsf[0][TTSL][leptype] = CRsftot[CR5];
                     bkgsferr[0][TTSL][leptype] = CRsftoterr[CR5];
                     //bkgsf[0][TTFA][leptype] = CRsftot[CR5];
                     //bkgsferr[0][TTFA][leptype] = CRsftoterr[CR5];
                     bkgsf[0][WJETS][leptype] = CRsftot[CR5];
                     bkgsferr[0][WJETS][leptype] = CRsftoterr[CR5];
                     bkgsf[0][TWSL][leptype] = CRsftot[CR5];
                     bkgsferr[0][TWSL][leptype] = CRsftoterr[CR5];
                     bkgsf[0][TWDL][leptype] = 1.;
                     bkgsferr[0][TWDL][leptype] = 0.;
                     bkgsf[0][DIBO][leptype] = 1.;
                     bkgsferr[0][DIBO][leptype] = 0.;
                     bkgsf[0][DY][leptype] = 1.;
                     bkgsferr[0][DY][leptype] = 0.;
                     bkgsf[0][DYTT][leptype] = 1.;
                     bkgsferr[0][DYTT][leptype] = 0.;
                     bkgsf[0][TTV][leptype] = 1.;
                     bkgsferr[0][TTV][leptype] = 0.;
                     bkgsf[0][VVV][leptype] = 1.;
                     bkgsferr[0][VVV][leptype] = 0.;

                     bkgsf[0][TTDL][leptype] = bkgsf[1][TTDL][leptype];
                     bkgsferr[0][TTDL][leptype] = bkgsferr[1][TTDL][leptype];
                     bkgsf[2][TTDL][leptype] = bkgsf[1][TTDL][leptype];
                     bkgsferr[2][TTDL][leptype] = bkgsferr[1][TTDL][leptype];
                     bkgsf[3][TTDL][leptype] = bkgsf[1][TTDL][leptype];
                     bkgsferr[3][TTDL][leptype] = bkgsferr[1][TTDL][leptype];
            }

        }

        for (int i = 0; i < nRegions; ++i)
        {   

            dt_dl_temp[DIEL][i] -> Close();
            dt_dl_temp[DIMU][i] -> Close();
            dt_dl_temp[MUEG][i] -> Close();

            for (int j = 0; j < MCID; ++j)
            {
                mc_dl_temp[j][i] -> Close();
            }
        }

    }

    if(scalebkgtoCRs || (scalettdiltodata && histtag_number != SIG) ) cout<<"Average SR sf for ttdl: "<<CRsftot[SIG]<<" +/- "<<CRsftoterr[SIG]<<endl;
    if(scalebkgtoCRs) cout<<"average CR5 sf for fakes: "<<CRsftot[CR5]<<" +/- "<<CRsftoterr[CR5]<<endl;
    if(scalebkgtoCRs) cout<<"average CR5v sf for fakes: "<<CRsftot[CR5v]<<" +/- "<<CRsftoterr[CR5v]<<endl;
    if(scalebkgtoCRs) cout<<"average CR6 sf for fakes: "<<CRsftot[CR6]<<" +/- "<<CRsftoterr[CR6]<<endl;
    if(scalebkgtoCRs) cout<<"average CR6v sf for fakes: "<<CRsftot[CR6v]<<" +/- "<<CRsftoterr[CR6v]<<endl;
    if(scalebkgtoCRs) cout<<"average CR1 sf for DYeemm:  "<<CRsftot[CR1]<<" +/- "<<CRsftoterr[CR1]<<endl;
    if(scalebkgtoCRs) cout<<"CR0 sf for DYeemm (exc. b modelling SF):  "<<CRsftot[CR0]<<" +/- "<<CRsftoterr[CR0]<<endl;
    if(scalebkgtoCRs) cout<<"composite sf for DYeemm (CR0*CR2/CR2v):  "<<CRsftot[CR0]*CRsftot[CR2]/CRsftot[CR2v]<<" +/- "<<(CRsftot[CR0]*CRsftot[CR2]/CRsftot[CR2v])*sqrt(pow(CRsftoterr[CR0]/CRsftot[CR0],2)+pow(CRsftoterr[CR2]/CRsftot[CR2],2)+pow(CRsftoterr[CR2v]/CRsftot[CR2v],2))<<endl;
    if(scalebkgtoCRs) cout<<"composite sf for DYeemm (CR1v*CR2/CR2v):  "<<CRsftot[CR1v]*CRsftot[CR2]/CRsftot[CR2v]<<" +/- "<<(CRsftot[CR1v]*CRsftot[CR2]/CRsftot[CR2v])*sqrt(pow(CRsftoterr[CR1v]/CRsftot[CR1v],2)+pow(CRsftoterr[CR2]/CRsftot[CR2],2)+pow(CRsftoterr[CR2v]/CRsftot[CR2v],2))<<endl;
    if(scalebkgtoCRs) cout<<"average CR3 sf:  "<<CRsftot[CR3]<<" +/- "<<CRsftoterr[CR3]<<endl;
    if(scalebkgtoCRs) cout<<"average CR3v sf:  "<<CRsftot[CR3v]<<" +/- "<<CRsftoterr[CR3v]<<endl;
    if(scalebkgtoCRs) cout<<"average CR2 sf for DYtautau:  "<<CRsftot[CR2]<<" +/- "<<CRsftoterr[CR2]<<endl;
    if(scalebkgtoCRs) cout<<"composite sf for DYtautau (CR1):  "<<(CRsftot[CR1]/CRsftot[CR1v])*CRsftot[CR2v]<<" +/- "<<(CRsftot[CR1]/CRsftot[CR1v])*CRsftot[CR2v]*sqrt(pow(CRsftoterr[CR1]/CRsftot[CR1],2)+pow(CRsftoterr[CR1v]/CRsftot[CR1v],2)+pow(CRsftoterr[CR2v]/CRsftot[CR2v],2))<<endl;
    if(scalebkgtoCRs) cout<<"composite sf for DYtautau (CR3):  "<<(CRsftot[CR3]/CRsftot[CR3v])*CRsftot[CR2v]<<" +/- "<<(CRsftot[CR3]/CRsftot[CR3v])*CRsftot[CR2v]*sqrt(pow(CRsftoterr[CR3]/CRsftot[CR3],2)+pow(CRsftoterr[CR3v]/CRsftot[CR3v],2)+pow(CRsftoterr[CR2v]/CRsftot[CR2v],2))<<endl;
    if(scalebkgtoCRs) cout<<"CR4v sf for DYtautau (exc. b modelling SF):  "<<CRsftot[CR4v]<<" +/- "<<CRsftoterr[CR4v]<<endl;
    //if(scalebkgtoCRs) cout<<"CR40 sf for DYtautau:  "<<CRsftot[CR40]<<" +/- "<<CRsftoterr[CR40]<<endl;


    TGraphErrors *graphDYeemm = new TGraphErrors();
    TGraphErrors *graphDYtt = new TGraphErrors();
    TGraphErrors *graphfakes = new TGraphErrors();


    if(scalebkgtoCRs) {

        cout<<"fakes:"<<endl;

        for (int leptype = 0; leptype < nCh; ++leptype)
        {
            cout<<CRsf[CR5][leptype]<<"+/-"<<CRsferr[CR5][leptype]<<endl;
            cout<<CRsf[CR5v][leptype]<<"+/-"<<CRsferr[CR5v][leptype]<<endl;
            cout<<CRsf[CR6][leptype]<<"+/-"<<CRsferr[CR6][leptype]<<endl;
        }

        for (int leptype = 0; leptype < nCh; ++leptype)
        {
            graphfakes->SetPoint(leptype*3+0, leptype*3+0+0.5, CRsf[CR5][leptype]); graphfakes->SetPointError( leptype*3+0, 0, CRsferr[CR5][leptype] );
            graphfakes->SetPoint(leptype*3+1, leptype*3+1+0.5, CRsf[CR5v][leptype]); graphfakes->SetPointError( leptype*3+1, 0, CRsferr[CR5v][leptype] );
            graphfakes->SetPoint(leptype*3+2, leptype*3+2+0.5, CRsf[CR6][leptype]); graphfakes->SetPointError( leptype*3+2, 0, CRsferr[CR6][leptype] );
        }






        cout<<"DY:"<<endl;

        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR1][leptype]<<"+/-"<<CRsferr[CR1][leptype]<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR0][leptype]*CRsf[CR1][leptype]/CRsf[CR1v][leptype]<<" +/- "<<(CRsf[CR0][leptype]*CRsf[CR1][leptype]/CRsf[CR1v][leptype])*sqrt(pow(CRsferr[CR0][leptype]/CRsf[CR0][leptype],2)+pow(CRsferr[CR1][leptype]/CRsf[CR1][leptype],2)+pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2))<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR0][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype]<<" +/- "<<(CRsf[CR0][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype])*sqrt(pow(CRsferr[CR0][leptype]/CRsf[CR0][leptype],2)+pow(CRsferr[CR2][leptype]/CRsf[CR2][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2))<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR0][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype]<<" +/- "<<(CRsf[CR0][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype])*sqrt(pow(CRsferr[CR0][leptype]/CRsf[CR0][leptype],2)+pow(CRsferr[CR3][leptype]/CRsf[CR3][leptype],2)+pow(CRsferr[CR3v][leptype]/CRsf[CR3v][leptype],2))<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR1v][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype]<<" +/- "<<(CRsf[CR1v][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype])*sqrt(pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2)+pow(CRsferr[CR2][leptype]/CRsf[CR2][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2))<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR1v][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype]<<" +/- "<<(CRsf[CR1v][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype])*sqrt(pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2)+pow(CRsferr[CR3][leptype]/CRsf[CR3][leptype],2)+pow(CRsferr[CR3v][leptype]/CRsf[CR3v][leptype],2))<<endl;


        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR2][leptype]<<"+/-"<<CRsferr[CR2][leptype]<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<(CRsf[CR1][leptype]/CRsf[CR1v][leptype])*CRsf[CR2v][leptype]<<" +/- "<<(CRsf[CR1][leptype]/CRsf[CR1v][leptype])*CRsf[CR2v][leptype]*sqrt(pow(CRsferr[CR1][leptype]/CRsf[CR1][leptype],2)+pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2))<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<(CRsf[CR3][leptype]/CRsf[CR3v][leptype])*CRsf[CR2v][leptype]<<" +/- "<<(CRsf[CR3][leptype]/CRsf[CR3v][leptype])*CRsf[CR2v][leptype]*sqrt(pow(CRsferr[CR3][leptype]/CRsf[CR3][leptype],2)+pow(CRsferr[CR3v][leptype]/CRsf[CR3v][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2))<<endl;
        for (int leptype = 0; leptype < 2; ++leptype) cout<<CRsf[CR4v][2]*CRsf[CR2][leptype]/CRsf[CR2v][leptype]<<" +/- "<<(CRsf[CR4v][2]*CRsf[CR2][leptype]/CRsf[CR2v][leptype])*sqrt(pow(CRsferr[CR4v][2]/CRsf[CR4v][2],2)+pow(CRsferr[CR2][leptype]/CRsf[CR2][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2))<<endl;


        for (int leptype = 0; leptype < 2; ++leptype) { graphDYeemm->SetPoint( 1+leptype-1, 1+leptype-1+0.5, CRsf[CR1][leptype] );    graphDYeemm->SetPointError( 1+leptype-1, 0, CRsferr[CR1][leptype] ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYeemm->SetPoint( 3+leptype-1, 3+leptype-1+0.5, CRsf[CR0][leptype]*CRsf[CR1][leptype]/CRsf[CR1v][leptype] );    graphDYeemm->SetPointError( 3+leptype-1, 0, (CRsf[CR0][leptype]*CRsf[CR1][leptype]/CRsf[CR1v][leptype])*sqrt(pow(CRsferr[CR0][leptype]/CRsf[CR0][leptype],2)+pow(CRsferr[CR1][leptype]/CRsf[CR1][leptype],2)+pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2)) ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYeemm->SetPoint( 5+leptype-1, 5+leptype-1+0.5, CRsf[CR0][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype] );    graphDYeemm->SetPointError( 5+leptype-1, 0, (CRsf[CR0][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype])*sqrt(pow(CRsferr[CR0][leptype]/CRsf[CR0][leptype],2)+pow(CRsferr[CR2][leptype]/CRsf[CR2][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2)) ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYeemm->SetPoint( 7+leptype-1, 7+leptype-1+0.5, CRsf[CR0][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype] );    graphDYeemm->SetPointError( 7+leptype-1, 0, (CRsf[CR0][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype])*sqrt(pow(CRsferr[CR0][leptype]/CRsf[CR0][leptype],2)+pow(CRsferr[CR3][leptype]/CRsf[CR3][leptype],2)+pow(CRsferr[CR3v][leptype]/CRsf[CR3v][leptype],2)) ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYeemm->SetPoint( 9+leptype-1, 9+leptype-1+0.5, CRsf[CR1v][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype] );    graphDYeemm->SetPointError( 9+leptype-1, 0, (CRsf[CR1v][leptype]*CRsf[CR2][leptype]/CRsf[CR2v][leptype])*sqrt(pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2)+pow(CRsferr[CR2][leptype]/CRsf[CR2][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2)) ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYeemm->SetPoint( 11+leptype-1, 11+leptype-1+0.5, CRsf[CR1v][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype] );    graphDYeemm->SetPointError( 11+leptype-1, 0, (CRsf[CR1v][leptype]*CRsf[CR3][leptype]/CRsf[CR3v][leptype])*sqrt(pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2)+pow(CRsferr[CR3][leptype]/CRsf[CR3][leptype],2)+pow(CRsferr[CR3v][leptype]/CRsf[CR3v][leptype],2)) ); }


        for (int leptype = 0; leptype < 2; ++leptype) { graphDYtt->SetPoint( 1+leptype-1, 1+leptype-1+0.5, CRsf[CR2][leptype] );    graphDYtt->SetPointError( 1+leptype-1, 0, CRsferr[CR2][leptype] ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYtt->SetPoint( 3+leptype-1, 3+leptype-1+0.5, (CRsf[CR1][leptype]/CRsf[CR1v][leptype])*CRsf[CR2v][leptype] );    graphDYtt->SetPointError( 3+leptype-1, 0, (CRsf[CR1][leptype]/CRsf[CR1v][leptype])*CRsf[CR2v][leptype]*sqrt(pow(CRsferr[CR1][leptype]/CRsf[CR1][leptype],2)+pow(CRsferr[CR1v][leptype]/CRsf[CR1v][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2)) ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYtt->SetPoint( 5+leptype-1, 5+leptype-1+0.5, (CRsf[CR3][leptype]/CRsf[CR3v][leptype])*CRsf[CR2v][leptype] );    graphDYtt->SetPointError( 5+leptype-1, 0, (CRsf[CR3][leptype]/CRsf[CR3v][leptype])*CRsf[CR2v][leptype]*sqrt(pow(CRsferr[CR3][leptype]/CRsf[CR3][leptype],2)+pow(CRsferr[CR3v][leptype]/CRsf[CR3v][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2)) ); }
        for (int leptype = 0; leptype < 2; ++leptype) { graphDYtt->SetPoint( 7+leptype-1, 7+leptype-1+0.5, CRsf[CR4v][2]*CRsf[CR2][leptype]/CRsf[CR2v][leptype] );    graphDYtt->SetPointError( 7+leptype-1, 0, (CRsf[CR4v][2]*CRsf[CR2][leptype]/CRsf[CR2v][leptype])*sqrt(pow(CRsferr[CR4v][2]/CRsf[CR4v][2],2)+pow(CRsferr[CR2][leptype]/CRsf[CR2][leptype],2)+pow(CRsferr[CR2v][leptype]/CRsf[CR2v][leptype],2)) ); }



        TH2F* hdummyDYeemm = new TH2F("hdummyDYeemm","",12,0,12,100,0.9,1.8);
        
        hdummyDYeemm->GetXaxis()->SetBinLabel(1,"CR1 ee");
        hdummyDYeemm->GetXaxis()->SetBinLabel(1+1,"CR1 #mu#mu");
        hdummyDYeemm->GetXaxis()->SetBinLabel(3,"CR0#times1 ee");
        hdummyDYeemm->GetXaxis()->SetBinLabel(3+1,"CR0#times1 #mu#mu");
        hdummyDYeemm->GetXaxis()->SetBinLabel(5,"CR0#times2 ee");
        hdummyDYeemm->GetXaxis()->SetBinLabel(5+1,"CR0#times2 #mu#mu");
        hdummyDYeemm->GetXaxis()->SetBinLabel(7,"CR0#times3 ee");
        hdummyDYeemm->GetXaxis()->SetBinLabel(7+1,"CR0#times3 #mu#mu");
        hdummyDYeemm->GetXaxis()->SetBinLabel(9,"CR1v#times2 ee");
        hdummyDYeemm->GetXaxis()->SetBinLabel(9+1,"CR1v#times2 #mu#mu");
        hdummyDYeemm->GetXaxis()->SetBinLabel(11,"CR1v#times3 ee");
        hdummyDYeemm->GetXaxis()->SetBinLabel(11+1,"CR1v#times3 #mu#mu");
        
        hdummyDYeemm->GetXaxis()->SetTitle("Control Region");
        hdummyDYeemm->GetXaxis()->SetTitleOffset(1.05);
        hdummyDYeemm->GetYaxis()->SetTitle("DY #rightarrow ee/#mu#mu SF");
        
        
        TH1F* hupDYeemm = new TH1F("hupDYeemm","",12,0,12);
        TH1F* hcnDYeemm = new TH1F("hcnDYeemm","",12,0,12);
        TH1F* hdnDYeemm = new TH1F("hdnDYeemm","",12,0,12);
        for (int ibin = 0; ibin < 12; ++ibin)
        {
            hupDYeemm->SetBinContent(ibin+1,CRsftot[CR1]+0.2);
            hcnDYeemm->SetBinContent(ibin+1,CRsftot[CR1]);
            hdnDYeemm->SetBinContent(ibin+1,CRsftot[CR1]-0.2);
        }
        
        hupDYeemm->SetLineWidth(3);
        hupDYeemm->SetLineColor(kRed);
        hcnDYeemm->SetLineWidth(3);
        hcnDYeemm->SetLineColor(kBlue);
        hdnDYeemm->SetLineWidth(3);
        hdnDYeemm->SetLineColor(kRed);
    



        TH2F* hdummyDYtt = new TH2F("hdummyDYtt","",8,0,8,100,1.0,1.35);
        
        hdummyDYtt->GetXaxis()->SetBinLabel(1,"CR2 ee");
        hdummyDYtt->GetXaxis()->SetBinLabel(1+1,"CR2 #mu#mu");
        hdummyDYtt->GetXaxis()->SetBinLabel(3,"CR1/v ee");
        hdummyDYtt->GetXaxis()->SetBinLabel(3+1,"CR1/v #mu#mu");
        hdummyDYtt->GetXaxis()->SetBinLabel(5,"CR3/v ee");
        hdummyDYtt->GetXaxis()->SetBinLabel(5+1,"CR3/v #mu#mu");
        hdummyDYtt->GetXaxis()->SetBinLabel(7,"CR4v#times2/v ee");
        hdummyDYtt->GetXaxis()->SetBinLabel(7+1,"CR4v#times2/v #mu#mu");

        
        hdummyDYtt->GetXaxis()->SetTitle("Control Region");
        hdummyDYtt->GetXaxis()->SetTitleOffset(1.05);
        hdummyDYtt->GetYaxis()->SetTitle("DY #rightarrow #tau#tau SF");
        
        
        TH1F* hupDYtt = new TH1F("hupDYtt","",8,0,8);
        TH1F* hcnDYtt = new TH1F("hcnDYtt","",8,0,8);
        TH1F* hdnDYtt = new TH1F("hdnDYtt","",8,0,8);
        for (int ibin = 0; ibin < 8; ++ibin)
        {
            hupDYtt->SetBinContent(ibin+1,CRsftot[CR2]+0.1);
            hcnDYtt->SetBinContent(ibin+1,CRsftot[CR2]);
            hdnDYtt->SetBinContent(ibin+1,CRsftot[CR2]-0.1);
        }
        
        hupDYtt->SetLineWidth(3);
        hupDYtt->SetLineColor(kRed);
        hcnDYtt->SetLineWidth(3);
        hcnDYtt->SetLineColor(kBlue);
        hdnDYtt->SetLineWidth(3);
        hdnDYtt->SetLineColor(kRed);
    



        TH2F* hdummyfakes = new TH2F("hdummyfakes","",9,0,9,100,0.5,4.0);
        
        hdummyfakes->GetXaxis()->SetBinLabel(1,"CR5 ee");
        hdummyfakes->GetXaxis()->SetBinLabel(2,"CR5v ee");
        hdummyfakes->GetXaxis()->SetBinLabel(3,"CR6 ee");
        hdummyfakes->GetXaxis()->SetBinLabel(1+3,"CR5 #mu#mu");
        hdummyfakes->GetXaxis()->SetBinLabel(2+3,"CR5v #mu#mu");
        hdummyfakes->GetXaxis()->SetBinLabel(3+3,"CR6 #mu#mu");
        hdummyfakes->GetXaxis()->SetBinLabel(1+6,"CR5 e#mu");
        hdummyfakes->GetXaxis()->SetBinLabel(2+6,"CR5v e#mu");
        hdummyfakes->GetXaxis()->SetBinLabel(3+6,"CR6 e#mu");

        
        hdummyfakes->GetXaxis()->SetTitle("Control Region");
        hdummyfakes->GetXaxis()->SetTitleOffset(1.05);
        hdummyfakes->GetYaxis()->SetTitle("Fakes SF");
        
        
        TH1F* hupfakes = new TH1F("hupfakes","",9,0,9);
        TH1F* hcnfakes = new TH1F("hcnfakes","",9,0,9);
        TH1F* hdnfakes = new TH1F("hdnfakes","",9,0,9);
        for (int ibin = 0; ibin < 9; ++ibin)
        {
            hupfakes->SetBinContent(ibin+1,CRsftot[CR5]+1.);
            hcnfakes->SetBinContent(ibin+1,CRsftot[CR5]);
            hdnfakes->SetBinContent(ibin+1,CRsftot[CR5]-1.);
        }
        
        hupfakes->SetLineWidth(3);
        hupfakes->SetLineColor(kRed);
        hcnfakes->SetLineWidth(3);
        hcnfakes->SetLineColor(kBlue);
        hdnfakes->SetLineWidth(3);
        hdnfakes->SetLineColor(kRed);
    











        TCanvas *canSFs = new TCanvas("canSFs","canSFs",1800, 600);
        canSFs->Divide(3,1);
        canSFs->cd(1);
        //gPad->SetGridx();
        //gPad->SetGridy();
        hdummyDYeemm->Draw();
        hupDYeemm->Draw("samehist");
        hcnDYeemm->Draw("samehist");
        hdnDYeemm->Draw("samehist");
        graphDYeemm->Draw("sameP");
        canSFs->cd(2);
        hdummyDYtt->Draw();
        hupDYtt->Draw("samehist");
        hcnDYtt->Draw("samehist");
        hdnDYtt->Draw("samehist");
        graphDYtt->Draw("sameP"); 
        canSFs->cd(3);
        hdummyfakes->Draw();
        hupfakes->Draw("samehist");
        hcnfakes->Draw("samehist");
        hdnfakes->Draw("samehist");
        graphfakes->Draw("sameP"); 
        canSFs->Print("SFs.pdf");


    }

    //initialize

    TFile *dt_dl[nCh];
    dt_dl[DIEL] = TFile::Open(Form("%soutput/data_diel_histos.root",region));
    dt_dl[DIMU] = TFile::Open(Form("%soutput/data_dimu_histos.root",region));
    dt_dl[MUEG] = TFile::Open(Form("%soutput/data_mueg_histos.root",region));

    TFile *mc_dl[MCID];

    for (int j = 0; j < MCID; ++j)
    {
        if (j < 2)
            mc_dl[j] = TFile::Open(Form("%soutput/%s_%s_histos.root", region,
                                        mcsample[j], ( j != TTDL && ttbar_string.Contains("FullLept") ? "mcatnlo" : ttbar_tag )));
        else
            mc_dl[j] = TFile::Open(Form("%soutput/%s_histos.root", region,
                                        mcsample[j]));
    }




    //for (int isr = 0; isr < NSAMPLE; ++isr)
    //{

    const int N1DHISTS = 60;
    TH1F *h_dt1d_comb[N1DHISTS];
    TH1F *h_mc1d_comb[N1DHISTS][MCID];
    vector<TH1F *> sorted_mc1d_comb[N1DHISTS];
    vector<TH1F *> unsorted_mc1d_comb[N1DHISTS];
    TH1F *h_mc1d_tot_comb[N1DHISTS];
    bool hasplot[nCh][N1DHISTS] = {{0}};
    bool hasplotall[N1DHISTS] = {0};

    // Output file names - don't change order of the first 14
    const char *file1dname[N1DHISTS] =
    {
        "lep_azimuthal_asymmetry2",
        "lep_azimuthal_asymmetry",
        "lep_charge_asymmetry",
        "top_rapiditydiff_Marco",
        "lepPlus_costheta_cms",
        "lepMinus_costheta_cms",
        "top_spin_correlation",
        "lep_cos_opening_angle",
        "m_top",
        "tt_mass",
        "tt_pT",
        "ttRapidity2",

        "channel",
        "channel_withttbarsol",
        "n_jets",
        "n_bjets",
        "met",
        "metphi",
        "met_smeared",
        "lep1b_mindR",
        "lep1b_mindPhi",
        "lep2b_mindR",
        "lep2b_mindPhi",
        "lep_dR",
        "lep_dEta",
        "Mll",
        "Etall",
        "Phill",
        "Ptll",
        "lepPlus_Pt",
        "lepMinus_Pt",
        "lepPt",
        //"lepPt_ele",
        //"lepPt_muo",
        "jet_Pt",
        "jet0_Pt",
        "jet1_Pt",
        "b_Pt",
        "b0_Pt",
        "b1_Pt",
        "nonb_Pt",
        "lepPlus_Eta",
        "lepMinus_Eta",
        "lepEta",
        //"lepEta_ele",
        //"lepEta_muo",
        "jet_Eta",
        "jet0_Eta",
        "jet1_Eta",
        "b_Eta",
        "b0_Eta",
        "b1_Eta",
        "nonb_Eta",
        //"pfcalo_metratio",
        //"pfcalo_metratio2",
        //"pfcalodPhi",
        //"pfcalo_deltamet",
        //"pfcalo_deltametx",
        //"pfcalo_deltamety",
        //"calomet",
        //"calometphi",
        "mlb",
        "mlb_min",
        "maxAMWTweight",
        "maxAMWTweight_closestApproach",
        "otherAMWTweights",
        "closestDeltaMET_bestcombo",
        "nvtx",
        "top1_pt",
        "top2_pt",
        "top1_p_CM",
        "lep_costheta_cms"
    };


    // List of Log scale plots:
    vector<int> logScale;
    //logScale.push_back(0);
    //logScale.push_back(1);
    //logScale.push_back(2);
    //logScale.push_back(4);
    //logScale.push_back(5);
    //logScale.push_back(6);
    //logScale.push_back(7);
    logScale.push_back(8);
    //logScale.push_back(16);
    //logScale.push_back(56);
    //logScale.push_back(33+11);
    //logScale.push_back(34+11);
    //logScale.push_back(35+11);
    //logScale.push_back(36+11);
    logScale.push_back(41+8);
    logScale.push_back(42+8);

    // List of rebin factors:
    int rebinFactor[N1DHISTS] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    const char *xtitle1d[N1DHISTS] =
    {
        "#Delta#phi_{l#lower[-0.4]{+}l#lower[-0.48]{-}}",
        "#Delta#phi_{l#lower[-0.4]{+}l#lower[-0.48]{-}}",
        "#Delta|#eta_{l}|",
        "#Delta|y_{t}|",
        "cos(^{}#theta_{l#lower[-0.4]{+}}#kern[-1.38]{*})",
        "cos(^{}#theta_{l#lower[-0.48]{-}}#kern[-1.0]{*})",
        "cos(^{}#theta_{l#lower[-0.4]{+}}#kern[-1.38]{*}) cos(^{}#theta_{l#lower[-0.48]{-}}#kern[-1.0]{*})",
        "cos(#phi)",
        "M_{t}",
        "M_{t#bar{t}}",
        "p_{T}^{t#bar{t}}",
        "y_{t#bar{t}}",
        
        "channel",
        "channel (after solution)",
        "n_jets",
        "n_bjets",
        "ME_{T} [GeV]",
        "metphi",
        "ME_{T}^{smeared} [GeV]",
        "lep1b_mindR",
        "lep1b_mindPhi",
        "lep2b_mindR",
        "lep2b_mindPhi",
        "lep_dR",
        "lep_dEta",
        "Mll",
        "Etall",
        "Phill",
        "Ptll",
        "lepPlus_Pt",
        "lepMinus_Pt",
        "lepPt",
        //"lepPt_ele",
        //"lepPt_muo",
        "jet_Pt",
        "jet0_Pt",
        "jet1_Pt",
        "b_Pt",
        "b0_Pt",
        "b1_Pt",
        "nonb_Pt",
        "lepPlus_Eta",
        "lepMinus_Eta",
        "lepEta",
        //"lepEta_ele",
        //"lepEta_muo",
        "jet_Eta",
        "jet0_Eta",
        "jet1_Eta",
        "b_Eta",
        "b0_Eta",
        "b1_Eta",
        "nonb_Eta",
        //"pfcalo_metratio",
        //"pfcalo_metratio2",
        //"pfcalodPhi",
        //"pfcalo_deltamet",
        //"pfcalo_deltametx",
        //"pfcalo_deltamety",
        //"calomet",
        //"calometphi",
        "mlb",
        "mlb_min",
        "maxAMWTweight",
        "maxAMWTweight_closestApproach",
        "otherAMWTweights",
        "closestDeltaMET_bestcombo",
        "nvtx",
        "top1_pt",
        "top2_pt",
        "top1_p_CM",
        "lep_costheta_cms"
    };

    const char *ytitle1d[N1DHISTS] =
    {
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "GeV",
        "GeV",
        "GeV",
        "",

        //"",
        //"",
        //"",
        //"",
        //"",
        //"",
        "",
        "",
        "",
        "",
        "GeV",
        "",
        "GeV",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        "",
        ""
    };


    //not used currently
    /*
    for (int leptype = 0; leptype < nCh; ++leptype)
    {
        TH1F *h_pv_mc_mt[MCID];
        TH1F *h_pv_mc_mt_sf;//single lepton sample
        TH1F *h_pv_mc_mt_nosf;//non single lepton nor dilepton
        TH1F *h_pv_mc_mt_tot;
        TH1F *h_pv_tmp;

        histoname = "h_sig_m_top";

        //data
        TH1F *h_pv_dt_mt = (TH1F *)dt_dl[leptype]->Get(Form("%s%s%s",
                           histoname, metcut[isr],
                           leptag[leptype]));
        h_pv_dt_mt->SetName("h_pv_dt_mt");

        //mc
        for (int j = 0; j < MCID; ++j)
        {

            h_pv_mc_mt[j] = (TH1F *)mc_dl[j]->Get(Form("%s%s%s", histoname, metcut[isr],
                                                  leptag[leptype]));
            if (h_pv_mc_mt[j] == 0)
            {
                h_pv_mc_mt[j] =
                    (TH1F *)mc_dl[0]->Get(Form("%s%s%s", histoname, metcut[isr],
                                               leptag[leptype]));
                h_pv_mc_mt[j]->SetName(Form("h_pv_mc_mt_%s", mcsample[j]));
                zeroHist(h_pv_mc_mt[j]);
            }
            h_pv_mc_mt[j]->SetName(Form("h_pv_mc_mt_%s", mcsample[j]));

            if (j == TTDL)
            {
                //default combined dilepton
                h_pv_mc_mt_tot = (TH1F *)h_pv_mc_mt[j]->Clone("h_pv_mc_mt_tot");
            }
            else
            {
                h_pv_mc_mt_tot->Add(h_pv_mc_mt[j]);

                //single lepton samples
                if (j == TTSL)
                    h_pv_mc_mt_sf = (TH1F *)h_pv_mc_mt[j]->Clone("h_pv_mc_mt_sf");
                else if (j == WJETS || j == TWSL)
                    h_pv_mc_mt_sf->Add(h_pv_mc_mt[j]);

                if (j == TWDL)
                    h_pv_mc_mt_nosf = (TH1F *)h_pv_mc_mt[j]->Clone("h_pv_mc_mt_nosf");
                else if (j == DIBO || j == DY || j == TTV || j == VVV)
                    h_pv_mc_mt_nosf->Add(h_pv_mc_mt[j]);
            }
            if (doverbose)
            {
                cout << "-------------------------------------------------" << endl;
                printf("%s MC Yields per region \n", mcsample[j]);
                printf("Control: %.2f pm %.2f (%.1f %%) \n",
                       h_pv_mc_mt[j]->GetBinContent(i_ctr), h_pv_mc_mt[j]->GetBinError(1),
                       h_pv_mc_mt[j]->GetBinError(i_ctr) / h_pv_mc_mt[j]->GetBinContent(i_ctr) * 100.);
                printf("Signal: %.2f pm %.2f (%.1f %%) \n",
                       h_pv_mc_mt[j]->GetBinContent(i_sig), h_pv_mc_mt[j]->GetBinError(i_sig),
                       h_pv_mc_mt[j]->GetBinError(i_sig) / h_pv_mc_mt[j]->GetBinContent(i_sig) * 100.);
            }
        }

    }//end loop over lepton types
    */




    for (int leptype = 0; leptype < nCh; ++leptype)
    {

        // End MC Samples definition and open files ///////////////////////////////
        // Open and book histograms ///////////////////////////////////////////////

        TH1F *h_dt1d[N1DHISTS];
        TH1F *h_mc1d[N1DHISTS][MCID];
        vector<TH1F *> sorted_mc1d[N1DHISTS];

        // Book histograms
        TH1F *h_mc1d_tot[N1DHISTS];

        cout << "#### LEPTON TAG " << leptag[leptype] << endl;

        for (int i = 0; i < N1DHISTS; ++i)
        {
            h_dt1d[i] = (TH1F *)dt_dl[leptype]->Get(Form("%s%s%s%s", histtag[histtag_number], file1dname[i], metcut[isr], leptag[leptype]));

            //cout << "DT: " << h_dt1d[i] << " " << Form("%s%s%s", histtag[histtag_number], file1dname[i], leptag[leptype]) << endl;

            if (h_dt1d[i] != 0)
            {
                hasplot[leptype][i] = true;
                //hasplotall[i] = true;
            }
            else
            {
                cout << Form("%s%s%s", histtag[histtag_number], file1dname[i], leptag[leptype]) << " plot not found, skipping." << endl;
                continue;
            }

            //cout << "gtl: " << __LINE__ << endl;

            h_dt1d[i]->SetName(Form("%s%s", histtag[histtag_number], file1dname[i]));

            h_dt1d[i]->Rebin(rebinFactor[i]);


            bool doinit = false;
            for (int j = 0; j < MCID; ++j)
            {

                h_mc1d[i][j] =
                    (TH1F *)mc_dl[j]->Get(Form("%s%s%s%s", histtag[histtag_number], file1dname[i], metcut[isr], leptag[leptype]));

                //cout << "MC" << j << " " << mcsample[j] << ": " << h_mc1d[i][j] << " "
                //     << Form("%s%s", histtag[histtag_number], file1dname[i]) << endl;

                if (h_mc1d[i][j] == 0)
                {
                    h_mc1d[i][j] =
                        //(TH1F *)mc_dl[0]->Get(Form("%s%s%s%s", histtag[histtag_number], file1dname[i], metcut[isr], leptag[leptype]));
                        (TH1F *)dt_dl[leptype]->Get(Form("%s%s%s%s", histtag[histtag_number], file1dname[i], metcut[isr], leptag[leptype]));
                    h_mc1d[i][j]->SetName(Form("%s_%s%s", mcsample[j], histtag[histtag_number], file1dname[i]));
                    zeroHist(h_mc1d[i][j]);
                }
                h_mc1d[i][j]->SetName(Form("%s_%s%s", mcsample[j], histtag[histtag_number], file1dname[i]));

                //cout << "MC" << j << " " << mcsample[j] << ": " << h_mc1d[i][j] << " "
                //     << Form("%s%s", histtag[histtag_number], file1dname[i]) << endl;

                //h_mc1d[i][j]->Sumw2();
                //if(j==TTDL||j==TTSL||j==TTFA) h_mc1d[i][j]->Scale(40033./4884387.); //scale to powheg
                //if(j==TTDL||j==TTSL||j==TTFA) h_mc1d[i][j]->Scale(303732. / 32852589.); //scale to xsec


                if(scalebkgtoCRs && j!=TTDL) h_mc1d[i][j]->Scale(bkgsf[SFsApply[histtag_number]][j][leptype]);

                // Rebin, set limits
                h_mc1d[i][j]->Rebin(rebinFactor[i]);

                if (!doinit)
                {
                    h_mc1d_tot[i] = (TH1F *)h_mc1d[i][j]->Clone(Form("mctot_%s%s", histtag[histtag_number], file1dname[i]));
                    doinit = true;
                }
                else h_mc1d_tot[i]->Add(h_mc1d[i][j]);

            }
        }


        //applying ttdil scaling (for SIG do exact scaling, for CRs do scaling based on inclusive SIG region)
        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (!hasplot[leptype][i]) continue;

            //cout << "-------------------------------------------------" << endl;
            float mcall = h_mc1d_tot[i]->Integral();
            float dtall = h_dt1d[i]->Integral();
            //cout << "MC ALL " << mcall << endl;
            //cout << "Data ALL " << dtall << endl;
            float mcsf_all = dtall / mcall;
            //cout << "RATIO ALL " << mcsf_all << endl;
            //cout << "-------------------------------------------------" << endl;
            mcsf_all = 1; //don't scale all MC to data, instead just the ttdil component is scaled below
            h_mc1d_tot[i]->Scale(mcsf_all);

            float mcttdil_all = h_mc1d[i][TTDL]->Integral();
            float ttdilsf_all = 1. + (dtall-mcall)/mcttdil_all;
            if(scalettdiltodata && histtag_number==0) h_mc1d[i][TTDL]->Scale(ttdilsf_all);
            else if(scalettdiltodata) h_mc1d[i][TTDL]->Scale(bkgsf[1][TTDL][leptype]); 

        }

        //(re)calculate _tot (sum of MC) histograms
        for (int i = 0; i < N1DHISTS; ++i)
        {

            bool combplotcreated = hasplotall[i];

            if ( !hasplot[leptype][i] ) continue;
            hasplotall[i] = true;

            //e+mu combined data
            if ( combplotcreated == false )
            {
                h_dt1d_comb[i] = (TH1F *)h_dt1d[i]->Clone();
                h_dt1d_comb[i]->SetName(Form("%s%s_comb", histtag[histtag_number], file1dname[i]));
            }
            else
                h_dt1d_comb[i]->Add(h_dt1d[i]);

            bool doinit = false;
            for (int j = 0; j < MCID; ++j)
            {
                //combined MC histograms
                if (combplotcreated == false)
                {
                    h_mc1d_comb[i][j] = (TH1F *)h_mc1d[i][j]->Clone();
                    h_mc1d_comb[i][j]->SetName(Form("%s_%s%s_comb", mcsample[j], histtag[histtag_number], file1dname[i]));
                }
                else
                    h_mc1d_comb[i][j]->Add(h_mc1d[i][j]);

                //get total
                if (!doinit && combplotcreated == false)
                {
                    h_mc1d_tot_comb[i] = (TH1F *)h_mc1d[i][j]->Clone(Form("mctot_comb_%s%s", histtag[histtag_number], file1dname[i]));
                    h_mc1d_tot_comb[i]->SetName(Form("mctot_comb_%s%s", histtag[histtag_number], file1dname[i]));
                }
                else
                    h_mc1d_tot_comb[i]->Add(h_mc1d[i][j]);

                if (!doinit)
                {
                    h_mc1d_tot[i] = (TH1F *)h_mc1d[i][j]->Clone(Form("mctot_%s%s", histtag[histtag_number], file1dname[i]));
                    doinit = true;
                }
                else h_mc1d_tot[i]->Add(h_mc1d[i][j]);

            }
        }

        //after applying all the scalings, fill vector for sorting
        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (!hasplot[leptype][i]) continue;

            /*
                            if (i == 1) //mt
                            {
                                int minbin = h_mc1d_tot[i]->FindBin(160);
                                int maxbin = h_mc1d_tot[i]->GetNbinsX();

                                mc_mttail  = h_mc1d_tot[i]->Integral(minbin, maxbin + 1);
                                emc_mttail = histError( h_mc1d_tot[i], minbin, maxbin + 1);
                                dt_mttail  = h_dt1d[i]->Integral(minbin, maxbin + 1);
                                edt_mttail = histError( h_dt1d[i], minbin, maxbin + 1);
                                cout << "MT: Going from " << minbin << " to " << maxbin
                                     << " nbin range " << h_mc1d_tot[i]->GetXaxis()->GetBinLowEdge(minbin)
                                     << " - " << h_mc1d_tot[i]->GetXaxis()->GetBinUpEdge(maxbin) << endl;
                                cout << "MC: " << mc_mttail << " pm " << emc_mttail << endl;
                                cout << "Data: " << dt_mttail << " pm " << edt_mttail << endl;
                                sf_dtmc_mttail  = dt_mttail / mc_mttail;
                                esf_dtmc_mttail = sqrt(pow( (edt_mttail / dt_mttail), 2) +
                                                       pow( (emc_mttail / mc_mttail), 2 )) * sf_dtmc_mttail;
                                cout << "SF: " << sf_dtmc_mttail << " pm " << esf_dtmc_mttail << endl;
                            }
                            else if (i == 2)  //lepton pt
                            {
                                int minbin = h_mc1d_tot[i]->FindBin(150);
                                int maxbin = h_mc1d_tot[i]->GetNbinsX();

                                mc_lpttail  = h_mc1d_tot[i]->Integral(minbin, maxbin + 1);
                                emc_lpttail = histError( h_mc1d_tot[i], minbin, maxbin + 1);
                                dt_lpttail  = h_dt1d[i]->Integral(minbin, maxbin + 1);
                                edt_lpttail = histError( h_dt1d[i], minbin, maxbin + 1);
                                cout << "LEPTON PT: Going from " << minbin << " to " << maxbin
                                     << " nbin range " << h_mc1d_tot[i]->GetXaxis()->GetBinLowEdge(minbin)
                                     << " - " << h_mc1d_tot[i]->GetXaxis()->GetBinUpEdge(maxbin) << endl;
                                cout << "MC: " << mc_lpttail << " pm " << emc_lpttail << endl;
                                cout << "Data: " << dt_lpttail << " pm " << edt_lpttail << endl;
                                sf_dtmc_lpttail  = dt_lpttail / mc_lpttail;
                                esf_dtmc_lpttail = sqrt(pow( (edt_lpttail / dt_lpttail), 2) +
                                                        pow( (emc_lpttail / mc_lpttail), 2 )) * sf_dtmc_lpttail;
                                cout << "SF: " << sf_dtmc_lpttail << " pm " << esf_dtmc_lpttail << endl;
                            }
                            else if (i == 3)  //mt count
                            {
                                mc_mtctail[isr][leptype]  = h_mc1d_tot[i]->GetBinContent(2);
                                emc_mtctail[isr][leptype] = h_mc1d_tot[i]->GetBinError(2);
                                dt_mtctail[isr][leptype]  = h_dt1d[i]->GetBinContent(2);
                                edt_mtctail[isr][leptype] = h_dt1d[i]->GetBinError(2);
                                cout << "MT TAIL COUNT: " << endl;
                                cout << "MC: " << mc_mtctail[isr][leptype] << " pm " << emc_mtctail[isr][leptype] << endl;
                                cout << "Data: " << dt_mtctail[isr][leptype] << " pm " << edt_mtctail[isr][leptype] << endl;
                                sf_dtmc_mtctail[isr][leptype]  = dt_mtctail[isr][leptype] / mc_mtctail[isr][leptype];
                                esf_dtmc_mtctail[isr][leptype] = sqrt(pow( (edt_mtctail[isr][leptype] / dt_mtctail[isr][leptype]), 2) +
                                                                      pow( (emc_mtctail[isr][leptype] / mc_mtctail[isr][leptype]), 2 )) * sf_dtmc_mtctail[isr][leptype];
                                cout << "SF: " << sf_dtmc_mtctail[isr][leptype] << " pm " << esf_dtmc_mtctail[isr][leptype] << endl;
                            }
            */
            for (int j = 0; j < MCID; ++j)
            {

                //version for sorting
                //h_mc1d[i][j]->Scale(mcsf_all);
                /*
                                    if (j >= TWDL && j != VVV)
                                    {
                                        h_mc1d[i][VVV]->Add(h_mc1d[i][j]);
                                        zeroHist(h_mc1d[i][j]);
                                        if (leptype == 1)
                                        {
                                            h_mc1d_comb[i][VVV]->Add(h_mc1d_comb[i][j]);
                                            zeroHist(h_mc1d_comb[i][j]);
                                        }
                                    }

                                    if (j == TWSL)
                                    {
                                        h_mc1d[i][TTSL]->Add(h_mc1d[i][j]);
                                        zeroHist(h_mc1d[i][j]);
                                        if (leptype == 1)
                                        {
                                            h_mc1d_comb[i][TTSL]->Add(h_mc1d_comb[i][j]);
                                            zeroHist(h_mc1d_comb[i][j]);
                                        }
                                    }
                */
                sorted_mc1d[i].push_back( h_mc1d[i][j] );
                if (leptype == nCh - 1) sorted_mc1d_comb[i].push_back( h_mc1d_comb[i][j] );
                if (leptype == nCh - 1) unsorted_mc1d_comb[i].push_back( h_mc1d_comb[i][j] );

            }
            //cout << "sorting" << endl;
            //sort mc histograms
            sort( sorted_mc1d[i].begin(), sorted_mc1d[i].end(), sortByIntegral );
            sort( sorted_mc1d_comb[i].begin(), sorted_mc1d_comb[i].end(), sortByIntegral );
        }

        // End of open and book histograms ////////////////////////////////////////

        // Style up histograms ////////////////////////////////////////////////////

        cout << "Styling 1D" << endl;
        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (!hasplot[leptype][i]) continue;

            h_dt1d[i]->SetLineColor(kBlack);
            h_dt1d[i]->SetMarkerColor(kBlack);
            h_dt1d[i]->SetFillColor(0);
            h_dt1d[i]->SetFillStyle(0);
            for (int j = 0; j < MCID; ++j)
            {
                h_mc1d[i][j]->SetLineColor(kBlack);
                h_mc1d[i][j]->SetMarkerColor(mccolor[j]);
                h_mc1d[i][j]->SetFillColor(mccolor[j]);
                //      h_mc1d[i][j]->SetFillColor(j!=TTDL?mccolor[j]:10);
                h_mc1d[i][j]->SetFillStyle(1001);
            }

        }
        // End of style up histograms /////////////////////////////////////////////

        // Add titles to axes /////////////////////////////////////////////////////

        cout << "Setting Titles" << endl;
        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (!hasplot[leptype][i]) continue;

            h_dt1d[i]->GetXaxis()->SetTitle(Form("%s",
                                                 xtitle1d[i]));
            h_dt1d[i]->GetXaxis()->SetTitleOffset(1.);
            h_dt1d[i]->GetYaxis()->SetTitle(Form("Entries / %3.1f %s",
                                                 h_dt1d[i]->GetBinWidth(1),
                                                 ytitle1d[i]));

            if (find(logScale.begin(), logScale.end(), i) != logScale.end())
            {
                h_dt1d[i]->GetYaxis()->SetTitleOffset(1.1);
            }
            else
            {
                h_dt1d[i]->GetYaxis()->SetTitleOffset(1.1);
            }

            for (int j = 0; j < MCID; ++j)
            {

                h_mc1d[i][j]->GetXaxis()->SetTitle(Form("%s",
                                                        xtitle1d[i]));
                h_mc1d[i][j]->GetXaxis()->SetTitleOffset(1.);
                h_mc1d[i][j]->GetYaxis()->SetTitle(Form("Entries / %3.1f %s",
                                                        h_dt1d[i]->GetBinWidth(1),
                                                        ytitle1d[i]));

                if (find(logScale.begin(), logScale.end(), i) != logScale.end())
                {
                    h_mc1d[i][j]->GetYaxis()->SetTitleOffset(1.1);
                }
                else
                {
                    h_mc1d[i][j]->GetYaxis()->SetTitleOffset(1.1);
                }

            }
        }
        // End of add titles to axes /////////////////////////////////////////////

        // Print expected number of MC events, after all normalizations
        /*
        for (int j = 0; j < MCID; ++j)
        {
            cout << endl << "Expected " << legend[j]
                 << ": " << h_mc1d[1][j]->Integral(1, h_mc1d[1][j]->GetNbinsX())
                 << endl;
        }
        cout << endl << "Expected Data " << h_dt1d[1]->Integral(1, h_dt1d[1]->GetNbinsX())
             << endl;
        */

        // Stack histograms ////////////////////////////////////////////////////
        THStack *s_mc1d[N1DHISTS];

        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (!hasplot[leptype][i]) continue;

            s_mc1d[i] = new THStack(file1dname[i], file1dname[i]);
            s_mc1d[i]->SetName(Form("stack_%s%s", histtag[histtag_number], file1dname[i]));

            for (int j = 0; j < MCID; ++j)
            {
                s_mc1d[i]->Add(sorted_mc1d[i].at(j));
                //    for (int j=MCID-1;j>=0;--j) {
                //      s_mc1d[i]->Add(h_mc1d[i][j]);
            }
            //cout << endl << "MC " << h_mc1d_tot[i]->Integral(1, h_mc1d_tot[i]->GetNbinsX())
            //     << endl;
            //cout << endl << "Data " << h_dt1d[1]->Integral(1, h_dt1d[1]->GetNbinsX())
            //     << endl;
            //cout << endl;
            // if (i==0) {
            //   TFile *outfile_sl=new TFile("mt_leadmuo_njge4.root","RECREATE");
            //   s_mc1d[i]->Write();
            //   h_dt1d[i]->Write();
            //   outfile_sl->Write();
            //   outfile_sl->Close();
            // }


        }
        // End of stack histograms /////////////////////////////////////////////

        // Legends!
        cout << "Doing Legends" << endl;
        TLegend *leg1d[N1DHISTS];

        // Will make them in same position, put "if(i==..)" to change them
        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (!hasplot[leptype][i]) continue;

            leg1d[i] = new TLegend(0.57, 0.60, 0.93, 0.931);//0.56,0.64,0.92,0.915);
            //  leg1d[i] = new TLegend(0.712, 0.541, 0.928, 0.931);//0.56,0.64,0.92,0.915);
            //    leg1d[i] = new TLegend(0.731544, 0.55507, 0.947987, 0.946241);//0.56,0.64,0.92,0.915);
            leg1d[i]->SetName(Form("leg_%s", h_dt1d[i]->GetName()));
            //    leg1d[i]->SetTextSize(0.03);
            leg1d[i]->SetFillColor(0);
            leg1d[i]->
            AddEntry(h_dt1d[i],
                     "data ", "lp");

            for (int j = 0; j < MCID; ++j)
            {
                //if (j == TWSL) continue;
                //if (j >= TWDL && j < VVV) continue;
                leg1d[i]->AddEntry(h_mc1d[i][j],
                                   legend[j], "f");
            }
        }

        // Canvases, at last!
        TCanvas *canv1d[N1DHISTS];

        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (!hasplot[leptype][i]) continue;

            canv1d[i] = new TCanvas(Form("c_%s_%s", file1dname[i], leptag[leptype]),
                                    Form("c_%s_%s", file1dname[i], leptag[leptype]),
                                    700, 700);
            canv1d[i]->SetLeftMargin(0.18);
            canv1d[i]->SetRightMargin(0.05);

            TPad *fullpad = new TPad();
            TPad *plotpad = new TPad();
            TPad *respad  = new TPad();

            fullpad = new TPad("fullpad", "fullpad", 0, 0, 1, 1);
            fullpad->Draw();
            fullpad->cd();

            plotpad = new TPad("plotpad", "plotpad", 0, 0, 1, 0.8);
            plotpad->Draw();
            plotpad->cd();

            if (find(logScale.begin(), logScale.end(), i) != logScale.end())
            {
                plotpad->SetLogy();
                float maxval = h_dt1d[i]->GetMaximum() * 300.;
                h_dt1d[i]->SetMaximum(maxval);
                h_dt1d[i]->SetMinimum(5e-2);
                s_mc1d[i]->SetMinimum(5e-2);
                // h_dt1d[i]->SetMinimum(5e-3);
                // s_mc1d[i]->SetMinimum(5e-3);
            }
            else
            {
                if (i == 5)
                {
                    h_dt1d[i]->SetMaximum(h_dt1d[i]->GetMaximum() * 2.);
                    s_mc1d[i]->SetMaximum(h_dt1d[i]->GetMaximum() * 2.);
                }
                else
                {
                    h_dt1d[i]->SetMaximum(h_dt1d[i]->GetMaximum() * 2.);
                    s_mc1d[i]->SetMaximum(h_dt1d[i]->GetMaximum() * 2.);
                }

                s_mc1d[i]->SetMinimum(0.);
                h_dt1d[i]->SetMinimum(0.);
            }

            TH1F *h_basic = (TH1F *)h_dt1d[i]->Clone();
            h_basic->SetName("h_basic");
            h_basic->SetMarkerColor(kWhite);
            h_basic->SetMarkerSize(0.000001);
            h_basic->SetMarkerStyle(1);
            h_basic->Draw("e3");
            s_mc1d[i]->Draw("hist,same");

            TGraphAsymmErrors *g_data = GetPoissonizedGraph(h_dt1d[i]);
            // g_data->SetMarkerStyle(20);
            // g_data->SetMarkerSize(0.9);
            // g_data->SetLineWidth(2);
            g_data->Draw("e1,p,z,same");
            // h_dt1d[i]->SetMarkerStyle(20);
            // h_dt1d[i]->SetMarkerSize(0.9);
            // h_dt1d[i]->Draw("p,same");
            h_dt1d[i]->Draw("axis,same");

            leg1d[i]->Draw("same");
            /*
              h_dt1d[i]->GetXaxis()->SetTitle(xtitle1d[i]);
              h_dt1d[i]->Draw("E1");
              s_mc1d[i]->Draw("samehist");
              h_dt1d[i]->Draw("E1SAME");
              h_dt1d[i]->Draw("E1AXISSAME");

              // s_mc1d[i]->Draw();
              // s_mc1d[i]->Draw("hist");

              // TGraphAsymmErrors* g_data = GetPoissonizedGraph(h_dt1d[i]);
              // g_data->SetMarkerStyle(20);
              // g_data->SetMarkerSize(0.9);
              // g_data->SetLineWidth(2);
              // g_data->Draw("z,same");
              // h_dt1d[i]->SetMarkerStyle(20);
              // h_dt1d[i]->SetMarkerSize(0.9);
              // h_dt1d[i]->Draw("p,same");
              // h_dt1d[i]->Draw("axis,same");

              leg1d[i]->Draw("same");
            */
            // Labels!
            //    myCMSLabel(0.2,0.88,CMSLabel);

            double chi2prob = compatibilityTest(h_dt1d[i], s_mc1d[i]);
            double KSprob = h_dt1d[i]->KolmogorovTest(h_mc1d_tot[i]);
            chi2sum_channels[leptype] += chi2prob; KSsum_channels[leptype] += KSprob;

            double asymdata = GetAfb(h_dt1d[i]);
            double asymMC = GetAfb(h_mc1d_tot[i]);

            TLatex *text = new TLatex();
            text->SetNDC();
            text->SetTextSize(0.04);
            float xtex = 0.2;//0.16;//used to be 0.2
            text->DrawLatex(xtex, 0.88, "CMS Preliminary");
            text->DrawLatex(xtex, 0.83, Form("#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = %.1f fb^{-1}", lumi));
            if(!printnumbersonplots) text->DrawLatex(xtex, 0.78, Form("%s", leplabel[leptype]));
            else if(scalettdiltodata) text->DrawLatex(xtex, 0.78, Form("%s, #chi^{2} prob: %.2f, KS: %.2f",leplabel[leptype], chi2prob, KSprob));
            else text->DrawLatex(xtex, 0.78, Form("%s, KS: %.2f",leplabel[leptype], KSprob));
            if(i<8 && printnumbersonplots) text->DrawLatex(xtex, 0.73, Form("A_{data}=%.2f, A_{MC}=%.2f", asymdata, asymMC));

            if (i == 3)
            {
                //text->DrawLatex(xtex, 0.68, Form("MC = %.0f #pm %.0f, Data = %.0f", mc_mtctail[isr][leptype], emc_mtctail[isr][leptype], dt_mtctail[isr][leptype]));
                //text->DrawLatex(xtex, 0.63, Form("SF = %.2f #pm %.2f", sf_dtmc_mtctail[isr][leptype], esf_dtmc_mtctail[isr][leptype]));
                // text->DrawLatex(0.68,0.45,Form("MC = %.0f #pm %.0f", mc_mtctail[isr][leptype],emc_mtctail[isr][leptype]));
                // text->DrawLatex(0.68,0.40,Form("Data = %.0f",dt_mtctail[isr][leptype]));
                // text->DrawLatex(0.68,0.35,Form("SF = %.2f #pm %.2f",
                //                     sf_dtmc_mtctail[isr][leptype],esf_dtmc_mtctail[isr][leptype]));
            }
            // if (i==3) {//mt distribution
            //   text->DrawLatex(xtex,0.73,Form("m_{T}>160: MC = %.1f #pm %.1f, Data = %.0f",
            //                   mc_mttail,emc_mttail,dt_mttail));
            //   text->DrawLatex(xtex,0.68,Form("SF = %.2f #pm %.2f",
            //                   sf_dtmc_mttail,esf_dtmc_mttail));
            // } else if (i==6) {//pt distribution
            //   text->DrawLatex(xtex,0.73,Form("p_{T}>150: MC = %.1f #pm %.1f, Data = %.0f",
            //                   mc_lpttail,emc_lpttail,dt_lpttail));
            //   text->DrawLatex(xtex,0.68,Form("SF = %.2f #pm %.2f",
            //                   sf_dtmc_lpttail,esf_dtmc_lpttail));
            // }

            fullpad->cd();

            respad = new TPad("respad", "respad", 0, 0.8, 1, 1);
            respad->Draw();
            respad->cd();

            //gPad->SetGridy();
            //cout << "----------------------------------------------" << endl;
            TH1F *ratio = (TH1F *) h_dt1d[i]->Clone(Form("%s_ratio", h_dt1d[i]->GetName()));
            ratio->Divide(h_mc1d_tot[i]);

            ratio->GetYaxis()->SetTitleOffset(0.3);
            ratio->GetYaxis()->SetTitleSize(0.2);
            ratio->GetYaxis()->SetNdivisions(5);
            ratio->GetYaxis()->SetLabelSize(0.2);
            //    if (dozoom) ratio->GetYaxis()->SetRangeUser(0.8,1.2);
            //    if (dozoom) ratio->GetYaxis()->SetRangeUser(0.7,1.3);
            //    else
            //    ratio->GetYaxis()->SetRangeUser(0.7,1.3);
            ratio->GetYaxis()->SetRangeUser(0., 2.);
            //if (i == 3) ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
            //ratio->GetYaxis()->SetRangeUser(0.5,1.5);
            ratio->GetYaxis()->SetTitle("data/SM  ");
            ratio->GetXaxis()->SetLabelSize(0);
            ratio->GetXaxis()->SetTitleSize(0);
            //  ratio->SetMarkerSize(1);
            ratio->Draw();

            // float xmin = ratio->GetXaxis()->GetBinLowEdge(1);
            // float xmax = ratio->GetXaxis()->GetBinUpEdge(ratio->GetNbinsX());

            // TF1 *f_line = new TF1("f_line", straightline, xmin, xmax, 2);
            // ratio->Fit("f_line");
            // f_line->SetLineColor(kRed);
            // f_line->SetLineWidth(2);
            // f_line->Draw("SAME");

            TLine line;
            line.SetLineWidth(1);
            line.DrawLine(h_dt1d[i]->GetXaxis()->GetXmin(), 1, h_dt1d[i]->GetXaxis()->GetXmax(), 1);

            gErrorIgnoreLevel = 2000; //suppress Info in <TCanvas::Print> messages
            canv1d[i]->Print(Form("%splots/%s%s%s%s.pdf", region,
                                  file1dname[i], metcut[isr],
                                  leptag[leptype],
                                  ""));
                                  //ttbar_tag));
            gErrorIgnoreLevel = 0;

            //cout << " Probability " << file1dname[i] << " : "
                 //<< chi2prob <<" ";
            //cout<< h_dt1d[i]->Chi2Test(h_mc1d_tot[i],"UW")<<" "<< KSprob<<endl;

            //delete g_data;
            //delete canv1d[i];
        }
    /*
        cout << "-------------------------------------------------" << endl;
        cout << "************"<<leplabel[leptype]<<" YIELDS****************************" << endl;
        cout << "-------------------------------------------------" << endl;
        double error_temp = 0.;
        //printf("%s MC \t\t \n", leptag[leptype]);
        for (int j = 0; j < MCID; ++j)
        {
            printf(" %s : %.1f ± ", legend[j], h_mc1d[0][j]->IntegralAndError(0,h_mc1d[0][j]->GetNbinsX()+1,error_temp));
            printf("%.1f \n", error_temp);
        }
        printf(" Total_MC : %.1f ± ", h_mc1d_tot[0]->IntegralAndError(0,h_mc1d_tot[0]->GetNbinsX()+1,error_temp));
        printf("%.1f \n", error_temp);
        printf(" Data : %.0f \n", h_dt1d[0]->Integral());
        cout << "-------------------------------------------------" << endl;
        cout << "*************************************************" << endl;
        cout << "-------------------------------------------------" << endl;
*/

    }//end loop over lepton types



    //ee+mumu+emu combined plots
    // Style up histograms ////////////////////////////////////////////////////

    cout << "Styling 1D" << endl;
    for (int i = 0; i < N1DHISTS; ++i)
    {
        if (!hasplotall[i]) continue;

        h_dt1d_comb[i]->SetLineColor(kBlack);
        h_dt1d_comb[i]->SetMarkerColor(kBlack);
        h_dt1d_comb[i]->SetFillColor(0);
        h_dt1d_comb[i]->SetFillStyle(0);
        for (int j = 0; j < MCID; ++j)
        {
            h_mc1d_comb[i][j]->SetLineColor(kBlack);
            h_mc1d_comb[i][j]->SetMarkerColor(mccolor[j]);
            h_mc1d_comb[i][j]->SetFillColor(mccolor[j]);
            //      h_mc1d_comb[i][j]->SetFillColor(j!=TTDL?mccolor[j]:10);
            h_mc1d_comb[i][j]->SetFillStyle(1001);
        }
    }
    // End of style up histograms /////////////////////////////////////////////

    // Add titles to axes /////////////////////////////////////////////////////

    cout << "Setting Titles" << endl;
    for (int i = 0; i < N1DHISTS; ++i)
    {
        if (!hasplotall[i]) continue;

        h_dt1d_comb[i]->GetXaxis()->SetTitle(Form("%s",
                                             xtitle1d[i]));
        h_dt1d_comb[i]->GetYaxis()->SetTitle(Form("Entries / %3.0f %s",
                                             h_dt1d_comb[i]->GetBinWidth(1),
                                             ytitle1d[i]));

        h_dt1d_comb[i]->GetYaxis()->SetTitleOffset(1.1);
    }
    // End of add titles to axes /////////////////////////////////////////////

    // Stack histograms ////////////////////////////////////////////////////
    THStack *s_mc1d_comb[N1DHISTS];

    for (int i = 0; i < N1DHISTS; ++i)
    {
        if (!hasplotall[i]) continue;

        s_mc1d_comb[i] = new THStack(file1dname[i], file1dname[i]);
        s_mc1d_comb[i]->SetName(Form("stack_%s%s", histtag[histtag_number], file1dname[i]));


        for (int j = 0; j < MCID; ++j)
        {
            s_mc1d_comb[i]->Add(sorted_mc1d_comb[i].at(j));
        }


        //s_mc1d_comb[i]->Add(sorted_mc1d_comb[i].at(VVV));
        //s_mc1d_comb[i]->Add(sorted_mc1d_comb[i].at(WJETS));
        //s_mc1d_comb[i]->Add(sorted_mc1d_comb[i].at(TTDL));
        //s_mc1d_comb[i]->Add(sorted_mc1d_comb[i].at(TTSL));

        /*
                TList *histos = s_mc1d_comb[i]->GetHists();
                double stacksum = 0;
                TIter next(histos);
                TH1F *hist;
                while ((hist = (TH1F *)next()))
                {
                    cout << "Adding " << hist->GetName() << endl;
                    stacksum += hist[i].Integral(1, hist->GetNbinsX());
                }

                cout << endl << "stack " << stacksum << endl;
        */

        //cout << endl << "MC " << h_mc1d_tot_comb[i]->Integral(1, h_mc1d_tot_comb[i]->GetNbinsX())
        //     << endl;
        //cout << endl << "Data " << h_dt1d_comb[i]->Integral(1, h_dt1d_comb[i]->GetNbinsX())
        //     << endl;
        //cout << endl;
    }
    // End of stack histograms /////////////////////////////////////////////

    // Legends!
    cout << "Doing Legends" << endl;
    TLegend *leg1d_comb[N1DHISTS];

    // Will make them in same position, put "if(i==..)" to change them
    for (int i = 0; i < N1DHISTS; ++i)
    {
        if (!hasplotall[i]) continue;

        leg1d_comb[i] = new TLegend(0.57, 0.60, 0.93, 0.931);
        leg1d_comb[i]->SetName(Form("leg_%s", h_dt1d_comb[i]->GetName()));
        leg1d_comb[i]->SetFillColor(0);
        leg1d_comb[i]->AddEntry(h_dt1d_comb[i], "data ", "lp");
        //leg1d_comb[i]->AddEntry(h_mc1d_comb[i][TTSL], legend[TTSL], "f");
        //leg1d_comb[i]->AddEntry(h_mc1d_comb[i][TTDL], legend[TTDL], "f");
        //leg1d_comb[i]->AddEntry(h_mc1d_comb[i][WJETS], legend[WJETS], "f");
        //leg1d_comb[i]->AddEntry(h_mc1d_comb[i][VVV], legend[VVV], "f");
        for (int j = 0; j < MCID; ++j)
        {
            leg1d_comb[i]->AddEntry(h_mc1d_comb[i][j], legend[j], "f");
        }

    }

    // Canvases, at last!
    TCanvas *canv1d_comb[N1DHISTS];

    for (int i = 0; i < N1DHISTS; ++i)
    {
        if (!hasplotall[i]) continue;

        canv1d_comb[i] = new TCanvas(Form("c_%s", file1dname[i]),
                                     Form("c_%s", file1dname[i]),
                                     700, 700);
        canv1d_comb[i]->SetLeftMargin(0.18);
        canv1d_comb[i]->SetRightMargin(0.05);

        TPad *fullpad = new TPad();
        TPad *plotpad = new TPad();
        TPad *respad  = new TPad();

        fullpad = new TPad("fullpad", "fullpad", 0, 0, 1, 1);
        fullpad->Draw();
        fullpad->cd();

        plotpad = new TPad("plotpad", "plotpad", 0, 0, 1, 0.8);
        plotpad->Draw();
        plotpad->cd();

        if (find(logScale.begin(), logScale.end(), i) != logScale.end())
        {
            plotpad->SetLogy();
            float maxval = h_dt1d_comb[i]->GetMaximum() * 300.;
            h_dt1d_comb[i]->SetMaximum(maxval);
            h_dt1d_comb[i]->SetMinimum(5e-2);
            s_mc1d_comb[i]->SetMinimum(5e-2);
        }
        else
        {
            if (i == 5)
            {
                h_dt1d_comb[i]->SetMaximum(h_dt1d_comb[i]->GetMaximum() * 2.);
                s_mc1d_comb[i]->SetMaximum(h_dt1d_comb[i]->GetMaximum() * 2.);
            }
            else
            {
                h_dt1d_comb[i]->SetMaximum(h_dt1d_comb[i]->GetMaximum() * 2.);
                s_mc1d_comb[i]->SetMaximum(h_dt1d_comb[i]->GetMaximum() * 2.);
            }
            s_mc1d_comb[i]->SetMinimum(0.);
            h_dt1d_comb[i]->SetMinimum(0.);
        }

        //cout << "data entries " << h_dt1d_comb[i]->Integral() << endl;
        //    h_dt1d_comb[i]->Draw();

        TH1F *h_basic = (TH1F *)h_dt1d_comb[i]->Clone();
        h_basic->SetName("h_basic");
        h_basic->SetMarkerColor(kWhite);
        h_basic->SetMarkerSize(0.000001);
        h_basic->SetMarkerStyle(1);
        h_basic->Draw("e3");
        s_mc1d_comb[i]->Draw("hist,same");

        TGraphAsymmErrors *g_data = GetPoissonizedGraph(h_dt1d_comb[i]);
        g_data->Draw("e1,p,z,same");
        h_dt1d_comb[i]->Draw("axis,same");
        leg1d_comb[i]->Draw("same");

        //if(i==57) h_dt1d_comb[i]->Print("all");
        //if(i==57) s_mc1d_comb[i]->Print("all");
        //if(i==57) h_mc1d_tot_comb[i]->Print("all");
        double chi2prob = compatibilityTest(h_dt1d_comb[i], s_mc1d_comb[i]);
        double KSprob = h_dt1d_comb[i]->KolmogorovTest(h_mc1d_tot_comb[i]);
        chi2sum += chi2prob; KSsum += KSprob;

        double asymdata = GetAfb(h_dt1d_comb[i]);
        double asymMC = GetAfb(h_mc1d_tot_comb[i]);

        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.04);
        float xtex = 0.2;//0.16;//used to be 0.2
        text->DrawLatex(xtex, 0.88, "CMS Preliminary");
        text->DrawLatex(xtex, 0.83, Form("#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = %.1f fb^{-1}", lumi));
        if(scalettdiltodata && printnumbersonplots) text->DrawLatex(xtex, 0.78, Form("#chi^{2} prob: %.2f, KS: %.2f", chi2prob, KSprob));
        else if(printnumbersonplots) text->DrawLatex(xtex, 0.78, Form("KS: %.2f", KSprob));
        if(i<8 && printnumbersonplots) text->DrawLatex(xtex, 0.73, Form("A_{data}=%.2f, A_{MC}=%.2f", asymdata, asymMC));

        fullpad->cd();

        respad = new TPad("respad", "respad", 0, 0.8, 1, 1);
        respad->Draw();
        respad->cd();

        //gPad->SetGridy();
        //cout << "----------------------------------------------" << endl;
        TH1F *ratio = (TH1F *) h_dt1d_comb[i]->Clone(Form("%s_ratio", h_dt1d_comb[i]->GetName()));
        ratio->Divide(h_mc1d_tot_comb[i]);

        ratio->GetYaxis()->SetTitleOffset(0.3);
        ratio->GetYaxis()->SetTitleSize(0.2);
        ratio->GetYaxis()->SetNdivisions(5);
        ratio->GetYaxis()->SetLabelSize(0.2);
        ratio->GetYaxis()->SetRangeUser(0., 2.);
        ratio->GetYaxis()->SetTitle("data/SM  ");
        ratio->GetXaxis()->SetLabelSize(0);
        ratio->GetXaxis()->SetTitleSize(0);
        //      ratio->SetMarkerSize(1);
        ratio->Draw();

        TLine line;
        line.SetLineWidth(1);
        line.DrawLine(h_dt1d_comb[i]->GetXaxis()->GetXmin(), 1, h_dt1d_comb[i]->GetXaxis()->GetXmax(), 1);

        gErrorIgnoreLevel = 2000; //suppress Info in <TCanvas::Print> messages
        canv1d_comb[i]->Print(Form("%splots/%s%s%s_combined.pdf", region,
                                   file1dname[i], metcut[isr],
                                   ""));
                                   //ttbar_tag));
        gErrorIgnoreLevel = 0;

        /*
        cout << " Probability " << file1dname[i] << " : "
             //<< compatibilityTest(h_dt1d_comb[i], s_mc1d_comb[i]) <<" ";
             << chi2prob <<" ";
        cout<< h_dt1d_comb[i]->Chi2Test(h_mc1d_tot_comb[i],"UW")<<" "<< KSprob <<endl;
        */

        //delete g_data;

    }
    //}  NSAMPLE loop



/*
    cout << "-------------------------------------------------" << endl;
    cout << "**********************YIELDS*********************" << endl;
    cout << "-------------------------------------------------" << endl;
    double error_temp = 0.;
    //printf("%s MC \t\t \n", leptag[leptype]);
    for (int j = 0; j < MCID; ++j)
    {
        printf(" %s : %.1f ± ", legend[j], h_mc1d_comb[0][j]->IntegralAndError(0,h_mc1d_comb[0][j]->GetNbinsX()+1,error_temp));
        printf("%.1f \n", error_temp);
    }
    printf(" Total_MC : %.1f ± ", h_mc1d_tot_comb[0]->IntegralAndError(0,h_mc1d_tot_comb[0]->GetNbinsX()+1,error_temp));
    printf("%.1f \n", error_temp);
    printf(" Data : %.0f \n", h_dt1d_comb[0]->Integral());
    cout << "-------------------------------------------------" << endl;
    cout << "*************************************************" << endl;
    cout << "-------------------------------------------------" << endl;
*/


    cout << "-------------------------------------------------" << endl;
    cout << "**********************ASYMS**********************" << endl;
    cout << "-------------------------------------------------" << endl;

    float asymmetry = 0.;
    float asymmetryerror = 0.;

for (int k = 0; k < 8; ++k)
{
    if( strncmp(file1dname[k],"lep_azimuthal_asymmetry", 1000 ) == 0 ){
        printf(" Variable: %s \n", file1dname[k]);
        for (int j = 0; j < MCID; ++j)
        {
            printf("  %s : %.4f %.4f %.4f \n", mcsample[j], GetAfb(h_mc1d_comb[k][j]), GetAfbLHalf(h_mc1d_comb[k][j]), GetAfbRHalf(h_mc1d_comb[k][j]) );
        }
        printf(" Total_MC : %.4f %.4f %.4f \n", GetAfb(h_mc1d_tot_comb[k]), GetAfbLHalf(h_mc1d_tot_comb[k]), GetAfbRHalf(h_mc1d_tot_comb[k]) );
        printf(" Data : %.4f %.4f %.4f \n", GetAfb(h_dt1d_comb[k]), GetAfbLHalf(h_dt1d_comb[k]), GetAfbRHalf(h_dt1d_comb[k]) );
        cout << " -------------------------------------------------" << endl;
        cout << " *************************************************" << endl;
        cout << " -------------------------------------------------" << endl;
    }
    else {
        printf(" Variable: %s \n", file1dname[k]);
        for (int j = 0; j < MCID; ++j)
        {
            GetAfberr(h_mc1d_comb[k][j], asymmetry, asymmetryerror);
            printf("  %s : %.4f \\pm %.4f \n", mcsample[j], asymmetry, asymmetryerror);
        }
        GetAfberr(h_mc1d_tot_comb[k], asymmetry, asymmetryerror);
        printf(" Total_MC : %.4f \\pm %.4f \n", asymmetry, asymmetryerror);
        GetAfberr(h_dt1d_comb[k], asymmetry, asymmetryerror);
        printf(" Data : %.4f \\pm %.4f \n", asymmetry, asymmetryerror);
        cout << " -------------------------------------------------" << endl;
        cout << " *************************************************" << endl;
        cout << " -------------------------------------------------" << endl;
    }
}

    cout << "-------------------------------------------------" << endl;
    cout << "******************Compatability******************" << endl;
    cout << "-------------------------------------------------" << endl;
    for (int leptype = 0; leptype < nCh; ++leptype) {
        printf(" %s, mean chi2 prob: %.3f, mean KS prob: %.3f \n", leplabel[leptype], chi2sum_channels[leptype]/N1DHISTS, KSsum_channels[leptype]/N1DHISTS);   
    }
    printf(" combined, mean chi2 prob: %.3f, mean KS prob: %.3f \n", chi2sum/N1DHISTS, KSsum/N1DHISTS);   


    cout << "-------------------------------------------------" << endl;
    cout << "**********************YIELDS*********************" << endl;
    cout << "-------------------------------------------------" << endl;

    cout << "Before ttbar solution: "<<endl;
    printYields(unsorted_mc1d_comb[12], legend, h_dt1d_comb[12], true);

    cout << "After ttbar solution: "<<endl;
    printYields(unsorted_mc1d_comb[13], legend, h_dt1d_comb[13], true);



    dt_dl[DIEL] -> Close();
    dt_dl[DIMU] -> Close();
    dt_dl[MUEG] -> Close();

    for (int j = 0; j < MCID; ++j)
    {
            mc_dl[j] -> Close();
    }




    //#endif
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors *GetPoissonizedGraph(TH1F *histo)
{

    TGraphAsymmErrors *graph = new TGraphAsymmErrors();
    graph->SetName(Form("g_%s", histo->GetName()));

    int j = 0;
    for (int i = 1; i <= histo->GetNbinsX(); ++i)
    {

        if (histo->GetBinContent(i) != 0)
        {

            graph->SetPoint(j, histo->GetBinCenter(i), histo->GetBinContent(i));

            graph->SetPointError(j,
                                 histo->GetBinWidth(i) / 2.,
                                 histo->GetBinWidth(i) / 2.,
                                 error_poisson_down(histo->GetBinContent(i)),
                                 error_poisson_up(histo->GetBinContent(i)));

            ++j;
        }
    }
    return graph;

}

double error_poisson_up(double data)
{
    double y1 = data + 1.0;
    double d = 1.0 - 1.0 / (9.0 * y1) + 1.0 / (3 * TMath::Sqrt(y1));
    return y1 * d * d * d - data;
}

double error_poisson_down(double data)
{
    double y = data;
    if (y == 0.0) return 0.0;
    double d = 1.0 - 1.0 / (9.0 * y) - 1.0 / (3.0 * TMath::Sqrt(y));
    return data - y * d * d * d;
}
//////////////////////////////////////////////////////////////////////////////


// New function for comparing errors

double compatibilityTest(TH1F *hist, THStack *hstack)
{

    // A bit of a disaster: need to loop on list of histograms in stack

    TList *stack_hists = hstack->GetHists();

    int ndf = 0;
    double chisq = 0.;
    double h_entries = 0.; // entries in histo
    double s_entries = 0.; // entries in stack
    double s_sqerror = 0.;
    for (int i = 1; i <= hist->GetNbinsX(); ++i)
    {

        h_entries += hist->GetBinContent(i);
        //TIter iter(stack_hists->MakeIterator());
        TIter iter(stack_hists);
        TH1F *stack_hist;
        //while (TH1F *stack_hist = dynamic_cast<TH1F *>(iter()))
        while ((stack_hist = (TH1F *)iter()))
        {
            s_entries += stack_hist->GetBinContent(i);
            s_sqerror += (stack_hist->GetBinError(i) *
                          stack_hist->GetBinError(i));
        }
        //minimum requirement used to be 20, now lowered to 10
        if (h_entries < 10 && i <= hist->GetNbinsX())
            continue;
        else
        {
            ndf++;

            chisq += (h_entries - s_entries) * (h_entries - s_entries) /
                     (s_sqerror + h_entries);


            h_entries = s_entries = s_sqerror = 0.;
        }

    }
    //cout << "Chisq, ndf: " << chisq << ", " << ndf << endl;
    return TMath::Prob(chisq, ndf);
}


//------------------------------------
// print yield table
//------------------------------------
void printYields( vector<TH1F *> h_mc , const char *labels[] , TH1F *h_data , bool latex )
{

    initSymbols( latex );//for latex plots

    TCanvas *ctemp = new TCanvas();

    printLine(latex);
    printHeader();
    printLine(latex);

    TH1F *hmctot;

    //----------------------
    // print SM MC samples
    //----------------------

    for (unsigned int imc = 0 ; imc < h_mc.size() ; imc++)
    {

        print( h_mc.at(imc) , labels[imc] , false );

        if ( imc == 0 ) hmctot = (TH1F *) h_mc.at(imc)->Clone();
        else           hmctot->Add(h_mc.at(imc));
    }

    printLine(latex);

    //-------------------------------
    // print sum of SM MC samples
    //-------------------------------

    print( hmctot , "total SM MC" );

    printLine(latex);

    print( h_data , "data" );

    printLine(latex);

    delete ctemp;
}

void initSymbols( bool latex )
{

    //-------------------------------------------------------
    // table format
    //-------------------------------------------------------

    width1      = 20;
    width2      = 4;
    extrawidth  = 26;
    linelength  = (width1 + width2) * 5 + 1 + extrawidth;

    //-------------------------------------------------------
    // symbols
    //-------------------------------------------------------

    if ( latex )
    {
        percent    = " \\% ";
        pm         = " $\\pm$ ";
        delim      = "&";
        delimstart = "";
        delimend   = "\\\\";
        ee         = "$ee$";
        mm         = "$\\mu\\mu$";
        em         = "$e\\mu$";
        e          = "$e$";
        m          = "$\\mu$";
    }
    else
    {
        percent    = " % ";
        pm         = " +/- ";
        delim      = "|";
        delimstart = "|";
        delimend   = "|";
        ee         = "ee";
        mm         = "mm";
        em         = "em";
        e          = "e";
        m          = "m";
    }

}

void printLine( bool latex )
{

    if ( latex )
    {
        cout << "\\hline" << endl;
    }
    else
    {
        for ( int i = 0 ; i < linelength ; ++i ) cout << "-";
        cout << endl;
    }
}

void printHeader()
{

    cout << delimstart << setw(width1+extrawidth)   << "Sample"  << setw(width2)
         << delim      << setw(width1) << ee        << setw(width2)
         << delim      << setw(width1) << mm        << setw(width2)
         << delim      << setw(width1) << em        << setw(width2)
         << delim      << setw(width1) << "total"   << setw(width2)
         << delimend   << endl;

}

void print( TH1F *h , string label , bool correlatedError )
{

    stringstream see;
    stringstream smm;
    stringstream sem;
    stringstream stot;

    if ( label == "data" )
    {
        see   << Form( "%.0f" , h->GetBinContent(1) );
        smm   << Form( "%.0f" , h->GetBinContent(2) );
        sem   << Form( "%.0f" , h->GetBinContent(3) );
        stot  << Form( "%.0f" , h->GetBinContent(4) );
    }
    else
    {
        //see  << Form( "%.1f" , h->GetBinContent(1) );
        //smm  << Form( "%.1f" , h->GetBinContent(2) );
        //stot << Form( "%.1f" , h->Integral()       );

        see   << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
        smm   << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
        sem   << Form( "%.1f" , h->GetBinContent(3) ) << pm << Form( "%.1f" , h->GetBinError(3) );
        stot  << Form( "%.1f" , h->GetBinContent(4) ) << pm << Form( "%.1f" , h->GetBinError(4) );

        float error = 0;
        if ( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2) + h->GetBinError(3);
        //else                  error = histError(h, 1, 4);
        else                  error = h->GetBinError(4);

        //    stot << Form( "%.1f" , h->Integral()       ) << pm << Form( "%.1f" , error  );
    }

    cout << delimstart << setw(width1+extrawidth) << label      << setw(width2)
         << delim      << setw(width1) << see.str()   << setw(width2)
         << delim      << setw(width1) << smm.str()   << setw(width2)
         << delim      << setw(width1) << sem.str()   << setw(width2)
         << delim      << setw(width1) << stot.str() << setw(width2)
         << delimend   << endl;


}


void GetAfberr(TH1F* h, Float_t &afb, Float_t  &afberr){

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

float GetAfb(TH1F* h){

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

  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  float afb = (event_plus-event_minus)/(event_plus+event_minus);
  float afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

  //cout<<afb<< "+/-" <<afberr<<endl;

  return afb;

}

float GetAfbLHalf(TH1F* h){

  Int_t nbins = h->GetNbinsX()/2;
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  //event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus  = h-> IntegralAndError(0, nbins/2, event_minus_err,"");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus   = h-> IntegralAndError(nbins/2+1, nbins, event_plus_err,"");
  event_total = event_plus + event_minus;

  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  float afb = (event_plus-event_minus)/(event_plus+event_minus);
  float afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

  return afb;

}

float GetAfbRHalf(TH1F* h){

  Int_t nbins = h->GetNbinsX()/2;
  Float_t event_minus;
  Float_t event_plus;
  Float_t event_total;
  Double_t event_plus_err;
  Double_t event_minus_err;

  //event_minus  = h-> IntegralAndError(0, nbins/2, event_plus_err,"");
  event_minus  = h-> IntegralAndError(nbins+1, nbins + nbins/2, event_minus_err,"");
  //event_plus   = h-> IntegralAndError(nbins/2+1, nbins+1, event_minus_err,"");
  event_plus   = h-> IntegralAndError(nbins + nbins/2+1, nbins + nbins + 1, event_plus_err,"");
  event_total = event_plus + event_minus;

  //cout<<event_minus<<" "<<event_minus_err<<" "<<event_plus<<" "<<event_plus_err<<" "<<event_total<<endl;

  float afb = (event_plus-event_minus)/(event_plus+event_minus);
  float afberr   = sqrt(4*(event_plus*event_plus*event_minus_err*event_minus_err 
    + event_minus*event_minus*event_plus_err*event_plus_err)/
    (event_total*event_total*event_total*event_total));

  return afb;

}