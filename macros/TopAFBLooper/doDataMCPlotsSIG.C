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
void print( TH1F *h , string label , bool correlatedError = false );
char *pm;
char *percent;
char *delim;
char *delimstart;
char *delimend;
char *e;
char *m;
int   width1;
int   width2;
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
void doDataMCPlotsSIG(const char *ttbar_tag = "")
{

    //derive scale factors

    //list of samples
    const int MCID = 9;
    const char *mcsample[MCID] =
    {
        "ttdl_powheg",
        "ttsl_powheg",
        "w1to4jets",
        "tW_lepsl",
        "tW_lepdl",
        "diboson",
        "DY1to4Jtot",
        "ttV",
        "triboson"
    };

    enum sample {TTDL = 0,
                 TTSL,
                 WJETS,
                 TWSL,
                 TWDL,
                 DIBO,
                 DY,
                 TTV,
                 VVV
                };

    const char *legend[MCID] =
    {
        "t#bar{t} #rightarrow #font[12]{l^{+}l^{-}}",
        "t#bar{t} #rightarrow #font[12]{l^{#pm}} + jets",
        //    "t#bar{t} #rightarrow #font[12]{l^{#pm}} + jets",
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
        "DY+jets",
        "ttW/Z/#gamma",
        "triboson"
    };
    //    "triboson"};
    //    "QCD"};

    //  const int mccolor[]={7,2,6,4,5,kOrange,9,kAzure-9, 8, kViolet,kGreen+1, 15,12,13,27};
    //  const int mccolor[]={7,2,6,4,kOrange,8,9,kAzure-9, 5, kViolet,kGreen+1, 15,12,13,27};
    const int mccolor[] = {7, 2, 6, 4, kOrange, 8, kViolet + 1, kAzure - 9, kGreen - 2, kViolet, kGreen + 1, 15, 12, 13, 27};
    //  const int mccolor[]={7,2,6,4,kOrange,8,5,kAzure-9, 9,5,kViolet,kGreen+1, 15,12,13,27};

    //-------------------------------
    // SINGLE LEPTON - MT SCALING
    //-------------------------------

    const int nCh = 3;
    char *leptag[nCh] = {"_diel", "_dimu", "_mueg"};
    const char *leplabel[nCh] = {"ee", "#mu#mu", "e#mu"};
    enum channel {DIEL = 0, DIMU, MUEG };
    char *histoname = "h_sig_m_top";


    //initialize

    TFile *dt_dl[nCh];
    dt_dl[DIEL] = TFile::Open("SIGoutput/data_diel_histos.root");
    dt_dl[DIMU] = TFile::Open("SIGoutput/data_dimu_histos.root");
    dt_dl[MUEG] = TFile::Open("SIGoutput/data_mueg_histos.root");

    TFile *mc_dl[MCID];

    for (int j = 0; j < MCID; ++j)
    {
        if (j < 2)
            mc_dl[j] = TFile::Open(Form("SIGoutput/%s%s_histos.root",
                                        mcsample[j], ttbar_tag));
        else
            mc_dl[j] = TFile::Open(Form("SIGoutput/%s_histos.root",
                                        mcsample[j]));
    }

    bool doverbose = true;

    double lumi = 19.5;

    int isr = 0;
    char *metcut[1] = {""};


    //for (int isr = 0; isr < NSAMPLE; ++isr)
    //{

    const int N1DHISTS = 12;
    TH1F *h_dt1d_comb[N1DHISTS];
    TH1F *h_mc1d_comb[N1DHISTS][MCID];
    vector<TH1F *> sorted_mc1d_comb[N1DHISTS];
    TH1F *h_mc1d_tot_comb[N1DHISTS];

    // Output file names
    const char *file1dname[N1DHISTS] =
    {
        "met",
        "lep_azimuthal_asymmetry_2",
        "lep_charge_asymmetry",
        "top_rapiditydiff_Marco",
        "lepPlus_costheta_cms",
        "lepMinus_costheta_cms",
        "top_spin_correlation",
        "lep_cos_opening_angle",
        "m_top",
        "tt_mass",
        "tt_pT",
        "ttRapidity2"
    };

    // Histogram names
    const char *hist1dname[N1DHISTS] =
    {
        "h_sig_met",
        "h_sig_lep_azimuthal_asymmetry_2",
        "h_sig_lep_charge_asymmetry",
        "h_sig_top_rapiditydiff_Marco",
        "h_sig_lepPlus_costheta_cms",
        "h_sig_lepMinus_costheta_cms",
        "h_sig_top_spin_correlation",
        "h_sig_lep_cos_opening_angle",
        "h_sig_m_top",
        "h_sig_tt_mass",
        "h_sig_tt_pT",
        "h_sig_ttRapidity2"
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
    //logScale.push_back(9);
    //logScale.push_back(10);
    //logScale.push_back(11);
    //logScale.push_back(56);

    // List of rebin factors:
    int rebinFactor[N1DHISTS] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    const char *xtitle1d[N1DHISTS] =
    {
        "ME_{T} [GeV]",
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
        "y_{t#bar{t}}"
    };

    const char *ytitle1d[N1DHISTS] =
    {
        "GeV",
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
            h_dt1d[i] = (TH1F *)dt_dl[leptype]->Get(Form("%s%s%s", hist1dname[i], metcut[isr], leptag[leptype]));

            cout << "DT: " << h_dt1d[i] << " " << Form("%s%s", hist1dname[i], leptag[leptype]) << endl;

            h_dt1d[i]->SetName(Form("%s", hist1dname[i]));

            h_dt1d[i]->Rebin(rebinFactor[i]);

            //e+mu combined data
            if (leptype == 0)
            {
                h_dt1d_comb[i] = (TH1F *)h_dt1d[i]->Clone();
                h_dt1d_comb[i]->SetName(Form("%s_comb", hist1dname[i]));
            }
            else
                h_dt1d_comb[i]->Add(h_dt1d[i]);

            bool doinit = false;
            for (int j = 0; j < MCID; ++j)
            {

                h_mc1d[i][j] =
                    (TH1F *)mc_dl[j]->Get(Form("%s%s%s", hist1dname[i], metcut[isr], leptag[leptype]));

                if (h_mc1d[i][j] == 0)
                {
                    h_mc1d[i][j] =
                        (TH1F *)mc_dl[0]->Get(Form("%s%s%s", hist1dname[i], metcut[isr], leptag[0]));
                    h_mc1d[i][j]->SetName(Form("%s_%s", mcsample[j], hist1dname[i]));
                    zeroHist(h_mc1d[i][j]);
                }
                h_mc1d[i][j]->SetName(Form("%s_%s", mcsample[j], hist1dname[i]));

                cout << "MC" << j << " " << mcsample[j] << ": " << h_mc1d[i][j] << " "
                     << Form("%s", hist1dname[i]) << endl;

                h_mc1d[i][j]->Sumw2();

                // Rebin, set limits
                h_mc1d[i][j]->Rebin(rebinFactor[i]);

                //combined MC histograms
                if (leptype == 0)
                {
                    h_mc1d_comb[i][j] = (TH1F *)h_mc1d[i][j]->Clone();
                    h_mc1d_comb[i][j]->SetName(Form("%s_%s_comb", mcsample[j], hist1dname[i]));
                }
                else
                    h_mc1d_comb[i][j]->Add(h_mc1d[i][j]);

                //get total
                if (!doinit && leptype == 0)
                {
                    h_mc1d_tot_comb[i] = (TH1F *)h_mc1d[i][j]->Clone(Form("mctot_comb_%s", hist1dname[i]));
                    h_mc1d_tot_comb[i]->SetName(Form("mctot_comb_%s", hist1dname[i]));
                }
                else
                    h_mc1d_tot_comb[i]->Add(h_mc1d[i][j]);

                if (!doinit)
                {
                    h_mc1d_tot[i] = (TH1F *)h_mc1d[i][j]->Clone(Form("mctot_%s", hist1dname[i]));
                    doinit = true;
                }
                else h_mc1d_tot[i]->Add(h_mc1d[i][j]);

            }
        }


        //after applying all the scalings, fill vector for sorting
        for (int i = 0; i < N1DHISTS; ++i)
        {
            cout << "-------------------------------------------------" << endl;
            float mcall = h_mc1d_tot[i]->Integral();
            float dtall = h_dt1d[i]->Integral();
            cout << "MC ALL " << mcall << endl;
            cout << "Data ALL " << dtall << endl;
            float mcsf_all = dtall / mcall;
            cout << "RATIO ALL " << mcsf_all << endl;
            cout << "-------------------------------------------------" << endl;
            mcsf_all = 1;
            h_mc1d_tot[i]->Scale(mcsf_all);
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
                h_mc1d[i][j]->Scale(mcsf_all);
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

            }
            cout << "sorting" << endl;
            //sort mc histograms
            sort( sorted_mc1d[i].begin(), sorted_mc1d[i].end(), sortByIntegral );
            sort( sorted_mc1d_comb[i].begin(), sorted_mc1d_comb[i].end(), sortByIntegral );
        }

        // End of open and book histograms ////////////////////////////////////////

        // Style up histograms ////////////////////////////////////////////////////

        cout << "Styling 1D" << endl;
        for (int i = 0; i < N1DHISTS; ++i)
        {
            h_dt1d[i]->SetLineColor(kBlack);
            h_dt1d[i]->SetMarkerColor(kBlack);
            h_dt1d[i]->SetFillColor(0);
            h_dt1d[i]->SetFillStyle(0);
            for (int j = 0; j < MCID; ++j)
            {
                h_mc1d[i][j]->SetLineColor(kBlack);
                h_mc1d[i][j]->SetMarkerColor(mccolor[j]);
                h_mc1d[i][j]->SetFillColor(mccolor[j]);
                //      h_mc1d[i][j]->SetFillColor(j>0?mccolor[j]:10);
                h_mc1d[i][j]->SetFillStyle(1001);
            }

        }
        // End of style up histograms /////////////////////////////////////////////

        // Add titles to axes /////////////////////////////////////////////////////

        cout << "Setting Titles" << endl;
        for (int i = 0; i < N1DHISTS; ++i)
        {
            if (i == 3)
            {
                h_dt1d[i]->GetXaxis()->SetBinLabel(1, "Peak");
                h_dt1d[i]->GetXaxis()->SetBinLabel(2, "Tail");
                h_dt1d[i]->GetYaxis()->SetTitle("Entries");
            }
            else
            {
                h_dt1d[i]->GetXaxis()->SetTitle(Form("%s",
                                                     xtitle1d[i]));
                h_dt1d[i]->GetXaxis()->SetTitleOffset(1.);
                h_dt1d[i]->GetYaxis()->SetTitle(Form("Entries / %3.0f %s",
                                                     h_dt1d[i]->GetBinWidth(1),
                                                     ytitle1d[i]));
            }
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

                if (i == 3)
                {
                    h_mc1d[i][j]->GetXaxis()->SetBinLabel(1, "Peak");
                    h_mc1d[i][j]->GetXaxis()->SetBinLabel(2, "Tail");
                    h_mc1d[i][j]->GetYaxis()->SetTitle("Entries");
                }
                else
                {
                    h_mc1d[i][j]->GetXaxis()->SetTitle(Form("%s",
                                                            xtitle1d[i]));
                    h_mc1d[i][j]->GetXaxis()->SetTitleOffset(1.);
                    h_mc1d[i][j]->GetYaxis()->SetTitle(Form("Entries / %3.0f %s",
                                                            h_dt1d[i]->GetBinWidth(1),
                                                            ytitle1d[i]));
                }
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
        for (int j = 0; j < MCID; ++j)
        {
            cout << endl << "Expected " << legend[j]
                 << ": " << h_mc1d[1][j]->Integral(1, h_mc1d[1][j]->GetNbinsX())
                 << endl;
        }
        cout << endl << "Expected Data " << h_dt1d[1]->Integral(1, h_dt1d[1]->GetNbinsX())
             << endl;

        // Stack histograms ////////////////////////////////////////////////////
        THStack *s_mc1d[N1DHISTS];

        for (int i = 0; i < N1DHISTS; ++i)
        {
            s_mc1d[i] = new THStack(hist1dname[i], hist1dname[i]);
            s_mc1d[i]->SetName(Form("stack_%s", hist1dname[i]));

            for (int j = 0; j < MCID; ++j)
            {
                s_mc1d[i]->Add(sorted_mc1d[i].at(j));
                //    for (int j=MCID-1;j>=0;--j) {
                //      s_mc1d[i]->Add(h_mc1d[i][j]);
            }
            cout << endl << "MC " << h_mc1d_tot[i]->Integral(1, h_mc1d_tot[i]->GetNbinsX())
                 << endl;
            cout << endl << "Data " << h_dt1d[1]->Integral(1, h_dt1d[1]->GetNbinsX())
                 << endl;
            cout << endl;
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

            canv1d[i] = new TCanvas(Form("c_%s", file1dname[i]),
                                    Form("c_%s", file1dname[i]),
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
                h_dt1d[i]->SetMinimum(5e-1);
                s_mc1d[i]->SetMinimum(5e-1);
                // h_dt1d[i]->SetMinimum(5e-3);
                // s_mc1d[i]->SetMinimum(5e-3);
            }
            else
            {
                if (i == 5)
                {
                    h_dt1d[i]->SetMaximum(h_dt1d[i]->GetMaximum() * 2.3);
                    s_mc1d[i]->SetMaximum(h_dt1d[i]->GetMaximum() * 2.3);
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
            TLatex *text = new TLatex();
            text->SetNDC();
            text->SetTextSize(0.04);
            float xtex = 0.2;//0.16;//used to be 0.2
            text->DrawLatex(xtex, 0.88, "CMS Preliminary");
            text->DrawLatex(xtex, 0.83, Form("#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = %.1f fb^{-1}", lumi));
            text->DrawLatex(xtex, 0.78, Form("%s", leplabel[leptype]));
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
            cout << "----------------------------------------------" << endl;
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
            if (i == 3) ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
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

            canv1d[i]->Print(Form("SIGplots/%s%s%s%s.pdf",
                                  file1dname[i], metcut[isr],
                                  leptag[leptype],
                                  ttbar_tag));

            cout << " Probability " << file1dname[i] << " : "
                 << compatibilityTest(h_dt1d[i], s_mc1d[i]) << endl;

            //delete canv1d[i];
        }
    }//end loop over lepton types



    //e+mu combined plots
    // Style up histograms ////////////////////////////////////////////////////

    cout << "Styling 1D" << endl;
    for (int i = 0; i < N1DHISTS; ++i)
    {
        h_dt1d_comb[i]->SetLineColor(kBlack);
        h_dt1d_comb[i]->SetMarkerColor(kBlack);
        h_dt1d_comb[i]->SetFillColor(0);
        h_dt1d_comb[i]->SetFillStyle(0);
        for (int j = 0; j < MCID; ++j)
        {
            h_mc1d_comb[i][j]->SetLineColor(kBlack);
            h_mc1d_comb[i][j]->SetMarkerColor(mccolor[j]);
            h_mc1d_comb[i][j]->SetFillColor(mccolor[j]);
            //      h_mc1d_comb[i][j]->SetFillColor(j>0?mccolor[j]:10);
            h_mc1d_comb[i][j]->SetFillStyle(1001);
        }
    }
    // End of style up histograms /////////////////////////////////////////////

    // Add titles to axes /////////////////////////////////////////////////////

    cout << "Setting Titles" << endl;
    for (int i = 0; i < N1DHISTS; ++i)
    {
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
        s_mc1d_comb[i] = new THStack(hist1dname[i], hist1dname[i]);
        s_mc1d_comb[i]->SetName(Form("stack_%s", hist1dname[i]));


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

        cout << endl << "MC " << h_mc1d_tot_comb[i]->Integral(1, h_mc1d_tot_comb[i]->GetNbinsX())
             << endl;
        cout << endl << "Data " << h_dt1d_comb[i]->Integral(1, h_dt1d_comb[i]->GetNbinsX())
             << endl;
        cout << endl;
    }
    // End of stack histograms /////////////////////////////////////////////

    // Legends!
    cout << "Doing Legends" << endl;
    TLegend *leg1d_comb[N1DHISTS];

    // Will make them in same position, put "if(i==..)" to change them
    for (int i = 0; i < N1DHISTS; ++i)
    {
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
            h_dt1d_comb[i]->SetMinimum(5e-1);
            s_mc1d_comb[i]->SetMinimum(5e-1);
            if (i == 7)
            {
                h_dt1d_comb[i]->SetMinimum(5);
                s_mc1d_comb[i]->SetMinimum(5);
            }
            if (isr > 4 || i == 0)
            {
                h_dt1d_comb[i]->SetMinimum(1e-2);
                s_mc1d_comb[i]->SetMinimum(1e-2);
            }
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

        cout << "data entries " << h_dt1d_comb[i]->Integral() << endl;
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
        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.04);
        float xtex = 0.2;//0.16;//used to be 0.2
        text->DrawLatex(xtex, 0.88, "CMS Preliminary");
        text->DrawLatex(xtex, 0.83, Form("#sqrt{s} = 8 TeV, #scale[0.6]{#int}Ldt = %.1f fb^{-1}", lumi));

        fullpad->cd();

        respad = new TPad("respad", "respad", 0, 0.8, 1, 1);
        respad->Draw();
        respad->cd();

        //gPad->SetGridy();
        cout << "----------------------------------------------" << endl;
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

        canv1d_comb[i]->Print(Form("SIGplots/%s%s%s_combined.pdf",
                                   file1dname[i], metcut[isr],
                                   ttbar_tag));

    }
    //}  NSAMPLE loop

    /*
        cout << "-------------------------------------------------" << endl;
        cout << "*************************************************" << endl;
        cout << "-------------------------------------------------" << endl;
        cout << "PRINT SUMMARY YIELDS TABLE" << endl;

        for (int leptype = 0; leptype < nCh; ++leptype)
        {
            cout << "\\hline" << endl;
            cout << "\\hline" << endl;
            printf("%s pre-veto \\mt-SF \t ", leptag[leptype]);
            for (int isr = 0; isr < NSAMPLE; ++isr)
            {
                printf(" & $%.2f \\pm %.2f$", mtsf_pv[isr][leptype], emtsf_pv[isr][leptype]);
            }
            printf(" \\\\\n");
        }
        cout << "\\hline" << endl;
        cout << "\\hline" << endl;
        for (int leptype = 0; leptype < nCh; ++leptype)
        {
            printf("%s post-veto \\mt-SF \t ", leptag[leptype]);
            for (int isr = 0; isr < NSAMPLE; ++isr)
            {
                printf(" & $%.2f \\pm %.2f$", mtsf[isr][leptype], emtsf[isr][leptype]);
            }
            printf(" \\\\\n");
        }
        cout << "\\hline" << endl;
        printf("Comb. post-veto \\mt-SF \t ");
        for (int isr = 0; isr < NSAMPLE; ++isr)
        {
            printf(" & $%.2f \\pm %.2f$", mtsf_comb[isr], emtsf_comb[isr]);
        }
        printf(" \\\\\n");
        cout << "\\hline" << endl;
        cout << "\\hline" << endl;

        cout << "-------------------------------------------------" << endl;
        cout << "*************************************************" << endl;
        cout << "-------------------------------------------------" << endl;
        for (int leptype = 0; leptype < nCh; ++leptype)
        {
            cout << "\\hline" << endl;
            cout << "\\hline" << endl;
            printf("%s MC \t\t ", leptag[leptype]);
            for (int isr = 0; isr < NSAMPLE; ++isr)
            {
                printf(" & $%.0f \\pm %.0f$", mc_mtctail[isr][leptype], emc_mtctail[isr][leptype]);
            }
            printf(" \\\\\n");
            printf("%s Data \t\t ", leptag[leptype]);
            for (int isr = 0; isr < NSAMPLE; ++isr)
            {
                printf(" & $%.0f$", dt_mtctail[isr][leptype]);
            }
            printf(" \\\\\n");
            cout << "\\hline" << endl;
            printf("%s Data/MC SF \t ", leptag[leptype]);
            for (int isr = 0; isr < NSAMPLE; ++isr)
            {
                printf(" & $%.2f \\pm %.2f$", sf_dtmc_mtctail[isr][leptype], esf_dtmc_mtctail[isr][leptype]);
            }
            printf(" \\\\\n");
        }
    */


    //#endif
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors *GetPoissonizedGraph(TH1F *histo)
{

    TGraphAsymmErrors *graph = new TGraphAsymmErrors();
    graph->SetName(Form("g_%s", histo->GetName()));

    int j = 1;
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

        TIter iter(stack_hists->MakeIterator());
        while (TH1F *stack_hist = dynamic_cast<TH1F *>(iter()))
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
    cout << "Chisq, ndf: " << chisq << ", " << ndf << endl;
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

    TH1F *hmctot = new TH1F("hmctot", "hmctot", 2, 0, 2);
    hmctot->Sumw2();

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
    linelength  = (width1 + width2) * 4 + 1;

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

    cout << delimstart << setw(width1) << "Sample"    << setw(width2)
         << delim      << setw(width1) << e           << setw(width2)
         << delim      << setw(width1) << m           << setw(width2)
         << delim      << setw(width1) << "total"     << setw(width2)
         << delimend   << endl;

}

void print( TH1F *h , string label , bool correlatedError )
{

    stringstream se;
    stringstream sm;
    stringstream stot;

    if ( label == "data" )
    {
        se   << Form( "%.0f" , h->GetBinContent(1) );
        sm   << Form( "%.0f" , h->GetBinContent(2) );
        stot << Form( "%.0f" , h->GetBinContent(3) );
    }
    else
    {
        //see  << Form( "%.1f" , h->GetBinContent(1) );
        //smm  << Form( "%.1f" , h->GetBinContent(2) );
        //stot << Form( "%.1f" , h->Integral()       );

        se   << Form( "%.1f" , h->GetBinContent(1) ) << pm << Form( "%.1f" , h->GetBinError(1) );
        sm   << Form( "%.1f" , h->GetBinContent(2) ) << pm << Form( "%.1f" , h->GetBinError(2) );
        stot << Form( "%.1f" , h->GetBinContent(3) ) << pm << Form( "%.1f" , h->GetBinError(3) );

        float error = 0;
        if ( correlatedError ) error = h->GetBinError(1) + h->GetBinError(2) + h->GetBinError(3);
        else                  error = histError(h, 1, 4);

        //    stot << Form( "%.1f" , h->Integral()       ) << pm << Form( "%.1f" , error  );
    }

    cout << delimstart << setw(width1) << label      << setw(width2)
         << delim      << setw(width1) << se.str()   << setw(width2)
         << delim      << setw(width1) << sm.str()   << setw(width2)
         << delim      << setw(width1) << stot.str() << setw(width2)
         << delimend   << endl;


}
