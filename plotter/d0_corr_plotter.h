/***********************************************************************************
 *  *
 *  * d0jet-corr-plotter
 *  *
 *  * Author: Leon He
 *  ***********************************************************************************
 *  *
 *  * Description: 
 *  *
 *  ***********************************************************************************
 *  *
 *  * Log:
 *  *
 *  ***********************************************************************************/
#include "fstream"
#include "iostream"
#include "TString.h"
class TFile;
class TH2F;

class D0CorrPlotter
{
  public:
    D0CorrPlotter() {}
    ~D0CorrPlotter() {}
    void init(TString inputFileNamu);
    void finish();
    void getJetCorrelation(std::pair<int,int> &,std::pair<double,double> &, int,pair<double,double> &);
    void getHadronCorrelation(std::pair<int,int> &,std::pair<double,double> &, int);
    void getHadronJetCorrelation(std::pair<int,int> &,std::pair<double,double> &, int,pair<double,double> &);
    void getCorrelation(int);
    void plotSignificance(TH2D *);
    // void plotCorrelation();
    void plotJetInvM();
    double getSBRatio() {return mSOverC;}
    void setCorrAxis(TH1D *);
    TH1D *getSignalJetCorrelation(){return signalJetCorrelation;};
    TH1D *getBkgJetCorrelation(){return bkgJetCorrelation;};
    TH1D *getCandJetCorrelation(){return candJetCorrelation;};
    TH1D *getSignalHadronCorrelation(){return signalHadronCorrelation;};
    TH1D *getBkgHadronCorrelation(){return bkgHadronCorrelation;};
    TH1D *getCandHadronCorrelation(){return candHadronCorrelation;};
    TH1D *getHadronJetCorrelation(){return hadronJetCorrelation;};


  private:
    void fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3]);
    double getSBRatio(TH1D *massHisto);
    TFile *inputFile;
    ofstream mLog;
    TH1D *signalJetCorrelation;
    TH1D *candJetCorrelation;
    TH1D *bkgJetCorrelation;
    TH1D *signalHadronCorrelation;
    TH1D *candHadronCorrelation;
    TH1D *bkgHadronCorrelation;
    TH1D *hadronJetCorrelation;
    double mSOverC;
};
