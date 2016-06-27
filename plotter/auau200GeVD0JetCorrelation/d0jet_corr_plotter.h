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

class D0jetCorrPlotter
{
  public:
    D0jetCorrPlotter() {}
    ~D0jetCorrPlotter() {}
    void init(TString inputFileNamu);
    void finish();
    void getPxCut(double);
    void getCorrelation(std::pair<int,int> &, int);
    void getCorrelation(int);
    void plotSignificance(TH2D *);
    void plotCorrelation();
    void plotJetInvM();
    double getSBRatio() {return mSOverC;}
    TH1D *getSignalCorrelation(){return signalCorrelation;};
    TH1D *getBkgCorrelation(){return bkgCorrelation;};
    TH1D *getCandCorrelation(){return candCorrelation;};
    TH1D *getCandCorrelationFar(){return candCorrelationFar;};
    TH1D *getBkgCorrelationFar(){return bkgCorrelationFar;};


  private:
    void fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3]);
    double getSBRatio(TH1D *massHisto);
    TFile *inputFile;
    ofstream mLog;
    TH1D *signalCorrelation;
    TH1D *candCorrelation;
    TH1D *bkgCorrelation;
    TH1D *candCorrelationClose;
    TH1D *bkgCorrelationClose;
    TH1D *candCorrelationFar;
    TH1D *bkgCorrelationFar;
    double mSOverC;
};
