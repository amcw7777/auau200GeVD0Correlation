#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TMath.h"
#include "TAxis.h"
#include "TLegend.h"
#include "Math/MinimizerOptions.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraph.h"
#include "THnSparse.h"

#include "d0_corr_plotter.h"
using namespace std;

void D0CorrPlotter::init(TString inputFileName)
{
	inputFile = new TFile(inputFileName.Data());
  mLog.open("d0_corr_plotter.log");
}

void D0CorrPlotter::finish()
{
}

//
// void D0CorrPlotter::getCorrelation(int rebinFactor)
// {
//   pair<int,int> ptBinCut(2,10);
// 	TH2F *corBkg2D = (TH2F *)inputFile->Get("corBkg")->Clone("corBkg2D");
// 	TH2F *corCand2D = (TH2F *)inputFile->Get("corCand")->Clone("corCand2D");
//   TH2F *candCount = (TH2F *)inputFile->Get("candCount")->Clone("candCount");
//   TH2F *bkgCount = (TH2F *)inputFile->Get("bkgCount")->Clone("bkgCount");
//   double nD0Cand = candCount->ProjectionX("candCount1D")->Integral();
//   double nD0Bkg= bkgCount->ProjectionX("bkgCount1D")->Integral();
//   cout<<"# of candidate D0 = "<<nD0Cand<<endl;
//   cout<<"# of background D0 = "<<nD0Bkg<<endl;
// 	candCorrelation = (TH1D *)corCand2D->ProjectionX("candCorrelation",2,2)->Clone("candCorrelation");
// 	bkgCorrelation = (TH1D *)corBkg2D->ProjectionX("bkgCorrelation",2,2)->Clone("bkgCorrelation");
// 	candCorrelation->Rebin(rebinFactor);
// 	candCorrelation->Scale(1./candCorrelation->GetBinWidth(1)/nD0Cand);
// 	bkgCorrelation->Rebin(rebinFactor);
// 	bkgCorrelation->Scale(1./bkgCorrelation->GetBinWidth(1)/nD0Bkg);
//
// 	TH3F *mass2d = (TH3F *)inputFile->Get("massPt")->Clone("mass2d");
// 	TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptBinCut.first,ptBinCut.second,1,10)->Clone("massAllTrigger");
// 	double f = getSBRatio(massAllTrigger);
// }
void D0CorrPlotter::getJetCorrelation(pair<int,int> &centralityBinCut,pair<double,double> &ptCut,int rebinFactor)
{
	// TH2F *corBkg2D = (TH2F *)inputFile->Get("corBkg")->Clone("corBkg2D");
	// TH2F *corCand2D = (TH2F *)inputFile->Get("corCand")->Clone("corCand2D");
  THnSparseD *d0JetCorrCand = (THnSparseD *)inputFile->Get("hD0JetCorrCand")->Clone("d0JetCorrCand");
  THnSparseD *d0JetCorrBkg = (THnSparseD *)inputFile->Get("hD0JetCorrBkg")->Clone("d0JetCorrBkg");
  TH2F *candCount = (TH2F *)inputFile->Get("candCount")->Clone("candCount");
  TH2F *bkgCount = (TH2F *)inputFile->Get("bkgCount")->Clone("bkgCount");

  double nD0Cand = candCount->ProjectionX("candCount1D",centralityBinCut.first,centralityBinCut.second)->Integral();
  double nD0Bkg= bkgCount->ProjectionX("bkgCount1D",centralityBinCut.first,centralityBinCut.second)->Integral();
  cout<<"# of candidate D0 = "<<nD0Cand<<endl;
  cout<<"# of background D0 = "<<nD0Bkg<<endl;

  d0JetCorrCand->GetAxis(1)->SetRange(centralityBinCut.first,centralityBinCut.second);
  d0JetCorrCand->GetAxis(2)->SetRangeUser(ptCut.first,ptCut.second);
  d0JetCorrBkg->GetAxis(1)->SetRange(centralityBinCut.first,centralityBinCut.second);
  d0JetCorrBkg->GetAxis(2)->SetRangeUser(ptCut.first,ptCut.second);
	// candCorrelation = (TH1D *)corCand2D->ProjectionX("candCorrelation",centralityBinCut.first,centralityBinCut.second)->Clone("candCorrelation");
	// bkgCorrelation = (TH1D *)corBkg2D->ProjectionX("bkgCorrelation",centralityBinCut.first,centralityBinCut.second)->Clone("bkgCorrelation");
	candJetCorrelation = (TH1D *)d0JetCorrCand->Projection(0)->Clone("candJetCorrelation");
	bkgJetCorrelation = (TH1D *)d0JetCorrBkg->Projection(0)->Clone("bkgJetCorrelation");
  setCorrAxis(candJetCorrelation);
  setCorrAxis(bkgJetCorrelation);
  // get candidate correlation
	candJetCorrelation->Rebin(rebinFactor);
	candJetCorrelation->Scale(1./candJetCorrelation->GetBinWidth(1)/nD0Cand);
	bkgJetCorrelation->Rebin(rebinFactor);
	bkgJetCorrelation->Scale(1./bkgJetCorrelation->GetBinWidth(1)/nD0Bkg);

	TH3F *mass2d = (TH3F *)inputFile->Get("massPt")->Clone("mass2d");
	// TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptBinCut.first,ptBinCut.second,centralityBinCut.first,centralityBinCut.second)->Clone("massAllTrigger");
	TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptCut.first,ptCut.second,centralityBinCut.first,centralityBinCut.second)->Clone("massAllTrigger");
	double f = getSBRatio(massAllTrigger);
}
void D0CorrPlotter::getHadronCorrelation(pair<int,int> &centralityBinCut,pair<double,double> &ptCut,int rebinFactor)
{
	// TH2F *corBkg2D = (TH2F *)inputFile->Get("corBkg")->Clone("corBkg2D");
	// TH2F *corCand2D = (TH2F *)inputFile->Get("corCand")->Clone("corCand2D");
  THnSparseD *d0HadronCorrCand = (THnSparseD *)inputFile->Get("hD0HadronCorrCand")->Clone("d0HadronCorrCand");
  THnSparseD *d0HadronCorrBkg = (THnSparseD *)inputFile->Get("hD0HadronCorrBkg")->Clone("d0HadronCorrBkg");
  TH2F *candCount = (TH2F *)inputFile->Get("candCount")->Clone("candCount");
  TH2F *bkgCount = (TH2F *)inputFile->Get("bkgCount")->Clone("bkgCount");

  double nD0Cand = candCount->ProjectionX("candCount1D",centralityBinCut.first,centralityBinCut.second)->Integral();
  double nD0Bkg= bkgCount->ProjectionX("bkgCount1D",centralityBinCut.first,centralityBinCut.second)->Integral();
  cout<<"# of candidate D0 = "<<nD0Cand<<endl;
  cout<<"# of background D0 = "<<nD0Bkg<<endl;

  d0HadronCorrCand->GetAxis(1)->SetRange(centralityBinCut.first,centralityBinCut.second);
  d0HadronCorrCand->GetAxis(2)->SetRangeUser(ptCut.first,ptCut.second);
  d0HadronCorrBkg->GetAxis(1)->SetRange(centralityBinCut.first,centralityBinCut.second);
  d0HadronCorrBkg->GetAxis(2)->SetRangeUser(ptCut.first,ptCut.second);
	// candCorrelation = (TH1D *)corCand2D->ProjectionX("candCorrelation",centralityBinCut.first,centralityBinCut.second)->Clone("candCorrelation");
	// bkgCorrelation = (TH1D *)corBkg2D->ProjectionX("bkgCorrelation",centralityBinCut.first,centralityBinCut.second)->Clone("bkgCorrelation");
	candHadronCorrelation = (TH1D *)d0HadronCorrCand->Projection(0)->Clone("candCorrelation");
	bkgHadronCorrelation = (TH1D *)d0HadronCorrBkg->Projection(0)->Clone("candCorrelation");
  setCorrAxis(candHadronCorrelation);
  setCorrAxis(bkgHadronCorrelation);
  // get candidate correlation
	candHadronCorrelation->Rebin(rebinFactor);
	candHadronCorrelation->Scale(1./candHadronCorrelation->GetBinWidth(1)/nD0Cand);
	bkgHadronCorrelation->Rebin(rebinFactor);
	bkgHadronCorrelation->Scale(1./bkgHadronCorrelation->GetBinWidth(1)/nD0Bkg);

	TH3F *mass2d = (TH3F *)inputFile->Get("massPt")->Clone("mass2d");
	// TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptBinCut.first,ptBinCut.second,centralityBinCut.first,centralityBinCut.second)->Clone("massAllTrigger");
	TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptCut.first,ptCut.second,centralityBinCut.first,centralityBinCut.second)->Clone("massAllTrigger");
	double f = getSBRatio(massAllTrigger);
}

void D0CorrPlotter::plotJetInvM()
{
  TH3F *invMCand = (TH3F *)inputFile->Get("invMCandJets")->Clone("invMCand");
  TH3F *invMBkg = (TH3F *)inputFile->Get("invMBkgJets")->Clone("invMBkg");
  TCanvas *cInvM[9];
  TH1D *invM1D[9][6];
  for(int i=0;i<9;i++)
  {
    cInvM[i] = new TCanvas();
    cInvM[i]->Divide(3,1);
    for(int j=0;j<3;j++)
    {
      cInvM[i]->cd(j+1);
      invM1D[i][j] = (TH1D *)invMCand->ProjectionX(Form("massCand_%i_%i",i,j),j*2,j*2+1,i+1,i+1)->Clone(Form("massCand_%i_%i",i,j));
      invM1D[i][j]->GetXaxis()->SetRangeUser(1.7,2.1);
      invM1D[i][j]->SetNormFactor(1);
      invM1D[i][j]->Draw();
      // cInvM[i]->cd(j+4);
      invM1D[i][j+3] = (TH1D *)invMBkg->ProjectionX(Form("massBkg_%i_%i",i,j),j*2,j*2+1,i+1,i+1)->Clone(Form("massBkg_%i_%i",i,j));
      invM1D[i][j+3]->GetXaxis()->SetRangeUser(1.7,2.1);
      invM1D[i][j+3]->SetNormFactor(1);
      invM1D[i][j+3]->SetLineColor(2);
      // invM1D[i][j+3]->Draw("same");
    }
  }
}
// void D0CorrPlotter::plotCorrelation()
// {
//
// 	TAxis* aDeltaTrigger = candCorrelation->GetXaxis();
// 	candCorrelation->Rebin(1);
//   double pi = 3.1415926;
// 	int bin_d1 = candCorrelation->FindBin(-0.5*pi);
// 	int bin_d2 = candCorrelation->FindBin(0);
// 	int bin_d3 = candCorrelation->FindBin(0.5*pi);
// 	int bin_d4 = candCorrelation->FindBin(1.0*pi);
// 	int bin_d5 = candCorrelation->FindBin(1.5*pi);
// //	aDeltaTrigger->SetBit(TAxis::kLabelsHori);
// 	aDeltaTrigger->SetBinLabel(bin_d1,"-#pi/2");
// 	aDeltaTrigger->SetBinLabel(bin_d2,"0");
// 	aDeltaTrigger->SetBinLabel(bin_d3,"#pi/2");
// 	aDeltaTrigger->SetBinLabel(bin_d4,"#pi");
// 	aDeltaTrigger->SetBinLabel(bin_d5,"3#pi/2");
// 	aDeltaTrigger->SetBit(TAxis::kLabelsHori);
// 	aDeltaTrigger->SetLabelSize(0.075);
//
//   candCorrelation->Draw();
//
// 	TCanvas *bkgCanvas = new TCanvas();
//   bkgCorrelation->Draw();
//
// 	TCanvas *candBkgCanvas = new TCanvas();
// 	candCorrelation->Draw();
// 	candCorrelation->GetYaxis()->SetTitle("(1/N_{trig})(dN_{trig}/d#Delta#phi)");
// 	candCorrelation->GetXaxis()->SetTitle("#Delta#phi");
// 	bkgCorrelation->Draw("same""P");
// 	TLegend *leg = new TLegend(0.2,0.6,0.8,0.8);
// 	leg->AddEntry(candCorrelation,"candidate correlation");
// 	leg->AddEntry(bkgCorrelation,"bkg correlation");
// 	leg->Draw("same");
//
// 	TCanvas *signalCorrelationCanvas = new TCanvas();
// 	signalCorrelation->Draw();
// 	bkgCorrelation->Draw("same");
// }

void D0CorrPlotter::plotSignificance(TH2D *mass2D)
{
	TH1D *massPt[9];
	double x[5] = {1,2,3,4,5};
	double y[5] = {0};
	double z[5] = {0};
	TCanvas *fitCanvas = new TCanvas();
	fitCanvas->Divide(3,2);
  double fitResult[3];
	for(int i=0;i<5;i++)
	{
		massPt[i] = (TH1D *)mass2D->ProjectionX(Form("massPt_%i",i),i+2,10)->Clone(Form("massPt_%i",i));
		fit_hist(massPt[i],fitCanvas,i+1,3,fitResult);
		y[i] = fitResult[0]/fitResult[1];
		z[i] = fitResult[0]/sqrt(fitResult[1]);
	}
	TGraph *sbplot = new TGraph(5,x,y);		
	TGraph *sigplot = new TGraph(5,x,z);		
	TCanvas *significanceCanvas = new TCanvas();
	significanceCanvas->Divide(2,1);
	significanceCanvas->cd(1);
	sbplot->Draw();
	significanceCanvas->cd(2);
	sigplot->Draw();
}
double D0CorrPlotter::getSBRatio(TH1D *massHisto)
{
  double fitResult[3];
  TCanvas *fitCanvas = new TCanvas();
	fit_hist(massHisto,fitCanvas,1,3,fitResult);
  mSOverC = fitResult[0]/fitResult[1];
  fitCanvas->Close();
  return mSOverC;
}




//pair<double,double> fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3])
void  D0CorrPlotter::fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3])
{

	ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000); 
	cfg->cd(iptbin);

	cout << "/////////////////////////////////////********        i             **************         " << iptbin << endl;
	histo->SetMarkerSize(0.8);
	histo->SetLineColor(2);
	histo->SetMarkerColor(2);
	histo->SetMarkerStyle(20);
	histo->GetXaxis()->SetNdivisions(505);
	histo->GetXaxis()->SetTitle("m_{#piK} (GeV)");
	histo->GetYaxis()->SetTitle("Counts");

	double fit_range_low = 1.7;//eff_fit_range_low[iptbin];
	double fit_range_high = 2.0;//eff_fit_range_high[iptbin];

	histo->GetXaxis()->SetRangeUser(fit_range_low, fit_range_high);

	//.. fit with a Gaussian and pol
	TF1* fit_fun = new TF1("fit_fun", "gausn(0) + pol2(3)", fit_range_low, fit_range_high);
	//TF1* fit_fun = new TF1("fit_fun", "gausn(0) + expo(3)", fit_range_low, fit_range_high);
	float max = histo->GetMaximum();
	histo->SetMaximum(1.1 * max);

	float p0 = 1000, p1 = 1.87, p2 = 0.02;
	float p0_L = 0, p1_L = 1.84, p2_L = 0;
	float p0_H = 2*max, p1_H = 1.9, p2_H = 0.05;

	float p3 = -1. * max, p4 = max, p5 = -1. * max;

	int pass = 0;
	int fittingtry = 0;

	char sig_print[100], chi2_print[100], mean_print[100], sigma_print[100],sb_ratio[100],sign_print[100];

	while (!pass) {

		fit_fun->SetParameter(0, p0);
		fit_fun->SetParameter(1, p1);
		fit_fun->SetParameter(2, p2);

		//.. fit constraint ..
		fit_fun->SetParLimits(0, p0_L, p0_H);
		fit_fun->SetParLimits(1, p1_L, p1_H);
		fit_fun->SetParLimits(2, p2_L, p2_H);

		//        fit_fun->SetParameter(3, p3);
		//		fit_fun->SetParameter(4, p4);
		//		fit_fun->SetParameter(5, p5);

		if( fittingtry == 0 )
			histo->Fit(fit_fun,"L","", fit_range_low, fit_range_high);
		else 
			histo->Fit(fit_fun,"L","", fit_range_low, fit_range_high);

		//.. draw foreground and background ..
		histo->Draw();

		TF1* fit_fun_1st = (TF1*)fit_fun->Clone("fit_fun_1st");
		fit_fun_1st->SetParameter(3, 0);
		fit_fun_1st->SetParameter(4, 0);
		fit_fun_1st->SetParameter(5, 0);
		//        fit_fun_1st->Draw("same");


		TF1* fit_fun_bg = (TF1*)fit_fun->Clone("fit_fun_bg");
		//        TF1* fit_fun_bg = new TF1("fit_fun_bg", fitfunction, cut_m_low, cut_m_high, 6);
		fit_fun_bg->SetParameter(0, 0);
		fit_fun_bg->SetParameter(1, 0);
		fit_fun_bg->SetParameter(2, 0);
		//		fit_fun_bg->SetParameter(3, fit_fun->GetParameter(3));
		//		fit_fun_bg->SetParameter(4, fit_fun->GetParameter(4));
		//		fit_fun_bg->SetParameter(5, fit_fun->GetParameter(5));


		fit_fun_bg->SetLineColor(8);
		fit_fun_bg->SetLineStyle(2);
		fit_fun_bg->Draw("same");


		fittingtry++;

		//    if( ptbins[iptbin] > lowrange && ptbins[iptbin+1] < highrange )
		{
			float binwidth = 0.01;//histo->GetBinWidth(10);
			//float ptbinwidth = ptbins[iptbin+1] - ptbins[iptbin];
			//counts->SetBinContent( iptbin+1, fit_fun->GetParameter(0)/( binwidth * ptbinwidth ));
			//counts->SetBinError( iptbin+1, fit_fun->GetParError(0)/( binwidth * ptbinwidth ));

			float Nsig = fit_fun->GetParameter(0)/( binwidth );
			float err_Nsig = fit_fun->GetParError(0)/( binwidth );
			float fitchi2 = fit_fun->GetChisquare();
			float fitmeanerror = fit_fun->GetParError(1);
			float fitsigmaerror = fit_fun->GetParError(2);
			int noffreepara = fit_fun->GetNumberFreeParameters();
			int noffitpoints = fit_fun->GetNumberFitPoints();

			float fitmean = fit_fun->GetParameter(1);
			float fitsigma = fit_fun->GetParameter(2);

			//hfg_masssigma->SetBinContent(iptbin+1, fitsigma);
			//hfg_masssigma->SetBinError(iptbin+1, fitsigmaerror);

			cout << " fitchi2: " << fitchi2 << "   noffreepara: " << noffreepara << "  noffitpoints: " << noffitpoints << endl;

			//      if( !isMC )
			//	sprintf( sig_print,"N_{sig}: %7.1f#pm%7.1f", Nsig, err_Nsig);
			//      else
			sprintf( sig_print,"N_{sig}: %7.2f#pm%7.2f", Nsig, err_Nsig);
			sprintf( chi2_print, "#chi^{2}#/d.o.f: %3.2f", fitchi2/( noffitpoints - noffreepara));
			sprintf( mean_print, "mean: %6.4f#pm%6.4f", fitmean,fitmeanerror);
			sprintf( sigma_print, "#sigma: %6.4f#pm%6.4f", fitsigma,fitsigmaerror);
			cout<<fitmean<<"\t"<<fitsigma<<endl;

			TF1 *g = new TF1("g","gausn(0)",1.6,2.1);
			g->SetParameters(Nsig*binwidth,fitmean,fitsigma);
			float inteSig = g->Integral(fitmean-nSigma*fitsigma,fitmean+nSigma*fitsigma) / (binwidth);
			float inteSig_err = err_Nsig * inteSig/Nsig;
			pair<double,double> fitResult (inteSig,inteSig_err);
			fitArray[0] = inteSig;
			fitArray[1] = fit_fun->Integral(fitmean-nSigma*fitsigma,fitmean+nSigma*fitsigma)/binwidth;
			fitArray[2] = fit_fun->IntegralError(fitmean-nSigma*fitsigma,fitmean+nSigma*fitsigma)/binwidth;
			//fout<<iptbin<<"fitmean = "<<fitmean<<endl;
			//fout<<iptbin<<"fitsigma = "<<fitsigma<<endl;
			double sbratio = fitArray[0]/fitArray[1];
      double sign = fitArray[0]/sqrt(fitArray[1]);
			sprintf( sb_ratio, "s/(s+b): %6.4f", sbratio);
			sprintf( sign_print, "significance: %6.4f", sign);


			if (fittingtry == 2)
			{
				TLatex Tl;
				Tl.SetNDC();
				Tl.SetTextAlign(12);
				Tl.SetTextSize(0.06);
				Tl.DrawLatex(0.15,0.8, sig_print);
				Tl.DrawLatex(0.15,0.7, chi2_print);
				Tl.DrawLatex(0.15,0.6, sb_ratio);
				Tl.DrawLatex(0.55,0.8, mean_print);
				Tl.DrawLatex(0.55,0.7, sigma_print);
				Tl.DrawLatex(0.55,0.6, sign_print);
			}

		}

		if (fittingtry == 2)  
		{
			pass = 1;

		}
		if(!pass) {
			p0 = fit_fun->GetParameter(0);
			p1 = fit_fun->GetParameter(1);
			p2 = fit_fun->GetParameter(2);
			//            p1_L = 1.84, p2_L = 0;
			//            p1_H = 1.9, p2_H = 0.05;
		}
	}
	double fitpar[6];
	for(int i=0;i<6;i++)
	{
		//fout<<fit_fun->GetParameter(i)<<",";
		fitpar[i] = fit_fun->GetParameter(i);
	}
	//fout<<fitpar[0]<<"*exp(-0.5*((x-"<<fitpar[1]<<")/"<<fitpar[2]<<")**2)/(sqrt(2*pi)*"<<fitpar[2]<<"))/("<<fitpar[3]<<"+x*"<<fitpar[4]<<"+x*x*"<<fitpar[5]<<")"<<endl;
	//fout<<endl;


	//return fitResult;
}

void D0CorrPlotter::setCorrAxis(TH1D *corrPlot)
{
  double pi = TMath::Pi();
  int bin1 = corrPlot->FindBin(-0.5*pi);
  int bin2 = corrPlot->FindBin(0);
  int bin3 = corrPlot->FindBin(0.5*pi);
  int bin4 = corrPlot->FindBin(1.0*pi);
  int bin5 = corrPlot->FindBin(1.5*pi);
  TAxis* xAxis = corrPlot->GetXaxis();
  TAxis* yAxis = corrPlot->GetYaxis();  
  yAxis->SetNdivisions(505);
  xAxis->SetNdivisions(505);
  xAxis->SetTitle("#Delta #phi");
  yAxis->SetTitle("(1/N_{D^{0}})(dN_{jet}/d#Delta#phi)");

  xAxis->SetTitleOffset(1.1);
  // xAxis->SetTitleSize(0.6);
  xAxis->SetTitleFont(42);
  xAxis->SetLabelOffset(0.01);
  // xAxis->SetLabelSize(0.6);
  xAxis->SetLabelFont(42);
  yAxis->SetTitleOffset(1.1);
  // yAxis->SetTitleSize(0.6);
  yAxis->SetTitleFont(42);
  yAxis->SetLabelOffset(0.01);
  // yAxis->SetLabelSize(0.6);
  yAxis->SetLabelFont(42);

  xAxis->SetBit(TAxis::kLabelsHori);
  xAxis->SetBinLabel(bin1,"-#pi/2");
  xAxis->SetBinLabel(bin2,"0");
  xAxis->SetBinLabel(bin3,"#pi/2");
  xAxis->SetBinLabel(bin4,"#pi");
  xAxis->SetBinLabel(bin5,"3#pi/2");
  xAxis->SetBit(TAxis::kLabelsHori);
  xAxis->SetLabelSize(0.075);
}
