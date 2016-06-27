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

#include "d0jet_corr_plotter.h"
using namespace std;

void D0jetCorrPlotter::init(TString inputFileName)
{
	inputFile = new TFile(inputFileName.Data());
  mLog.open("d0jet_corr_plotter.log");
}

void D0jetCorrPlotter::finish()
{
}

void D0jetCorrPlotter::getPxCut(double pxCutPercentage)
{
	TH2F *dCorPxPlus = (TH2F *)inputFile->Get("dCorPxPlus")->Clone("dCorPxPlus");
	TH2F *dCorPxMinus = (TH2F *)inputFile->Get("dCorPxMinus")->Clone("dCorPxMinus");
	TH1F *corPxPlus[9];
	TH1F *corPxMinus[9];
	TCanvas *mPxCanvas = new TCanvas();
	mPxCanvas->Divide(3,3);
  double pxPlusCut[9],pxMinusCut[9];
	for(int i=0;i<9;++i)
	{
		corPxPlus[i] = (TH1F *)(dCorPxPlus->ProjectionX(Form("pxplus_%i",i),i+1,i+1)->Clone(Form("pxPlus_%i",i)));
		mPxCanvas->cd(i+1);
		corPxPlus[i]->Draw();
		double sumPxPlus = 0;
		double intePxPlus = corPxPlus[i]->Integral();
		int jPxBin;
		for(jPxBin=1 ; sumPxPlus<pxCutPercentage*intePxPlus ; jPxBin++)
			sumPxPlus += corPxPlus[i]->GetBinContent(jPxBin);
		pxPlusCut[i] = corPxPlus[i]->GetBinCenter(jPxBin);
		corPxMinus[i] = (TH1F *)dCorPxMinus->ProjectionX(Form("pxminus_%i",i),i+1,i+1)->Clone(Form("pxMinus_%i",i));
		corPxMinus[i]->Draw("same");
		corPxMinus[i]->SetLineColor(2);
		double sumPxMinus = 0;
		double intePxMinus = corPxMinus[i]->Integral();
		for(jPxBin=1 ; sumPxMinus<pxCutPercentage*intePxMinus ; jPxBin++)
			sumPxMinus += corPxMinus[i]->GetBinContent(jPxBin);
		pxMinusCut[i] = corPxMinus[i]->GetBinCenter(jPxBin);
	}
  cout<<"minus px cut:";
	for(int i=0;i<9;++i)
		cout<<pxMinusCut[i]<<",";
  cout<<endl;
  cout<<"plus px cut:";
	for(int i=0;i<9;++i)
		cout<<pxPlusCut[i]<<",";
  cout<<endl;
}


void D0jetCorrPlotter::getCorrelation(int rebinFactor)
{
  pair<int,int> ptBinCut(2,10);
	TH2F *corBkg2D = (TH2F *)inputFile->Get("corBkg")->Clone("corBkg2D");
	TH2F *corCand2D = (TH2F *)inputFile->Get("corCand")->Clone("corCand2D");
  TH2F *candCount = (TH2F *)inputFile->Get("candCount")->Clone("candCount");
  TH2F *bkgCount = (TH2F *)inputFile->Get("bkgCount")->Clone("bkgCount");
  double nD0Cand = candCount->ProjectionX("candCount1D")->Integral();
  double nD0Bkg= bkgCount->ProjectionX("bkgCount1D")->Integral();
  cout<<"# of candidate D0 = "<<nD0Cand<<endl;
  cout<<"# of background D0 = "<<nD0Bkg<<endl;
	candCorrelation = (TH1D *)corCand2D->ProjectionX("candCorrelation",2,2)->Clone("candCorrelation");
	bkgCorrelation = (TH1D *)corBkg2D->ProjectionX("bkgCorrelation",2,2)->Clone("bkgCorrelation");
	candCorrelation->Rebin(rebinFactor);
	candCorrelation->Scale(1./candCorrelation->GetBinWidth(1)/nD0Cand);
	bkgCorrelation->Rebin(rebinFactor);
	bkgCorrelation->Scale(1./bkgCorrelation->GetBinWidth(1)/nD0Bkg);

	TH3F *mass2d = (TH3F *)inputFile->Get("massPt")->Clone("mass2d");
	TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptBinCut.first,ptBinCut.second,1,10)->Clone("massAllTrigger");
	double f = getSBRatio(massAllTrigger);
}
void D0jetCorrPlotter::getCorrelation(pair<int,int> &centralityBinCut,int rebinFactor)
{
  pair<int,int> ptBinCut(2,10);
	TH2F *corBkg2D = (TH2F *)inputFile->Get("corBkg")->Clone("corBkg2D");
	TH2F *corCand2D = (TH2F *)inputFile->Get("corCand")->Clone("corCand2D");
  TH2F *candCount = (TH2F *)inputFile->Get("candCount")->Clone("candCount");
  TH2F *bkgCount = (TH2F *)inputFile->Get("bkgCount")->Clone("bkgCount");

  double nD0Cand = candCount->ProjectionX("candCount1D",centralityBinCut.first,centralityBinCut.second)->Integral();
  double nD0Bkg= bkgCount->ProjectionX("bkgCount1D",centralityBinCut.first,centralityBinCut.second)->Integral();
  cout<<"# of candidate D0 = "<<nD0Cand<<endl;
  cout<<"# of background D0 = "<<nD0Bkg<<endl;

	candCorrelation = (TH1D *)corCand2D->ProjectionX("candCorrelation",centralityBinCut.first,centralityBinCut.second)->Clone("candCorrelation");
	bkgCorrelation = (TH1D *)corBkg2D->ProjectionX("bkgCorrelation",centralityBinCut.first,centralityBinCut.second)->Clone("bkgCorrelation");
  // get candidate correlation
	candCorrelation->Rebin(rebinFactor);
	candCorrelation->Scale(1./candCorrelation->GetBinWidth(1)/nD0Cand);
	bkgCorrelation->Rebin(rebinFactor);
	bkgCorrelation->Scale(1./bkgCorrelation->GetBinWidth(1)/nD0Bkg);

	TH3F *mass2d = (TH3F *)inputFile->Get("massPt")->Clone("mass2d");
	// TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptBinCut.first,ptBinCut.second,centralityBinCut.first,centralityBinCut.second)->Clone("massAllTrigger");
	TH1D *massAllTrigger = (TH1D *)mass2d->ProjectionX("massAllTrigger",ptBinCut.first,ptBinCut.second,centralityBinCut.first,centralityBinCut.second)->Clone("massAllTrigger");
	double f = getSBRatio(massAllTrigger);

}
void D0jetCorrPlotter::plotJetInvM()
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
void D0jetCorrPlotter::plotCorrelation()
{
	candCorrelationClose->SetLineColor(2);
	double pi = TMath::Pi();
	int bin1 = candCorrelationClose->FindBin(-0.5*pi);
	int bin2 = candCorrelationClose->FindBin(0);
	int bin3 = candCorrelationClose->FindBin(0.5*pi);
	int bin4 = candCorrelationClose->FindBin(1.0*pi);
	int bin5 = candCorrelationClose->FindBin(1.5*pi);

	TAxis* aCloseTrigger = candCorrelationClose->GetXaxis();
	aCloseTrigger->SetBit(TAxis::kLabelsHori);
	aCloseTrigger->SetBinLabel(bin1,"-#pi/2");
	aCloseTrigger->SetBinLabel(bin2,"0");
	aCloseTrigger->SetBinLabel(bin3,"#pi/2");
	aCloseTrigger->SetBinLabel(bin4,"#pi");
	aCloseTrigger->SetBinLabel(bin5,"3#pi/2");
	aCloseTrigger->SetBit(TAxis::kLabelsHori);
	aCloseTrigger->SetLabelSize(0.075);

	TAxis* aFarTrigger = candCorrelationFar->GetXaxis();
	aFarTrigger->SetBit(TAxis::kLabelsHori);
	aFarTrigger->SetBinLabel(bin1,"-#pi/2");
	aFarTrigger->SetBinLabel(bin2,"0");
	aFarTrigger->SetBinLabel(bin3,"#pi/2");
	aFarTrigger->SetBinLabel(bin4,"#pi");
	aFarTrigger->SetBinLabel(bin5,"3#pi/2");
	aFarTrigger->SetBit(TAxis::kLabelsHori);
	aFarTrigger->SetLabelSize(0.075);

	TAxis* aCloseBkgTrigger = bkgCorrelationClose->GetXaxis();
	aCloseBkgTrigger->SetBit(TAxis::kLabelsHori);
	aCloseBkgTrigger->SetBinLabel(bin1,"-#pi/2");
	aCloseBkgTrigger->SetBinLabel(bin2,"0");
	aCloseBkgTrigger->SetBinLabel(bin3,"#pi/2");
	aCloseBkgTrigger->SetBinLabel(bin4,"#pi");
	aCloseBkgTrigger->SetBinLabel(bin5,"3#pi/2");
	aCloseBkgTrigger->SetBit(TAxis::kLabelsHori);
	aCloseBkgTrigger->SetLabelSize(0.075);

	TAxis* aDeltaTrigger = candCorrelation->GetXaxis();
	candCorrelation->Rebin(1);
	int bin_d1 = candCorrelation->FindBin(-0.5*pi);
	int bin_d2 = candCorrelation->FindBin(0);
	int bin_d3 = candCorrelation->FindBin(0.5*pi);
	int bin_d4 = candCorrelation->FindBin(1.0*pi);
	int bin_d5 = candCorrelation->FindBin(1.5*pi);
//	aDeltaTrigger->SetBit(TAxis::kLabelsHori);
	aDeltaTrigger->SetBinLabel(bin_d1,"-#pi/2");
	aDeltaTrigger->SetBinLabel(bin_d2,"0");
	aDeltaTrigger->SetBinLabel(bin_d3,"#pi/2");
	aDeltaTrigger->SetBinLabel(bin_d4,"#pi");
	aDeltaTrigger->SetBinLabel(bin_d5,"3#pi/2");
	aDeltaTrigger->SetBit(TAxis::kLabelsHori);
	aDeltaTrigger->SetLabelSize(0.075);

	TCanvas *candCanvas = new TCanvas();
  candCanvas->Divide(2,1);
  candCanvas->cd(1);
	candCorrelationClose->GetYaxis()->SetTitle("(1/N_{trig})(dN_{trig}/d#Delta#phi)");
	candCorrelationClose->GetXaxis()->SetTitle("#Delta#phi");
	candCorrelationClose->Draw();
	candCorrelationFar->Draw("same");
	candCorrelationClose->SetLineColor(2);
	TLegend *leg2 = new TLegend(0.2,0.6,0.8,0.8);
	leg2->AddEntry(candCorrelationClose,"Close Region");
	leg2->AddEntry(candCorrelationFar,"Far Region");
	leg2->Draw("same");
  candCanvas->cd(2);
  candCorrelation->Draw();

	TCanvas *bkgCanvas = new TCanvas();
  bkgCanvas->Divide(2,1);
  bkgCanvas->cd(1);
	bkgCorrelationClose->GetYaxis()->SetTitle("(1/N_{trig})(dN_{trig}/d#Delta#phi)");
	bkgCorrelationClose->GetXaxis()->SetTitle("#Delta#phi");
	bkgCorrelationClose->Draw();
	bkgCorrelationFar->Draw("same");
	bkgCorrelationClose->SetLineColor(2);
	leg2->Draw("same");
  bkgCanvas->cd(2);
  bkgCorrelation->Draw();

	TCanvas *candBkgCanvas = new TCanvas();
	candCorrelation->Draw();
	candCorrelation->GetYaxis()->SetTitle("(1/N_{trig})(dN_{trig}/d#Delta#phi)");
	candCorrelation->GetXaxis()->SetTitle("#Delta#phi");
	bkgCorrelation->Draw("same""P");
	TLegend *leg = new TLegend(0.2,0.6,0.8,0.8);
	leg->AddEntry(candCorrelation,"candidate correlation");
	leg->AddEntry(bkgCorrelation,"bkg correlation");
	leg->Draw("same");

	TCanvas *signalCorrelationCanvas = new TCanvas();
	signalCorrelation->Draw();
	bkgCorrelation->Draw("same");
}

void D0jetCorrPlotter::plotSignificance(TH2D *mass2D)
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
double D0jetCorrPlotter::getSBRatio(TH1D *massHisto)
{
  double fitResult[3];
  TCanvas *fitCanvas = new TCanvas();
	fit_hist(massHisto,fitCanvas,1,3,fitResult);
  mSOverC = fitResult[0]/fitResult[1];
  return mSOverC;
}




//pair<double,double> fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3])
void  D0jetCorrPlotter::fit_hist(TH1D *histo, TCanvas *cfg, int iptbin ,double nSigma,double fitArray[3])
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
