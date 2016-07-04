#include "d0_corr_plotter.cc"
#include "TStyle.h"
#include "TLegend.h"
void plotD0Jet()
{
  int inputRebinFactor = 50;
  D0CorrPlotter *plotter = new D0CorrPlotter();
  // plotter->init("root-files/dJCorrPt16.root");
  plotter->init("root-files/d0-corr-trig-3GeV.root");

  TCanvas *cJet[6];
  TLegend *legJet[6][9];
  TLegend *legCorr[6];
  TH1D *candJet[6][9];
  TH1D *bkgJet[6][9];
  TH1D *diff[6][9];
  TH1D *corrC[6];
  TH1D *corrB[6];
  TH1D *corrS[6];

  TCanvas *cMB = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cMB->Divide(3,2);
  int centBin[10] = {80,70,60,50,40,30,20,10,5,0};
  for(int j=0;j<6;j++)
  {
    cJet[j]  = new TCanvas();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetHistMinimumZero(kTRUE);
    cJet[j]->Divide(3,2);
    for(int i=0;i<6;i++)
    {
      pair<int,int> inputCentralityBin(i+4,i+4);
      pair<double,double> inputPtCut(j+3,100);
      plotter->getJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor);
      candJet[j][i] = plotter->getCandJetCorrelation();
      bkgJet[j][i] = plotter->getBkgJetCorrelation();
      cJet[j]->cd(i+1);

      candJet[j][i]->SetLineColor(1);
      candJet[j][i]->SetLineStyle(2);
      bkgJet[j][i]->SetLineColor(1);
      candJet[j][i]->Draw();
      bkgJet[j][i]->Draw("same");
      double sOverC = plotter->getSBRatio();
      cout<<"sOverC = "<<sOverC<<endl;
      TH1D *bkgTemp = (TH1D *)bkgJet[j][i]->Clone("bkgTemp");
      bkgTemp->Scale(1-sOverC);
      diff[j][i] = (TH1D *)candJet[j][i]->Clone(Form("diff%i",i));
      diff[j][i]->Add(bkgTemp,-1);
      diff[j][i]->Scale(1/sOverC);
      diff[j][i]->SetLineColor(2);
      diff[j][i]->SetLineStyle(1);
      diff[j][i]->Draw("same");
      legJet[j][i] = new TLegend(0.2,0.2,0.8,0.4);
      legJet[j][i]->SetHeader(Form("jet p_{T}>%iGeV/c,%i<centrality<%i",j+3,centBin[i+4],centBin[i+3]));
      legJet[j][i]->AddEntry(candJet[j][i],"D^{0} candidate");
      legJet[j][i]->AddEntry(bkgJet[j][i],"side-band bkg");
      legJet[j][i]->AddEntry(diff[j][i],"D^{0} signal");
      legJet[j][i]->Draw("same");
    }
    corrS[j] = (TH1D *)diff[j][0]->Clone("corrS");
    corrB[j] = (TH1D *)bkgJet[j][0]->Clone("corrB");
    corrC[j] = (TH1D *)candJet[j][0]->Clone("corrC");
    for(int i=1;i<6;i++)
    {
      corrS[j]->Add(diff[j][i]);
      corrB[j]->Add(bkgJet[j][i]);
      corrC[j]->Add(candJet[j][i]);
    }
    cMB->cd(j+1);
    corrS[j]->Draw("pe");
    corrB[j]->Draw("pe,same");
    // corrC[j]->Draw("pe,same");
    corrS[j]->SetLineColor(2);
    corrB[j]->SetLineColor(1);
    corrC[j]->SetLineColor(1);
    corrC[j]->SetLineStyle(2);
    legCorr[j] = new TLegend(0.3,0.2,0.7,0.35);
    legCorr[j]->SetHeader(Form("jet p_{T}>%iGeV/c,MB events",j+3));
    legCorr[j]->AddEntry(corrS[j],"D^{0} signal");
    legCorr[j]->AddEntry(corrB[j],"Side-band");
    legCorr[j]->Draw("same");
  }
}
void plotD0Hadron()
{
  int inputRebinFactor = 100;
  D0CorrPlotter *plotter = new D0CorrPlotter();
  // plotter->init("root-files/dJCorrPt16.root");
  plotter->init("root-files/d0-corr-trig-3GeV.root");

  TCanvas *cHadron[6];
  TLegend *legHadron[6][9];
  TLegend *legCorr[6];
  TH1D *candHadron[6][9];
  TH1D *bkgHadron[6][9];
  TH1D *diff[6][9];
  TH1D *corrC[6];
  TH1D *corrB[6];
  TH1D *corrS[6];

  TCanvas *cMB = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  cMB->Divide(3,2);
  int centBin[10] = {80,70,60,50,40,30,20,10,5,0};
  for(int j=0;j<6;j++)
  {
    cHadron[j]  = new TCanvas();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetHistMinimumZero(kTRUE);
    cHadron[j]->Divide(3,2);
    for(int i=0;i<6;i++)
    {
      pair<int,int> inputCentralityBin(i+4,i+4);
      pair<double,double> inputPtCut(0.5*j+1,100);
      plotter->getHadronCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor);
      candHadron[j][i] = plotter->getCandHadronCorrelation();
      bkgHadron[j][i] = plotter->getBkgHadronCorrelation();
      cHadron[j]->cd(i+1);

      candHadron[j][i]->SetLineColor(1);
      candHadron[j][i]->SetLineStyle(2);
      bkgHadron[j][i]->SetLineColor(1);
      candHadron[j][i]->Draw();
      bkgHadron[j][i]->Draw("same");
      double sOverC = plotter->getSBRatio();
      cout<<"sOverC = "<<sOverC<<endl;
      TH1D *bkgTemp = (TH1D *)bkgHadron[j][i]->Clone("bkgTemp");
      bkgTemp->Scale(1-sOverC);
      diff[j][i] = (TH1D *)candHadron[j][i]->Clone(Form("diff%i",i));
      diff[j][i]->Add(bkgTemp,-1);
      diff[j][i]->Scale(1/sOverC);
      diff[j][i]->SetLineColor(2);
      diff[j][i]->SetLineStyle(1);
      diff[j][i]->Draw("same");
      legHadron[j][i] = new TLegend(0.2,0.2,0.8,0.4);
      legHadron[j][i]->SetHeader(Form("hadron p_{T}>%fGeV/c,%i<centrality<%i",0.5*j+1,centBin[i+4],centBin[i+3]));
      legHadron[j][i]->AddEntry(candHadron[j][i],"D^{0} candidate");
      legHadron[j][i]->AddEntry(bkgHadron[j][i],"side-band bkg");
      legHadron[j][i]->AddEntry(diff[j][i],"D^{0} signal");
      legHadron[j][i]->Draw("same");
    }
    corrS[j] = (TH1D *)diff[j][0]->Clone("corrS");
    corrB[j] = (TH1D *)bkgHadron[j][0]->Clone("corrB");
    corrC[j] = (TH1D *)candHadron[j][0]->Clone("corrC");
    for(int i=1;i<6;i++)
    {
      corrS[j]->Add(diff[j][i]);
      corrB[j]->Add(bkgHadron[j][i]);
      corrC[j]->Add(candHadron[j][i]);
    }
    cMB->cd(j+1);
    corrS[j]->Draw("pe");
    corrB[j]->Draw("pe,same");
    // corrC[j]->Draw("pe,same");
    corrS[j]->SetLineColor(2);
    corrB[j]->SetLineColor(1);
    corrC[j]->SetLineColor(1);
    corrC[j]->SetLineStyle(2);
    legCorr[j] = new TLegend(0.3,0.2,0.7,0.35);
    legCorr[j]->SetHeader(Form("hadron p_{T}>%iGeV/c,MB events",j+3));
    legCorr[j]->AddEntry(corrS[j],"D^{0} signal");
    legCorr[j]->AddEntry(corrB[j],"Side-band");
    legCorr[j]->Draw("same");
  }
}

