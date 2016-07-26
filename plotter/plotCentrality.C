#include "d0_corr_plotter.cc"
#include "TStyle.h"
#include "TLegend.h"

TH1D *plotHadronJetInclusive(int );
TH1D *plotD0JetInclusiveME(double ,double );
void plotD0Jet()
{
  int inputRebinFactor = 50;
  D0CorrPlotter *plotter = new D0CorrPlotter();
  // plotter->init("root-files/dJCorrPt16.root");
  // plotter->init("root-files/d0-corr-trig-3GeV.root");
   plotter->init("root-files/test.root");

  TCanvas *cJet[6];
  TLegend *legJet[6][9];
  TLegend *legCorr[6];
  TH1D *candJet[6][9];
  TH1D *hadronJet[6][9];
  TH1D *bkgJet[6][9];
  TH1D *diff[6][9];
  TH1D *corrC[6];
  TH1D *corrB[6];
  TH1D *corrS[6];
  TH1D *corrH[6];

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
    // gStyle->SetHistMinimumZero(kTRUE);
    cJet[j]->Divide(3,3);
    for(int i=0;i<9;i++)
    {
      pair<int,int> inputCentralityBin(i+1,i+1);
      pair<double,double> inputPtCut(j+5,100);
      pair<double,double> inputD0PtCut(1,3);
      plotter->getJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor,inputD0PtCut);
      plotter->getHadronJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor,inputD0PtCut);
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
      legJet[j][i]->SetHeader(Form("jet p_{T}>%iGeV/c,%i<centrality<%i",j+3,centBin[i+1],centBin[i+0]));
      legJet[j][i]->AddEntry(candJet[j][i],"D^{0} candidate");
      legJet[j][i]->AddEntry(bkgJet[j][i],"side-band bkg");
      legJet[j][i]->AddEntry(diff[j][i],"D^{0} signal");
      legJet[j][i]->Draw("same");
    }
    corrS[j] = (TH1D *)diff[j][0]->Clone("corrS");
    corrB[j] = (TH1D *)bkgJet[j][0]->Clone("corrB");
    corrC[j] = (TH1D *)candJet[j][0]->Clone("corrC");
    corrH[j] = (TH1D *)hadronJet[j][0]->Clone("corrC");
    for(int i=1;i<9;i++)
    {
      corrS[j]->Add(diff[j][i]);
      corrB[j]->Add(bkgJet[j][i]);
      corrC[j]->Add(candJet[j][i]);
      corrH[j]->Add(hadronJet[j][i]);
    }
    cMB->cd(1);
    // corrS[j]->Scale(0.1);
    // corrC[j]->Scale(0.1);
    // corrB[j]->Scale(0.1);
    // corrH[j]->Scale(0.1);
    corrS[j]->Draw("pe");
    corrB[j]->Draw("pe,same");
    // corrH[j]->Draw("pe,same");
    corrC[j]->Draw("pe,same");
    corrS[j]->SetLineColor(2);
    corrB[j]->SetLineColor(1);
    corrC[j]->SetLineColor(1);
    corrC[j]->SetLineStyle(2);

    corrH[j]->SetLineColor(4);
    corrH[j]->SetMarkerColor(4);
    corrH[j]->SetMarkerStyle(20);
    legCorr[j] = new TLegend(0.3,0.2,0.7,0.35);
    legCorr[j]->SetHeader(Form("jet p_{T}>%iGeV/c,MB events",j+5));
    legCorr[j]->AddEntry(corrS[j],"D^{0} signal");
    legCorr[j]->AddEntry(corrB[j],"Side-band");
    legCorr[j]->AddEntry(corrC[j],"Candidates");
    legCorr[j]->AddEntry(corrH[j],"Hadron-jet correlation");
    legCorr[j]->Draw("same");
  }
}
void plotD0Hadron()
{
  int inputRebinFactor = 50;
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
  int centBin[10] = {80,70,60,50,40,30,20,10,5,0};
  for(int j=0;j<6;j++)
  {
    cHadron[j]  = new TCanvas();
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    // gStyle->SetHistMinimumZero(kTRUE);
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
    corrC[j]->Draw("pe,same");
    corrS[j]->SetLineColor(2);
    corrB[j]->SetLineColor(1);
    corrC[j]->SetLineColor(1);
    corrC[j]->SetLineStyle(2);
    legCorr[j] = new TLegend(0.3,0.2,0.7,0.35);
    legCorr[j]->SetHeader(Form("hadron p_{T}>%iGeV/c,MB events",j+3));
    legCorr[j]->AddEntry(corrS[j],"D^{0} signal");
    legCorr[j]->AddEntry(corrC[j],"Candidates");
    legCorr[j]->AddEntry(corrB[j],"Side-band");
    legCorr[j]->Draw("same");
  }
}
void plotD0JetInclusive(double dPtLow,double dPtHigh)
{
  int inputRebinFactor = 50;
  D0CorrPlotter *plotter = new D0CorrPlotter();
  plotter->init("root-files/d0CorrelationGlobal.root");

  TCanvas *cJet;
  TLegend *legJet[9];
  TLegend *legCorr;
  TH1D *candJet[9];
  TH1D *hadronJet[9];
  TH1D *bkgJet[9];
  TH1D *diff[9];
  TH1D *corrC;
  TH1D *corrB;
  TH1D *corrS;
  TH1D *corrH;

  TCanvas *cMB = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int centBin[10] = {80,70,60,50,40,30,20,10,5,0};
  cJet = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  // gStyle->SetHistMinimumZero(kTRUE);
  cJet->Divide(3,3);
  for(int i=0;i<9;i++)
  {
    pair<int,int> inputCentralityBin(i+1,i+1);
    pair<double,double> inputPtCut(8,100);
    pair<double,double> inputD0PtCut(dPtLow,dPtHigh);
    pair<double,double> inputHadronPtCut(dPtLow,dPtHigh);
    plotter->getJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor,inputD0PtCut);
    candJet[i] = plotter->getCandJetCorrelation();
    bkgJet[i] = plotter->getBkgJetCorrelation();

    cJet->cd(i+1);

    candJet[i]->SetLineColor(1);
    candJet[i]->SetLineStyle(2);
    bkgJet[i]->SetLineColor(1);
    candJet[i]->Draw();
    bkgJet[i]->Draw("same");
    double sOverC = plotter->getSBRatio();
    cout<<"sOverC = "<<sOverC<<endl;
    TH1D *bkgTemp = (TH1D *)bkgJet[i]->Clone("bkgTemp");
    bkgTemp->Scale(1-sOverC);
    diff[i] = (TH1D *)candJet[i]->Clone(Form("diff%i",i));
    diff[i]->Add(bkgTemp,-1);
    diff[i]->Scale(1/sOverC);
    diff[i]->SetLineColor(2);
    diff[i]->SetLineStyle(1);
    diff[i]->Draw("same");
    legJet[i] = new TLegend(0.2,0.2,0.8,0.4);
    legJet[i]->SetHeader(Form("jet p_{T}>%iGeV/c,%i<centrality<%i",8,centBin[i+1],centBin[i+0]));
    legJet[i]->AddEntry(candJet[i],"D^{0} candidate");
    legJet[i]->AddEntry(bkgJet[i],"side-band bkg");
    legJet[i]->AddEntry(diff[i],"D^{0} signal");
    legJet[i]->Draw("same");
  }
  corrS = (TH1D *)diff[1]->Clone("corrS");
  corrB = (TH1D *)bkgJet[1]->Clone("corrB");
  corrC = (TH1D *)candJet[1]->Clone("corrC");
  for(int i=1;i<9;i++)
  {
    corrS->Add(diff[i]);
    corrB->Add(bkgJet[i]);
    corrC->Add(candJet[i]);
  }
  cMB->cd();
  // corrS->Scale(0.1);
  // corrC->Scale(0.1);
  // corrB->Scale(0.1);
  corrS->Draw("pe");
  corrB->Draw("pe,same");
  corrC->Draw("pe,same");
  corrS->SetLineColor(2);
  corrB->SetLineColor(1);
  corrC->SetLineColor(4);

  legCorr = new TLegend(0.3,0.2,0.7,0.35);
  legCorr->SetHeader(Form("jet p_{T}>%iGeV/c,MB events",8));
  legCorr->AddEntry(corrS,"D^{0} signal");
  legCorr->AddEntry(corrB,"Side-band");
  legCorr->AddEntry(corrC,"Candidates");
  legCorr->Draw("same");

  TCanvas *temp = new TCanvas();
  // temp->Divide(2,1);
  // temp->cd(1);
  corrS->Draw();
  corrS->Add(plotHadronJetInclusive(dPtLow-1),-1);
  // temp->cd(2);
  // plotHadronJetInclusive(dPtLow-1)->Draw();
}

TH1D *plotHadronJetInclusive(int ptBin)
{
  int inputRebinFactor = 50;
  D0CorrPlotter *plotter = new D0CorrPlotter();
  plotter->init("root-files/hadronJet.root");
  TH1D *hHadronJet[6];
  TH1D *hHadronJetCent[6][9];
  double hadronPtCut[6] = {0.6,1.1,2,5,0.65,2.0};
  // TCanvas *cHadronJet = new TCanvas();
  // cHadronJet->Divide(3,2,0.01,0,0);
  // cHadronJet->Divide(3,2);
  // for(int i=0;i<6;i++)
  for(int i=ptBin;i<ptBin+1;i++)
  {
    pair<int,int> inputCentralityBin(1,1);
    pair<double,double> inputPtCut(8,100);
    pair<double,double> inputD0PtCut(hadronPtCut[i]-0.1,hadronPtCut[i]+0.1);
    plotter->getHadronJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor,inputD0PtCut);
    hHadronJet[i] = (TH1D *)(plotter->getHadronJetCorrelation())->Clone(Form("hadronJetPt%i",i));
    for(int j=1;j<9;j++)
    {
      pair<int,int> inputCentralityBin(j+1,j+1);
      pair<double,double> inputPtCut(8,100);
      pair<double,double> inputD0PtCut(hadronPtCut[i]-0.1,hadronPtCut[i]+0.1);
      plotter->getHadronJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor,inputD0PtCut);
      hHadronJet[i]->Add(plotter->getHadronJetCorrelation());
    }
    // cHadronJet->cd(i+1);
    // hHadronJet[i]->Draw();
  }
  delete plotter;
  return hHadronJet[ptBin];
}

TH1D *plotD0JetInclusiveME(double dPtLow,double dPtHigh)
{
  int inputRebinFactor = 50;
  D0CorrPlotter *plotter = new D0CorrPlotter();
  plotter->init("root-files/d0CorrelationGlobal.root");

  TH1D *candJet[9];
  TH1D *bkgJet[9];

  TCanvas *cMB = new TCanvas();
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  int centBin[10] = {80,70,60,50,40,30,20,10,5,0};
  double sOverC[9];

  for(int i=0;i<9;i++)
  {
    pair<int,int> inputCentralityBin(i+1,i+1);
    pair<double,double> inputPtCut(8,100);
    pair<double,double> inputD0PtCut(dPtLow,dPtHigh);
    pair<double,double> inputHadronPtCut(dPtLow,dPtHigh);
    plotter->getJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor,inputD0PtCut);
    candJet[i] = plotter->getCandJetCorrelation();
    bkgJet[i] = plotter->getBkgJetCorrelation();
    sOverC[i] = plotter->getSBRatio();
  }

  plotter->init("root-files/mix-event.root");
  TH1D *candJetME[9];
  TH1D *bkgJetME[9];
  for(int i=0;i<9;i++)
  {
    pair<int,int> inputCentralityBin(i+1,i+1);
    pair<double,double> inputPtCut(8,100);
    pair<double,double> inputD0PtCut(dPtLow,dPtHigh);
    pair<double,double> inputHadronPtCut(dPtLow,dPtHigh);
    plotter->getJetCorrelation(inputCentralityBin,inputPtCut,inputRebinFactor,inputD0PtCut);
    candJetME[i] = plotter->getCandJetCorrelation();
    bkgJetME[i] = plotter->getBkgJetCorrelation();
    candJetME[i]->Scale(6.28/candJetME[i]->Integral()/candJetME[i]->GetBinWidth(3));
    bkgJetME[i]->Scale(6.28/bkgJetME[i]->Integral()/bkgJetME[i]->GetBinWidth(3));
  }
  TH1D *diff[9];
  for(int i=0;i<9;i++)
  {
    candJet[i]->Divide(candJetME[i]);
    bkgJet[i]->Divide(bkgJetME[i]);
    TH1D *bkgTemp = (TH1D *)bkgJet[i]->Clone("bkgTemp");
    bkgTemp->Scale(1-sOverC[i]);
    diff[i] = (TH1D *)candJet[i]->Clone(Form("diff%i",i));
    diff[i]->Add(bkgTemp,-1);
    diff[i]->Scale(1/sOverC[i]);
  }
  TH1D *corrS = (TH1D *)diff[3]->Clone("corrS");
  TH1D *corrB = (TH1D *)bkgJet[3]->Clone("corrB");
  TH1D *corrC = (TH1D *)candJet[3]->Clone("corrC");
  for(int i=4;i<9;i++)
  {
    corrS->Add(diff[i]);
    corrB->Add(bkgJet[i]);
    corrC->Add(candJet[i]);
  }
  TCanvas *c = new TCanvas();
  // corrS->Draw("pe");
  // corrS->SetLineColor(2);
  // corrS->SetMarkerColor(2);
  // corrC->Draw("same");
  // corrB->Draw("same");
  // corrB->SetLineColor(1);
  // corrC->SetLineColor(4);

  c->Divide(2,1);
  c->cd(1);
  corrS->Draw();
  plotHadronJetInclusive(dPtLow-1)->Draw("same");

  c->cd(2);
  TH1D *corrS_cp = (TH1D *)corrS->Clone("corrS_cp");
  corrS_cp->Add(plotHadronJetInclusive(dPtLow-1),-1);
  corrS_cp->Draw();

  // TLegend *legCorr = new TLegend(0.3,0.2,0.7,0.35);
  // legCorr->AddEntry(corrS,"D^{0} signal");
  // legCorr->AddEntry(corrB,"Side-band");
  // legCorr->AddEntry(corrC,"Candidates");
  // legCorr->Draw("same");
  delete plotter;
  return corrS_cp;
}
void plotCentrality()
{
  double dPtBin[6] = {1,2,3,4,5,10};
  TH1D *realS = plotD0JetInclusiveME(dPtBin[0],dPtBin[0+1]);
  for(int i=1;i<5;i++)
    realS->Add(plotD0JetInclusiveME(dPtBin[i],dPtBin[i+1]));

  TCanvas *leon = new TCanvas();
  realS->Draw();
}
