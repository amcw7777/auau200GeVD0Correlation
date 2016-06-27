#include "d0jet_corr_plotter.cc"
void plotCentrality()
{
  int inputRebinFactor = 50;
  D0jetCorrPlotter *plotter = new D0jetCorrPlotter();
  // plotter->init("root-files/dJCorrPt16.root");
  plotter->init("root-files/dHCorr_Pt3_P16.root");
  TCanvas *c = new TCanvas();
  c->Divide(3,2);
  TH1D *candFar[9];
  TH1D *bkgFar[9];
  TH1D *diff[9];
  for(int i=0;i<6;i++)
  {
    pair<int,int> inputCentralityBin(i+4,i+4);
    plotter->getCorrelation(inputCentralityBin,inputRebinFactor);
    candFar[i] = plotter->getCandCorrelation();
    bkgFar[i] = plotter->getBkgCorrelation();
    c->cd(i+1);
    // candFar[i]->Draw();
    // bkgFar[i]->Draw("same");
    // candFar[i]->SetLineColor(2);
    double sOverC = plotter->getSBRatio();
    cout<<"sOverC = "<<sOverC<<endl;
    bkgFar[i]->Scale(1-sOverC);
    diff[i] = (TH1D *)candFar[i]->Clone(Form("diff%i",i));
    diff[i]->Add(bkgFar[i],-1);
    diff[i]->Scale(1/sOverC);
    diff[i]->SetLineColor(4);
    diff[i]->Draw();
  }
  TH1D *corrS = (TH1D *)diff[0]->Clone("corrS");
  TH1D *corrB = (TH1D *)bkgFar[0]->Clone("corrB");
  TH1D *corrC = (TH1D *)bkgFar[0]->Clone("corrC");
  for(int i=1;i<6;i++)
  {
    corrS->Add(diff[i]);
    corrB->Add(bkgFar[i]);
    corrC->Add(candFar[i]);
  }
  TCanvas *c3 = new TCanvas();
  corrS->Draw("pe");
  // corrB->Draw("pe,same");
  // corrC->Draw("pe,same");
  corrS->SetLineColor(1);
  corrB->SetLineColor(2);
  corrC->SetLineColor(3);
    
}

