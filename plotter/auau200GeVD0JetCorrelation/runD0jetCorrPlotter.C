#include "d0jet_corr_plotter.cc"
void runD0jetCorrPlotter()
{
  int inputRebinFactor = 100;
  pair<int,int> inputCentralityBin(0,10);
  D0jetCorrPlotter *plotter = new D0jetCorrPlotter();
  plotter->init("root-files/dHCorr_Pt2_P16.root");
  // plotter->getCorrelation(inputRebinFactor);
  plotter->getCorrelation(inputCentralityBin,inputRebinFactor);
  double sOverC = plotter->getSBRatio();
  // plotter->plotCorrelation();
  TCanvas *c = new TCanvas();
  c->cd();
  // TH1D *candFar = plotter->getCandCorrelationFar();
  // TH1D *bkgFar = plotter->getBkgCorrelationFar();
  TH1D *candFar = plotter->getCandCorrelation();
  TH1D *bkgFar = plotter->getBkgCorrelation();
  candFar->Draw();
  // bkgFar->Scale(1-sOverC);
  bkgFar->Draw("same");
  candFar->SetLineColor(2);

  TCanvas *c2 = new TCanvas();
  // c2->Divide(2,1);
  // c2->cd(1);
  TH1D *diff = (TH1D *)candFar->Clone("diff");
  diff->Add(bkgFar,-1);
  diff->Scale(1/sOverC);
  diff->SetLineColor(4);
  // diff->Draw("same");
  diff->Draw();
  TLegend *leg = new TLegend(0.1,0.5,0.5,0.9);
  leg->AddEntry(bkgFar,"Bkg correlation");
  leg->AddEntry(candFar,"candidate correlation");
  // leg->AddEntry(diff,"signal correlation");
  // leg->Draw("same");
  // c2->cd(2);
  // TH1D *ratio = (TH1D *)candFar->Clone("ratio");
  // ratio->Divide(bkgFar);
  // ratio->Draw();




  // plotter->getPxCut(0.1);
}

