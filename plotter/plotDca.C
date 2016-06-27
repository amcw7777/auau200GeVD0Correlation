void plotDca()
{
  // TFile *file = new TFile("root-files/dcaPt.root");
  TFile *file = new TFile("root-files/dcaPt.root");
  TH2F *globalDcaTracks = (TH2F *)file->Get("dcaGlobalTracks")->Clone("globalDcaTracks");
  TH2F *primaryDcaTracks = (TH2F *)file->Get("dcaPrimaryTracks")->Clone("primaryDcaTracks");
  TH1D *globalDca[9];
  TH1D *primaryDca[9];
  TCanvas *c = new TCanvas();
  c->Divide(3,3);
  for(int i=0;i<9;i++)
  {
    globalDca[i] = globalDcaTracks->ProjectionX(Form("globalDca%i",i+1),i+1,i+1);
    primaryDca[i] = primaryDcaTracks->ProjectionX(Form("primaryDca%i",i+1),i+1,i+1);
    primaryDca[i]->SetLineColor(2);
    c->cd(i+1);
    gPad->SetLogy();
    globalDca[i]->Draw();
    primaryDca[i]->Draw("same");
  }
}



