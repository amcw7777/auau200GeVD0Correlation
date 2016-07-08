#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <set>

#include "TFile.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TLorentzVector.h"

#include "phys_constants.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoD0EventMaker/StPicoD0Event.h"
#include "StPicoD0EventMaker/StKaonPion.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StCuts.h"
#include "../StPicoPrescales/StPicoPrescales.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
//////Refit include lib
#include "PhysicalConstants.h"
#include "StThreeVectorF.hh"
#include "StLorentzVectorD.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3D.h"
#include "TFile.h"
#include "StEvent/StDcaGeometry.h"
//
#include <vector>
//
#include <stdio.h>
#include <time.h>
#include <algorithm>
// Jet reconstruction lib
#include "TStarJetVectorContainer.h"
#include "TStarJetVector.h"
#include "fastjet/PseudoJet.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "StPicoD0AnaMaker.h"

ClassImp(StPicoD0AnaMaker)

  StPicoD0AnaMaker::StPicoD0AnaMaker(char const * name,char const * inputFilesList, 
      char const * outName,StPicoDstMaker* picoDstMaker,StRefMultCorr* grefmultCorrUtil): 
    StMaker(name),mPicoDstMaker(picoDstMaker),mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
    mOutFileName(outName), mInputFileList(inputFilesList),mOutputFile(NULL), mChain(NULL), mEventCounter(0){}

Int_t StPicoD0AnaMaker::Init()
{
  mPicoD0Event = new StPicoD0Event();

  mChain = new TChain("T");
  std::ifstream listOfFiles(mInputFileList.Data());
  if (listOfFiles.is_open())
  {
    std::string file;
    while (getline(listOfFiles, file))
    {
      LOG_INFO << "StPicoD0AnaMaker - Adding :" << file <<endm;
      mChain->Add(file.c_str());
      LOG_INFO<<" Entries = "<<mChain->GetEntries()<< endm; 
    }
  }
  else
  {
    LOG_ERROR << "StPicoD0AnaMaker - Could not open list of files. ABORT!" << endm;
    return kStErr;
  }

  mPrescales = new StPicoPrescales(mycuts::prescalesFilesDirectoryName);

  mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
  mChain->SetBranchAddress("dEvent", &mPicoD0Event);

  mOutputFile = new TFile(mOutFileName.Data(), "RECREATE");
//  phiFix = new TFile("phiFix.root");
  // phiWeight2D = (TH2D *)phiFix->Get("phiWeight")->Clone("phiWeight2D");

  mOutputFile->cd();

  float xbin[7] = {0,1,2,3,4,5,10};                                     
  float binMass[2001];
  float binCent[10];
  for(int i=0;i<2001;i++)                                               
    binMass[i] = 0.01*i; 
  for(int i=0;i<10;i++)
    binCent[i] = 1.0*i;
  massPt = new TH3F("massPt","",2000,binMass,6,xbin,9,binCent);
  massPt->Sumw2();

  vtxz = new TH1F("vtxz","",100,-10,10);
  dCount = new TH1F("dCount","dCount",10,0,10);
  hJetPt = new TH1F("hJetPt","",100,0,10);
  hHadronPt = new TH1F("hHadronPt","",100,0,10);
  candCount = new TH2F("candCount","",1,0.5,1.5,10,0,10);                   
  bkgCount = new TH2F("bkgCount","",1,0.5,1.5,10,0,10);
  hadronCount = new TH2F("hadronCount","",1,0.5,1.5,10,0,10);
  // phiRun = new TH2F("phiRun","",60018,15107000,15167018,1000,-1.*pi,1.*pi);
  // mOutputFile->cd();

  const int NDim = 3;// dPhi, centrality, pT
  const int NBinNumber[NDim] = {1000,10,100};
  const double XMin[NDim] = {-1.6,0,0};
  const double XMax[NDim] = {4.8,10,100};

  hD0JetCorrCand = new THnSparseD("hD0JetCorrCand","(1/N_{trig})(dN_{trig}/d#Delta#phi);#Delta #phi;centrality bin ; pT (GeV/c)",NDim,NBinNumber,XMin,XMax);
  hD0JetCorrBkg = new THnSparseD("hD0JetCorrBkg","(1/N_{trig})(dN_{trig}/d#Delta#phi);#Delta #phi;centrality bin;pT (GeV/c)",NDim,NBinNumber,XMin,XMax);
  hD0JetCorrCand->Sumw2();
  hD0JetCorrBkg->Sumw2();

  hD0HadronCorrCand = new THnSparseD("hD0HadronCorrCand","(1/N_{trig})(dN_{trig}/d#Delta#phi);#Delta #phi;centrality bin;pT (GeV/c)",NDim,NBinNumber,XMin,XMax);
  hD0HadronCorrBkg = new THnSparseD("hD0HadronCorrBkg","(1/N_{trig})(dN_{trig}/d#Delta#phi);#Delta #phi;centrality bin;pT (GeV/c)",NDim,NBinNumber,XMin,XMax);
  hD0HadronCorrCand->Sumw2();
  hD0HadronCorrBkg->Sumw2();

  const int NDimD0 = 4;// dPhi, centrality, pT
  const int NBinNumberD0[NDimD0] = {1000,10,100,5};
  const double XMinD0[NDimD0] = {-1.6,0,0,0};
  const double XMaxD0[NDimD0] = {4.8,10,100,5};
  hD0D0Corr = new THnSparseD("hD0D0CorrCand","(1/N_{trig})(dN_{trig}/d#Delta#phi);#Delta #phi;centrality bin;pT (GeV/c);corIndex",NDimD0,NBinNumberD0,XMinD0,XMaxD0);
  hD0D0Corr->Sumw2();

  hHadronJetCorr = new THnSparseD("hHadronJetCorr","(1/N_{trig})(dN_{trig}/d#Delta#phi);#Delta #phi;centrality bin;pT (GeV/c)",NDim,NBinNumber,XMin,XMax);
  hHadronJetCorr->Sumw2();


  jetPtPhi = new TH2F("jetPtPhi","jet-pt-phi;p_{T};#phi",1000,0,100,1000,-1.6,4.8);

  dcaCandJets = new TH3F("dcaCandJets","",10000,0,10,6,0,3.1416,10,0,10);
  dcaBkgJets = new TH3F("dcaBkgJets","",10000,0,10,6,0,3.1416,10,0,10);
  invMCandJets = new TH3F("invMCandJets","",1000,0,10,6,0,3.1416,10,0,10);
  invMBkgJets = new TH3F("invMBkgJets","",1000,0,10,6,0,3.1416,10,0,10);

  // -------------- USER VARIABLES -------------------------
  mGRefMultCorrUtil = new StRefMultCorr("grefmult");

  return kStOK;
}
//-----------------------------------------------------------------------------
StPicoD0AnaMaker::~StPicoD0AnaMaker()
{
  /*  */
  delete mGRefMultCorrUtil;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Finish()
{
  LOG_INFO << " StPicoD0AnaMaker - writing data and closing output file " <<endm;
  mOutputFile->cd();

  vtxz->Write();
  dCount->Write();
  hJetPt->Write();
  hHadronPt->Write();
  jetPtPhi->Write();

  massPt->Write();

  candCount->Write();
  bkgCount->Write();
  hadronCount->Write();

  hD0JetCorrCand->Write();
  hD0JetCorrBkg->Write();
  hD0HadronCorrCand->Write();
  hD0HadronCorrBkg->Write();
  hD0D0Corr->Write();
  hHadronJetCorr->Write();

  invMCandJets->Write();
  invMBkgJets->Write();

  mOutputFile->Close();
  delete mPrescales;

  return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoD0AnaMaker::Make()
{
  readNextEvent();
  if (!mPicoDstMaker)
  {
    LOG_WARN << " StPicoD0AnaMaker - No PicoDstMaker! Skip! " << endm;
    return kStWarn;
  }

  //StPicoDst const* picoDst = mPicoDstMaker->picoDst();
  picoDst = mPicoDstMaker->picoDst();

  if (!picoDst)
  {
    LOG_WARN << "StPicoD0AnaMaker - No PicoDst! Skip! " << endm;
    return kStWarn;
  }
  if(mPicoD0Event->runId() != picoDst->event()->runId() ||
      mPicoD0Event->eventId() != picoDst->event()->eventId())
  {
    LOG_ERROR <<" StPicoD0AnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!"<<endm;
    LOG_ERROR <<" StPicoD0AnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync."<<endm;
    exit(1);
  }

  // -------------- USER ANALYSIS -------------------------
  TClonesArray const * aKaonPion = mPicoD0Event->kaonPionArray();


  StThreeVectorF pVtx(-999.,-999.,-999.);
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  float field = event->bField();
  if(!(isGoodEvent()) )//minBias trigger requires
  {
    LOG_WARN << " Not Good Event! Skip! " << endm;
    return kStWarn;
  }
  if(event) {
    pVtx = event->primaryVertex();
  }
  vtxz->Fill(pVtx.z());
  if(!mGRefMultCorrUtil) {
    LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
    return kStWarn;
  }
  mGRefMultCorrUtil->init(picoDst->event()->runId());
  mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(),pVtx.z(),picoDst->event()->ZDCx()) ;
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  int centBin = 0;
  // int runBin = phiWeight2D->GetXaxis()->FindBin(event->runId());
  // phiWeight = (TH1D *)phiWeight2D->ProjectionY("phiWeight",runBin,runBin);
  // cout<<"runid = "<<event->runId()<<"\t run bin = "<<runBin<<endl;
  if(centrality>=7) centBin=1;
  else if(centrality>=4)  centBin=2;
  else centBin=3;
  fillHadronJetCorr();

  //Jet definition
  double jet_R = 0.2;
  fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm,jet_R);
  //Area defintion
  // double ghost_maxrap = 1.0;// fiducial cut for background estimation
  // fastjet::AreaDefinition areaDefinition(fastjet::active_area_explicit_ghosts,fastjet::GhostedAreaSpec(ghost_maxrap,1,0.01));
  // container for selected tracks

  double reweight = mGRefMultCorrUtil->getWeight();
  double pi = TMath::Pi();
  double fitmean[6] = {1.85921,1.8633,1.86403,1.86475,1.86252,1.86534};
  double fitsigma[6] = {0.018139,0.0139476,0.0158346,0.0169282,0.0199567,0.0189131};
  double d_counting = 0;
  for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
  {
    StKaonPion const* kp = (StKaonPion*)aKaonPion->At(idx);
    StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
    StPicoTrack const* pion = picoDst->track(kp->pionIdx());

    if (!isGoodTrack(kaon) || !isGoodTrack(pion)) continue;
    if (!isTpcPion(pion)) continue;
    bool tpcKaon = isTpcKaon(kaon,&pVtx);
    float kBeta = getTofBeta(kaon,&pVtx);
    bool tofAvailable = kBeta>0;
    bool tofKaon = tofAvailable && isTofKaon(kaon,kBeta);
    bool goodKaon = (tofAvailable && tofKaon) || (!tofAvailable && tpcKaon);
    if(!goodKaon) continue;
    int charge=isD0Pair(kp);
    if(-2==charge) continue;
    if(charge!=0)  d_counting++;
    // dcaD0Daughters->Fill(kp->pionDca(),pion->gPt());
    // dcaD0Daughters->Fill(kp->kaonDca(),kaon->gPt());
    set<unsigned int> dDaughters;
    dDaughters.insert(kp->kaonIdx());
    dDaughters.insert(kp->pionIdx());
    std::vector<fastjet::PseudoJet> selectedTracks;
    // std::vector<fastjet::PseudoJet,int> trackAdd;
    cout<<"number of tracks = "<<picoDst->numberOfTracks()<<endl;
    int hCount[10] = {0};
    for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
    {
      StPicoTrack const* mTrack = picoDst->track(i);
      if(!mTrack) continue;
      if(mTrack->nHitsFit()<20) continue;
      if(1.0*mTrack->nHitsFit()/mTrack->nHitsMax()<0.52) continue;
      if(mTrack->gPt()<0.2) continue;
      double trackDca = getDca(mTrack);
      if(trackDca>3)  continue;
      if(dDaughters.find(i) != dDaughters.end())  continue;
      double trackPx = mTrack->gMom(pVtx,field).x();   
      double trackPy = mTrack->gMom(pVtx,field).y();   
      double trackPz = mTrack->gMom(pVtx,field).z();   
      double trackM = 0.13957;
      if(isTpcKaon(mTrack,&pVtx))
        trackM = 0.493667;

      double trackE = sqrt( trackPx*trackPx + trackPy*trackPy + trackPz*trackPz + trackM);
      fastjet::PseudoJet pseudoJet(trackPx,trackPy,trackPz,trackE);
      selectedTracks.push_back(pseudoJet);
      // trackAdd[pseudoJet] = i;
      // std::pair<fastjet::PseudoJet,int> pairT(pseudoJet,i);
      // trackAdd.push_back(pairT);
    }
    fastjet::ClusterSequence cs(selectedTracks,jetDefinition);
    // double ptmin = 0.2;
    std::vector<fastjet::PseudoJet> mInclusiveJets;
    mInclusiveJets = cs.inclusive_jets();

    // for(unsigned int i=0;i<mInclusiveJets.size();i++)
    // {
    //   double jetPhi = mInclusiveJets[i].phi();
    //   if(jetPhi>pi)  jetPhi = 2.*pi-jetPhi;
    //   jetPtPhi->Fill(mInclusiveJets[i].pt(),jetPhi);
    // }


    float d0Pt = kp->pt();
    double d0Mass = kp->m();
    if(d0Pt>10) continue;
    if(d0Pt<3)  continue;
    int fitindex = 5;
    if(d0Pt<5)
      fitindex = static_cast<int>(d0Pt);
    int ptIdx = 5;
    if(kp->pt()<5)
      ptIdx = static_cast<int>(kp->pt());
    double mean = fitmean[ptIdx];
    double sigma = fitsigma[ptIdx];

    double reweight_eff = (efficiency[0][fitindex]/efficiency[centBin][fitindex]);
    reweight*=reweight_eff;
    if(charge!=0)
      massPt->Fill(kp->m(),d0Pt,centrality,reweight);
    bool isCand = d0Mass>mean-3*sigma && d0Mass<mean+3*sigma && charge!=0;
    bool isBkg = (d0Mass>mean-9*sigma && d0Mass<mean-4*sigma) ||
      (d0Mass<mean+9*sigma && d0Mass>mean+4*sigma) ||
      (d0Mass>mean-3*sigma && d0Mass<mean+3*sigma && charge==0);
    if(isCand)
      candCount->Fill(1,centrality);
    if(isBkg)
      bkgCount->Fill(1,centrality);

    cout<<"is tracks = "<<hCount[0]<<endl;
    cout<<"pass # hits = "<<hCount[1]<<endl;
    cout<<"pass # prob = "<<hCount[2]<<endl;
    cout<<"pass pT = "<<hCount[3]<<endl;
    cout<<"pass DCA = "<<hCount[4]<<endl;

    dCount->Fill(d_counting);
///// Loop to do D0-hadron correlation
    double hadronPtCut = 1.;
    double jetPtCut = 3.;
    for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
    {
      StPicoTrack const* mTrack = picoDst->track(i);
      if(!mTrack) continue;
      if(mTrack->nHitsFit()<20) continue;
      if(1.0*mTrack->nHitsFit()/mTrack->nHitsMax()<0.52) continue;
      if(mTrack->gPt()<hadronPtCut) continue;
      double trackDca = getDca(mTrack);
      if(trackDca>3)  continue;
      hHadronPt->Fill(mTrack->gPt());
      if(dDaughters.find(i) != dDaughters.end())  continue;
      double hadron_phi = mTrack->gMom(pVtx,field).phi();
      double trackPt = mTrack->gPt();

      double deltaPhi = (hadron_phi-kp->phi());
      if(deltaPhi<-0.5*pi)  deltaPhi += 2*pi;
      if(deltaPhi>1.5*pi)  deltaPhi -= 2*pi;
      double hadronFill[] = {deltaPhi,(double)centrality,trackPt};
      if(isCand)
        hD0HadronCorrCand->Fill(hadronFill,reweight);
      if(isBkg)
        hD0HadronCorrBkg->Fill(hadronFill,reweight);
    }
//////////////////////////////////////
///// Loop to do D0-jetcorrelation
//
    for(unsigned int i=0;i<mInclusiveJets.size();i++)
    {
      if(mInclusiveJets[i].pt()<jetPtCut)  continue;
      double deltaPhi = (mInclusiveJets[i].phi()-kp->phi());
      if(deltaPhi<-0.5*pi)  deltaPhi += 2*pi;
      if(deltaPhi>1.5*pi)  deltaPhi -= 2*pi;
      double jetFill[] = {deltaPhi,(double)centrality,mInclusiveJets[i].pt()};
      if(isCand)
        hD0JetCorrCand->Fill(jetFill,reweight);
      if(isBkg)
        hD0JetCorrBkg->Fill(jetFill,reweight);
      std::vector<fastjet::PseudoJet> constituents = mInclusiveJets[i].constituents();
      for(unsigned int iC=0;iC<constituents.size();iC++)
      {
         for(unsigned int jC = iC;jC<constituents.size();jC++)
         {
           fastjet::PseudoJet cPair = constituents[iC]+constituents[jC];
           if(isCand)
             invMCandJets->Fill(cPair.m(),deltaPhi,centrality,reweight);
           if(isBkg)
             invMBkgJets->Fill(cPair.m(),deltaPhi,centrality,reweight);
         }
      }
    }
///// Loop to do D0-jet correlation
//////////////////////////////////////
///// Loop to do D0-D0 correlation
    for (int jdx = idx; jdx < aKaonPion->GetEntries(); ++jdx)
    {
      StKaonPion const* kp_acc = (StKaonPion*)aKaonPion->At(jdx);
      StPicoTrack const* kaon_acc = picoDst->track(kp_acc->kaonIdx());
      StPicoTrack const* pion_acc = picoDst->track(kp_acc->pionIdx());

      if (!isGoodTrack(kaon_acc) || !isGoodTrack(pion_acc)) continue;
      if (!isTpcPion(pion_acc)) continue;
      bool tpcKaon_acc = isTpcKaon(kaon_acc,&pVtx);
      float kBeta_acc = getTofBeta(kaon_acc,&pVtx);
      bool tofAvailable_acc = kBeta_acc>0;
      bool tofKaon_acc = tofAvailable_acc && isTofKaon(kaon_acc,kBeta_acc);
      bool goodKaon_acc = (tofAvailable_acc && tofKaon_acc) || (!tofAvailable_acc && tpcKaon_acc);
      if(!goodKaon_acc) continue;
      int charge_acc=isD0Pair150(kp_acc);
      if(charge_acc==-2) continue;
      if(charge_acc*charge>=0) continue;

      float d0PtAcc = kp_acc->pt();
      double d0MassAcc = kp_acc->m();
      int ptIdxAcc = 5;
      if(kp_acc->pt()<5)
        ptIdxAcc = static_cast<int>(kp_acc->pt());
      double meanAcc = fitmean[ptIdxAcc];
      double sigmaAcc = fitsigma[ptIdxAcc];

      bool isCandAcc = d0MassAcc>meanAcc-3*sigmaAcc && d0MassAcc<meanAcc+3*sigmaAcc && charge!=0;
      bool isBkgAcc = (d0MassAcc>meanAcc-9*sigmaAcc && d0MassAcc<meanAcc-4*sigmaAcc) ||
        (d0MassAcc<meanAcc+9*sigmaAcc && d0MassAcc>meanAcc+4*sigmaAcc) ||
        (d0MassAcc>meanAcc-3*sigmaAcc && d0MassAcc<meanAcc+3*sigmaAcc && charge==0);
      int indexCor = 0;
      if(isCand && isCandAcc)  indexCor = 1;
      if(isCand && isBkgAcc)  indexCor = 2;
      if(isBkg && isCandAcc)  indexCor = 3;
      if(isBkg && isBkgAcc)  indexCor = 4;
      double deltaPhi = kp->phi()-kp_acc->phi();
      if(fabs(deltaPhi-0.)<0.000001)  continue;
      if(deltaPhi<-0.5*pi)  deltaPhi += 2*pi;
      if(deltaPhi>1.5*pi)  deltaPhi -= 2*pi;
      double d0Fill[] = {deltaPhi,(double)centrality,d0PtAcc,(double)indexCor};
      hD0D0Corr->Fill(d0Fill,reweight);
    }

  }
  return kStOK;
}
//-----------------------------------------------------------------------------

bool StPicoD0AnaMaker::isGoodPair(StKaonPion const* const kp) const
{
  if(!kp) return false;

  StPicoTrack const* kaon = mPicoDstMaker->picoDst()->track(kp->kaonIdx());
  StPicoTrack const* pion = mPicoDstMaker->picoDst()->track(kp->pionIdx());

  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  pVtx = event->primaryVertex();
  int charge = kaon->charge() * pion->charge();
  bool pairCuts = kp->m()>1.6 && kp->m()<2.1 &&
    charge==-1;

  return (isTpcKaon(kaon,&pVtx) && isTpcPion(pion) && 
      pairCuts);
}


int StPicoD0AnaMaker::isD0Pair(StKaonPion const* const kp) const // return -1:unlike sign; 1:like-sign 
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0061 &&
      kp->pionDca() > 0.0110 && kp->kaonDca() > 0.0103 &&
      kp->dcaDaughters() < 0.0084 && kp->decayLength()>0.0145;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0049 &&
      kp->pionDca() > 0.0111 && kp->kaonDca() > 0.0091 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0181;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0095 &&
      kp->dcaDaughters() < 0.0057 && kp->decayLength()>0.0212;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0038 &&
      kp->pionDca() > 0.0081 && kp->kaonDca() > 0.0079 &&
      kp->dcaDaughters() < 0.0050 && kp->decayLength()>0.0247;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0040 &&
      kp->pionDca() > 0.0062 && kp->kaonDca() > 0.0058 &&
      kp->dcaDaughters() < 0.0060 && kp->decayLength()>0.0259;  
  }

  int charge = (kaon->charge() * pion->charge()>0)?0:kaon->charge();
  if(pairCuts)
    return charge;
  else
    return -2;
}
/*
*/

int StPicoD0AnaMaker::isD0Pair150(StKaonPion const* const kp) const
{

  StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
  StPicoTrack const* pion = picoDst->track(kp->pionIdx());
  bool pairCuts = false;
  if(kp->pt()<1)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0072 &&
      kp->pionDca() > 0.0092 && kp->kaonDca() > 0.0105 &&
      kp->dcaDaughters() < 0.0077 && kp->decayLength()>0.0110;  
  }
  else if(kp->pt()<2)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0053 &&
      kp->pionDca() > 0.0078 && kp->kaonDca() > 0.0068 &&
      kp->dcaDaughters() < 0.0078 && kp->decayLength()>0.0168;  
  }
  else if(kp->pt()<3)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0047 &&
      kp->pionDca() > 0.0086 && kp->kaonDca() > 0.0080 &&
      kp->dcaDaughters() < 0.0074 && kp->decayLength()>0.0187;  
  }
  else if(kp->pt()<5)
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0042 &&
      kp->pionDca() > 0.0065 && kp->kaonDca() > 0.0066 &&
      kp->dcaDaughters() < 0.0068 && kp->decayLength()>0.0199;  
  }
  else 
  {
    pairCuts =  sin(kp->pointingAngle())*kp->decayLength() < 0.0062 &&
      kp->pionDca() > 0.0047 && kp->kaonDca() > 0.0041 &&
      kp->dcaDaughters() < 0.0066 && kp->decayLength()>0.0180;  
  }

  int charge = (kaon->charge() * pion->charge()>0)?0:kaon->charge();
  if(pairCuts)
    return charge;
  else
    return -2;
}

bool StPicoD0AnaMaker::isGoodEvent()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  return (isMBTrigger()) &&
    fabs(event->primaryVertex().z()) < mycuts::vz &&
    fabs(event->primaryVertex().z() - event->vzVpd()) < mycuts::vzVpdVz;
  //  return event->triggerWord() & mycuts::triggerWord;
}
bool StPicoD0AnaMaker::isMBTrigger()
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  return (event->isTrigger(450050) ||
      event->isTrigger(450060)||
      event->isTrigger(450005)||
      event->isTrigger(450015)||
      event->isTrigger(450025));
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodTrack(StPicoTrack const * const trk) const
{
  // Require at least one hit on every layer of PXL and IST.
  // It is done here for tests on the preview II data.
  // The new StPicoTrack which is used in official production has a method to check this
  return trk->gPt() > mycuts::minPt && trk->nHitsFit() >= mycuts::nHitsFit && trk->isHFTTrack();
  //return  trk->nHitsFit() >= mycuts::nHitsFit;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcPion(StPicoTrack const * const trk) const
{
  return fabs(trk->nSigmaPion()) < mycuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTpcKaon(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{
  return fabs(trk->nSigmaKaon()) < mycuts::nSigmaKaon;
  //      || tofKaon;
}
//-----------------------------------------------------------------------------
float StPicoD0AnaMaker::getTofBeta(StPicoTrack const * const trk, StThreeVectorF const* const pVtx) const
{

  int index2tof = trk->bTofPidTraitsIndex();

  float beta = std::numeric_limits<float>::quiet_NaN();

  if(index2tof >= 0)
  {
    StPicoBTofPidTraits *tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

    if(tofPid)
    {
      beta = tofPid->btofBeta();

      if (beta < 1e-4)
      {
        StThreeVectorF const btofHitPos = tofPid->btofHitPos();
        StPhysicalHelixD helix = trk->helix();

        float L = tofPathLength(pVtx, &btofHitPos, helix.curvature());
        float tof = tofPid->btof();
        if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
        else beta = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  return beta;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isTofKaon(StPicoTrack const * const trk, float beta) const
{
  bool tofKaon = false;

  if(beta>0)
  {
    double ptot = trk->dcaGeometry().momentum().mag();
    float beta_k = ptot/sqrt(ptot*ptot+M_KAON_PLUS*M_KAON_PLUS);
    tofKaon = fabs(1/beta - 1/beta_k) < mycuts::kTofBetaDiff ? true : false;
  }
  return tofKaon;
}
//-----------------------------------------------------------------------------
bool StPicoD0AnaMaker::isGoodHadron(StPicoTrack const * const trk) const
{
  return trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->pMom().pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}
bool StPicoD0AnaMaker::isGoodGlobalHadron(StPicoTrack const * const trk) const
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  float field = event->bField();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  pVtx = event->primaryVertex();
  return trk->nHitsFit() >= mycuts::nHitsFit &&fabs(trk->gMom(pVtx,field).pseudoRapidity())<1. && (1.0*trk->nHitsFit()/trk->nHitsMax())>0.52;
}
double StPicoD0AnaMaker::getDca(StPicoTrack const * const trk) const 
{
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  StThreeVectorF pVtx(-999.,-999.,-999.);
  pVtx = event->primaryVertex();
  StPhysicalHelixD hlx = trk->helix();
  double pathlength = hlx.pathLength( pVtx, false ); // do not scan periods
  return (hlx.at(pathlength)-pVtx).mag();
}

void StPicoD0AnaMaker::fillHadronJetCorr()
{
  StThreeVectorF pVtx(-999.,-999.,-999.);
  StPicoEvent *event = (StPicoEvent *)picoDst->event();
  float field = event->bField();
  pVtx = event->primaryVertex();
  int centrality  = mGRefMultCorrUtil->getCentralityBin9();
  double reweight = mGRefMultCorrUtil->getWeight();

  std::vector<fastjet::PseudoJet> selectedTracks;
  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* mTrack = picoDst->track(i);
    if(!mTrack) continue;
    if(mTrack->nHitsFit()<20) continue;
    if(1.0*mTrack->nHitsFit()/mTrack->nHitsMax()<0.52) continue;
    if(mTrack->gPt()<0.2) continue;
    double trackDca = getDca(mTrack);
    if(trackDca>3)  continue;
    double trackPx = mTrack->gMom(pVtx,field).x();   
    double trackPy = mTrack->gMom(pVtx,field).y();   
    double trackPz = mTrack->gMom(pVtx,field).z();   
    double trackM = 0.13957;
    if(isTpcKaon(mTrack,&pVtx))
      trackM = 0.493667;

    double trackE = sqrt( trackPx*trackPx + trackPy*trackPy + trackPz*trackPz + trackM);
    fastjet::PseudoJet pseudoJet(trackPx,trackPy,trackPz,trackE);
    selectedTracks.push_back(pseudoJet);
    // trackAdd[pseudoJet] = i;
    // std::pair<fastjet::PseudoJet,int> pairT(pseudoJet,i);
    // trackAdd.push_back(pairT);
  }
  double jet_R = 0.2;
  fastjet::JetDefinition jetDefinition(fastjet::antikt_algorithm,jet_R);
  fastjet::ClusterSequence cs(selectedTracks,jetDefinition);
  // double ptmin = 0.2;
  std::vector<fastjet::PseudoJet> mInclusiveJets;
  mInclusiveJets = cs.inclusive_jets();

  for(unsigned int i=0;i<picoDst->numberOfTracks();++i)
  {
    StPicoTrack const* mTrack = picoDst->track(i);
    if(!isGoodHadron(mTrack)) continue;
    if(mTrack->gPt()<3) continue;

    double hadronPhi = mTrack->gMom(pVtx,field).phi();
    hadronCount->Fill(1,centrality);

    for(unsigned int i=0;i<mInclusiveJets.size();i++)
    {
      if(mInclusiveJets[i].pt()<3)  continue;
      double deltaPhi = (mInclusiveJets[i].phi()-hadronPhi);
      if(deltaPhi<-0.5*pi)  deltaPhi += 2*pi;
      if(deltaPhi>1.5*pi)  deltaPhi -= 2*pi;
      double jetFill[] = {deltaPhi,(double)centrality,mInclusiveJets[i].pt()};
      hHadronJetCorr->Fill(jetFill,reweight);
    }
  }
}

