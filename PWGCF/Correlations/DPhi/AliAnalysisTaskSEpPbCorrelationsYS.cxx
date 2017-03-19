/**************************************************************************************************                                                                           
 *      Leading Charged Track+VZERO and VZERO-VZERO Correlation.(Works for Real,Monte Carlo Data)  *                                                         
 *                                Yuko Sekiguchi                                                   *                                                                              
 *               Center for Nuclear Study (CNS) , University of Tokyo                              *                                                                              
 *                       Email:y_sekiguchi@cns.s.u-tokyo.ac.jp                                     *                                                                              
 **************************************************************************************************/
#include <TROOT.h>
#include <TList.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TRandom.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <TMap.h>
#include "AliLog.h"
#include "AliAnalysisManager.h"
//#include "AliAnalysisUtils.h"
#include "AliGenEventHeader.h"
#include "AliCFContainer.h"
#include "AliTHn.h"
#include "AliAODTrack.h"

#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVAD.h"
#include "AliESDAD.h"

#include "AliAODv0.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliCentrality.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliInputEventHandler.h"
#include "AliVParticle.h"
#include "AliMultiplicity.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliExternalTrackParam.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliAODcascade.h"

#include "AliMultSelection.h"

//#include "AliAODForwardMult.h"
//#include "AliForwardUtil.h"

#include "AliAnalysisTaskSEpPbCorrelationsYS.h"

ClassImp(AliAnalysisTaskSEpPbCorrelationsYS)
ClassImp(AliAssociatedTrackYS)
ClassImp(AliMixTrackYS)
ClassImp(AliAssociatedVZEROYS)

AliAnalysisTaskSEpPbCorrelationsYS::AliAnalysisTaskSEpPbCorrelationsYS():
  AliAnalysisTaskSE(),
  fcollisiontype     ("pPb"),
  fDataType          (kTRUE),
  fQA                (kTRUE),
  fAnaMode           ( "V0AV0C" ),
  fasso              ("Phi"),
  fPID               (kFALSE),
  fCentType          ( "ZNA" ),
  lCentrality        (0),
  bSign              (0),
  fZVertex           (10.),
  fOutputList        (0),
  fOutputList1       (0),
  fOutputList2       (0),
  fPIDResponse       (0),
  ffilterbit(5),
  fPtMin             (0.3),
  fEtaMax            (0.8),
  fEtaMinV0(0.8),
  fdcaDaughtersToPrimVtx(0.06),
  fdcaBetweenDaughters(1.0),
  fRadiMin(0.5),
  fRadiMax(100),
  fcutcTauLam(30),
  fcutcTauK0(20),
  fcosMinK0s(0.97),
  fcosMinLambda(0.995),
  fMaxnSigmaTPCV0(5),
  fclustermin(70),
  fratiocluster(0.8),
  fEtaMaxDaughter(0.8),
  fHistMass_K0s(0),
  fHistMass_Lambda(0),
  fHistMass_ALambda(0),
  fHist_V0QA(0),
  fHist_AP(0),
  fEvent             (0),
  lPrimaryBestVtx(0),
  fPrimaryZVtx(0),
  fvzero             (0),
  fPoolMgr           (0),
  fPoolMgr1           (0),
  poolmin            (0),
  poolmax            (0),
  fPoolMaxNEvents                              (2000),
  fPoolMinNTracks                              (50000),
  fMinEventsToMix                              (5),
  fNzVtxBins                                   (0),
  fNCentBins                                   (0),
  fMaxnSigmaTPCTOF   (3.),
  fHistzvertex    (0),
  fHistCentrality    (0),
  fHistPIDQA(0),
  fHist_vzeromult          (0),
  fHist_vzeromultEqweighted(0),
  fHist2dmult(0),
  fHistVZERO(0),
  fHist_Stat         (0),
  fHistPhiDTPCNSig(0),
  fHistPhiDTOFNSig(0),
  fHistPhiDTPCTOFNSig(0),
  fHistMass_PhiMeson(0),
  fHistMass_PhiMeson_MIX(0),
  fHist_PhiQA(0),
   fHistTriggerTrack  (0),
  fHistReconstTrack(0),
  fHistTriggerTrackMix(0),
  fHistReconstTrackMix(0)
  
{  
  for(Int_t iBin = 0; iBin < 100; iBin++){
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }
  for(Int_t i=0;i<3;i++){
    tPrimaryVtxPosition[i]=0;
    fHistPosNsig[i]=NULL;
    fHistNegNsig[i]=NULL;
  }
  
}

AliAnalysisTaskSEpPbCorrelationsYS::AliAnalysisTaskSEpPbCorrelationsYS(const char *name):
  AliAnalysisTaskSE(name),
  fcollisiontype     ("pPb"),
  fDataType          (kTRUE),
  fQA                (kTRUE),
  fAnaMode           ( "V0AV0C" ),
  fasso              ("Phi"),
  fPID               (kFALSE),
  fCentType          ( "ZNA" ),
  lCentrality        (0),
  bSign               (0),
  fZVertex           (10.),
  fOutputList        (0),
  fOutputList1       (0),
  fOutputList2       (0),
  fPIDResponse       (0),
  ffilterbit(5),
  fPtMin             (0.3),
  fEtaMax            (0.8),
  fEtaMinV0(0.8),
  fdcaDaughtersToPrimVtx(0.06),
  fdcaBetweenDaughters(1.0),
  fRadiMin(0.5),
  fRadiMax(100),
  fcutcTauLam(30),
  fcutcTauK0(20),
  fcosMinK0s(0.97),
  fcosMinLambda(0.995),
  fMaxnSigmaTPCV0(5),
  fclustermin(70),
  fratiocluster(0.8),
  fEtaMaxDaughter(0.8),
  fHistMass_K0s(0),
  fHistMass_Lambda(0),
  fHistMass_ALambda(0),
  fHist_V0QA(0),
  fHist_AP(0),
  fEvent             (0),
  lPrimaryBestVtx(0),
  fPrimaryZVtx(0),
  fvzero             (0),
  fPoolMgr           (0),  
  fPoolMgr1           (0),
  poolmin            (0),
  poolmax            (0),
  fPoolMaxNEvents                              (2000),
  fPoolMinNTracks                              (50000),
  fMinEventsToMix                              (5),
  fNzVtxBins                                   (0),
  fNCentBins                                   (0),
  fMaxnSigmaTPCTOF   (3.),
  fHistzvertex   (0),
  fHistCentrality    (0),
  fHistPIDQA(0),
  fHist_vzeromult          (0),
  fHist_vzeromultEqweighted(0),
  fHist2dmult(0),
  fHistVZERO(0),
  fHist_Stat         (0),
  fHistPhiDTPCNSig(0),
  fHistPhiDTOFNSig(0),
  fHistPhiDTPCTOFNSig(0),
  fHistMass_PhiMeson(0),
  fHistMass_PhiMeson_MIX(0),
  fHist_PhiQA(0),
  fHistTriggerTrack  (0),
  fHistReconstTrack(0),
  fHistTriggerTrackMix(0),
  fHistReconstTrackMix(0)
  
{
  for(Int_t iBin = 0; iBin < 100; iBin++){
    fZvtxBins[iBin] = 0.;
    fCentBins[iBin] = 0.;
  }
  for(Int_t i=0;i<3;i++){
    tPrimaryVtxPosition[i]=0;
    fHistPosNsig[i]=NULL;
    fHistNegNsig[i]=NULL;
  }
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TList::Class());
}

AliAnalysisTaskSEpPbCorrelationsYS::~AliAnalysisTaskSEpPbCorrelationsYS()
{
  if (fOutputList && fOutputList1 && fOutputList2 && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    delete fOutputList;
    delete fOutputList1;
    delete fOutputList2;
  }
}

void AliAnalysisTaskSEpPbCorrelationsYS::UserCreateOutputObjects()
{
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->SetName("global");
  DefineGeneralOutput();
  PostData(1, fOutputList);

  fOutputList1 = new TList();
  fOutputList1->SetOwner();
  fOutputList1->SetName("anahistos");
  DefineCorrOutput();
  DefineVZEROOutput();
  PostData(2, fOutputList1);

  fOutputList2 = new TList();
  fOutputList2 -> SetOwner();
  fOutputList2 -> SetName("QA");
  DefinedQAHistos();
  PostData(3,fOutputList2);
  
  fPoolMgr = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins, fNzVtxBins, fZvtxBins);
  if(!fPoolMgr) return;
  fPoolMgr->SetTargetValues(fPoolMinNTracks,0.1,5);
  
  fPoolMgr1 = new AliEventPoolManager(fPoolMaxNEvents, fPoolMinNTracks, fNCentBins, fCentBins, fNzVtxBins, fZvtxBins);
  if(!fPoolMgr1) return;
  fPoolMgr1->SetTargetValues(fPoolMinNTracks,0.1,5);
}

void AliAnalysisTaskSEpPbCorrelationsYS::DefineGeneralOutput(){
  fHist_Stat= new TH1F("fHist_Stat","Stat Histogram",11, -0.5, 10.5);
  fHist_Stat->GetXaxis()->SetBinLabel(1,"All Events");
  fHist_Stat->GetXaxis()->SetBinLabel(2,"Analyzed Events");
  fHist_Stat->GetXaxis()->SetBinLabel(3,"Trigger OK");
  fHist_Stat->GetXaxis()->SetBinLabel(4,"Vertex OK");
  fHist_Stat->GetXaxis()->SetBinLabel(5,"Centrality OK");
  fHist_Stat->GetXaxis()->SetBinLabel(6,"Pile-Up event rejected");
  fOutputList->Add(fHist_Stat);
  
  fHistzvertex=new TH1F("fHistzvertex", ";VZ;count", 60, -15, 15);
  fOutputList->Add(fHistzvertex);
  
  fHistCentrality=new TH1F("fHistCentrality", ";centrality;count", 20, 0, 100);
  fOutputList->Add(fHistCentrality);
  
  TTree *settingsTree   = new TTree("UEAnalysisSettings","Analysis Settings in UE estimation");
  settingsTree->Branch("fZVertex", &fZVertex,"fZVertex/D");
  settingsTree->Branch("fEtaMax",&fEtaMax,"fEtaMax/D");
  settingsTree->Branch("fPtMin",&fPtMin,"fPtMin/D");
  settingsTree->Branch("fMaxnSigmaTPCTOF",&fMaxnSigmaTPCTOF,"fMaxnSigmaTPCTOF/D");
  
  //V0 Particle
  settingsTree->Branch("fEtaMinV0",&fEtaMinV0,"fEtaMinV0/D");
  settingsTree->Branch("fdcaDaughtersToPrimVtx",&fdcaDaughtersToPrimVtx,"fdcaDaughtersToPrimVtx/D");
  settingsTree->Branch("fdcaBetweenDaughters",& fdcaBetweenDaughters," fdcaBetweenDaughters/D");
  settingsTree->Branch("fRadiMin",&fRadiMin,"fRadiMin/D");
  settingsTree->Branch("fRadiMax",&fRadiMax,"fRadiMax/D");
  settingsTree->Branch("fcutcTauK0",&fcutcTauK0,"fcutcTauK0");
  settingsTree->Branch("fcutcTauLam",&fcutcTauLam,"fcutcTauLam");
  settingsTree->Branch("fcosMinK0s",&fcosMinK0s,"fcosMinK0s");
  settingsTree->Branch("fcosMinLambda",&fcosMinLambda,"fcosMinLambda");
  settingsTree->Branch("fMaxnSigmaTPCV0",&fMaxnSigmaTPCV0,"fMaxnSigmaTPCV0");
  //Phi
  settingsTree->Branch("ffilterbit",&ffilterbit,"ffilterbit/I");
  

  //  settingsTree->Branch("fanamode",&fAnaMode,"fAnaMode/B");
  //  settingsTree->Branch("fanalysisasso",&fanalysisasso,"fanalysisasso/I");
  //settingsTree->Branch("fanalysiscent",&fanalysiscent,"fanalysiscent/I");
  settingsTree->Fill();
  fOutputList->Add(settingsTree);

  return;
}
void AliAnalysisTaskSEpPbCorrelationsYS::DefineVZEROOutput(){
  fHist_vzeromult=new TH2F("fHist_vzeromult","fHist_vzeromult",64,-0.5,63.5,500,0,500);
  fOutputList1->Add(fHist_vzeromult);
  
  fHist_vzeromultEqweighted=new TH2F("fHist_vzeromultEqweighted","fHist_vzeromultEqweighted",64,-0.5,63.5,500,0,500);
  fOutputList1->Add(fHist_vzeromultEqweighted);
  
  fHist2dmult=new TH3F("fHist2dmult","fHist2dmult",64,-0.5,63.5,500,0,500,500,0,500);
  fOutputList1->Add(fHist2dmult);
  
  const  Int_t nVZEROBins[3]={10,8,15};
  Double_t binning_eta_vzero[11]={-3.7,-3.2,-2.7,-2.2,-1.7,0.,2.8,3.4,3.9,4.5,5.1};
  Double_t binning_phi_vzero[9]={0.,0.7853,1.5707,2.3561,3.1415,3.9269,4.7123,5.4977,6.2831};
  Double_t binning_cent[16]={0.,1.,2.,3.,4.,5.,10.,20., 30., 40., 50., 60., 70., 80., 90., 100.1};
  
  const Int_t ipidBin[4]={11,40,72,15};
  Double_t binning_pt_lead[12]={0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,  8.0};
  Double_t binning_eta[41]={-1.,-0.95,-0.9,-0.85,-0.8,-0.75,-0.7,-0.65,-0.6,-0.55,-0.5,-0.45,-0.4,-0.35,-0.3,-0.25,-0.2,-0.15,-0.1,-0.05,0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0};
  Double_t binning_dphi[73]={-1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0, 0.087266, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.047198, 1.134464, 1.221730, 1.308997, 1.396263, 1.483530, 1.570796, 1.658063, 1.745329, 1.832596, 1.919862, 2.007129, 2.094395, 2.181662, 2.268928, 2.356194, 2.443461, 2.530727, 2.617994, 2.705260, 2.792527, 2.879793, 2.967060, 3.054326, 3.141593, 3.228859, 3.316126, 3.403392, 3.490659, 3.577925, 3.665191, 3.752458, 3.839724, 3.926991, 4.014257, 4.101524, 4.188790, 4.276057, 4.363323, 4.450590, 4.537856, 4.625123, 4.712389};
  
  if(fQA){
    fHistVZERO=new AliTHn("fHistVZERO","fHistVZERO",1,3,nVZEROBins);
    fHistVZERO->SetBinLimits(0,binning_eta_vzero);
    fHistVZERO->SetBinLimits(1,binning_phi_vzero);
    fHistVZERO->SetBinLimits(2,binning_cent);
    fOutputList1->Add(fHistVZERO);
    
    fHistPIDQA=new AliTHn("fHistPIDQA", "fHistPIDQA", 1, 4,ipidBin);
    fHistPIDQA->SetBinLimits(0,binning_pt_lead);
    fHistPIDQA->SetBinLimits(1,binning_eta);
    fHistPIDQA->SetBinLimits(2,binning_dphi);
    fHistPIDQA->SetBinLimits(3,binning_cent);
    fHistPIDQA->SetVarTitle(0,"pt");
    fHistPIDQA->SetVarTitle(1,"eta");
    fHistPIDQA->SetVarTitle(2,"phi");
    fHistPIDQA->SetVarTitle(3,"centrality");
    fOutputList1->Add(fHistPIDQA);
  }

  return;
  
}

void AliAnalysisTaskSEpPbCorrelationsYS::DefinedQAHistos(){
  for(Int_t i=0;i<3;i++){
    fHistPosNsig[i]=new TH2D(Form("fHistPosNsig_%d",i),"fHistPosNsig",80,0.,8.,200,-10.,10);
    fHistNegNsig[i]=new TH2D(Form("fHistNegNsig_%d",i),"fHistNegNsig",80,0.,8.,200,-10.,10);
    fOutputList2->Add(fHistPosNsig[i]);
    fOutputList2->Add(fHistNegNsig[i]);
  }
  
  
  fHistPhiDTPCNSig=new TH2D("fHistPhiDTPCNSig","fHistPhiDTPCNSig",150,0.,15.,200,-10.,10);
  fOutputList2->Add(fHistPhiDTPCNSig);
  fHistPhiDTOFNSig=new TH2D("fHistPhiDTOFNSig","fHistPhiDTOFNSig",150,0.,15.,200,-10.,10);
  fOutputList2->Add(fHistPhiDTOFNSig);
  fHistPhiDTPCTOFNSig=new TH2D("fHistPhiDTPCTOFNSig","fHistPhiDTPCTOFNSig",150,0.,15.,200,-10.,10);
  fOutputList2->Add(fHistPhiDTPCTOFNSig);
  
  Int_t nBins = 400;
  Double_t mphiMin=1.02-0.1;
  Double_t mphiMax=1.02+0.1;

  Int_t nCentralityBins  = 20;
  Double_t centBins1[16]={0.,1.,2.,3.,4.,5.,10.,20., 30., 40., 50., 60., 70., 80., 90., 100.0};
  const Double_t* centralityBins = centBins1;
  
  Int_t nPtBinsV0 = 150; 
  const Double_t PtBinsV0[12] = {0, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 8.0, 15.0};

  Double_t mk0sMin=0.5-0.1;
  Double_t mk0sMax=0.5+0.1;
  const Int_t spBins[3] = {nBins, nPtBinsV0, nCentralityBins};
  //v0
  const Double_t spMink0s[3] = {mk0sMin, PtBinsV0[0], centralityBins[0]};
  const Double_t spMaxk0s[3] = {mk0sMax, PtBinsV0[11], centralityBins[15]};
  fHistMass_K0s=new THnSparseF("fHistMass_K0s","mass for K0s",3, spBins, spMink0s, spMaxk0s);
  fOutputList2->Add(fHistMass_K0s);

  Double_t mlambdaMin=1.15-0.1;
  Double_t mlambdaMax=1.15+0.1;
  const Double_t spMinLambda[3] = {mlambdaMin, PtBinsV0[0], centralityBins[0]};
  const Double_t spMaxLambda[3] = {mlambdaMax, PtBinsV0[11], centralityBins[15]};
  fHistMass_Lambda=new THnSparseF("fHistMass_Lambda","mass for Lambda",3, spBins, spMinLambda, spMaxLambda);
  fOutputList2->Add(fHistMass_Lambda);
  fHistMass_ALambda=new THnSparseF("fHistMass_ALambda","mass for Anti Lambda",3, spBins, spMinLambda, spMaxLambda);
  fOutputList2->Add(fHistMass_ALambda);
  
  const Double_t spMinPhi[3] = {mphiMin, PtBinsV0[0], centralityBins[0]};
  const Double_t spMaxPhi[3] = {mphiMax, PtBinsV0[11], centralityBins[15]};
    //Phimeson
  fHistMass_PhiMeson=new THnSparseF("fHistMass_PhiMeson","mass for phi meson",3, spBins, spMinPhi, spMaxPhi);
  fOutputList2->Add(fHistMass_PhiMeson);

  fHistMass_PhiMeson_MIX=new THnSparseF("fHistMass_PhiMeson_MIX","mass for phi meson of mixed events",3, spBins, spMinPhi, spMaxPhi);
  fOutputList2->Add(fHistMass_PhiMeson_MIX);
  
  const Int_t spBinsQA[4] = {40, 72, 7,20};
  const Double_t spMinPhiQA[4] = {-1. ,0, -0.5,0.};
  const Double_t spMaxPhiQA[4] = {1.  ,TMath::TwoPi(), 6.5,100.0};
  fHist_PhiQA=new THnSparseF("fHist_PhiQA","QA for Phimeson", 4, spBinsQA,spMinPhiQA,spMaxPhiQA);
  fOutputList2->Add(fHist_PhiQA);

  fHist_V0QA=new THnSparseF("fHist_V0QA","QA for V0 particle", 4, spBinsQA,spMinPhiQA,spMaxPhiQA);
  fOutputList2->Add(fHist_V0QA);
  
  fHist_AP=new TH2D("fHist_AP","fHist_AP", 200, -1, 1, 200, 0, 0.4);
  fOutputList2->Add(fHist_AP);

  return;
}

void AliAnalysisTaskSEpPbCorrelationsYS::DefineCorrOutput(){

  Double_t binning_pt_assoc[12]={0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,  8.0};
  Double_t binning_pt_lead[12]={0.3, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0,  8.0};
  Double_t binning_cent[16]={0.,1.,2.,3.,4.,5.,10.,20., 30., 40., 50., 60., 70., 80., 90., 100.1};
  Double_t binning_deta[49]={-2.4, -2.3, -2.2, -2.1, -2.0, -1.9, -1.8, -1.7, -1.6, -1.5, -1.4, -1.3, -1.2, -1.1, -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,2.1, 2.2, 2.3, 2.4};
  
  Double_t binning_dphi[73]={-1.570796, -1.483530, -1.396263, -1.308997, -1.221730, -1.134464, -1.047198, -0.959931, -0.872665, -0.785398, -0.698132, -0.610865, -0.523599, -0.436332, -0.349066, -0.261799, -0.174533, -0.087266, 0.0, 0.087266, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, 0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.047198, 1.134464, 1.221730, 1.308997, 1.396263, 1.483530, 1.570796, 1.658063, 1.745329, 1.832596, 1.919862, 2.007129, 2.094395, 2.181662, 2.268928, 2.356194, 2.443461, 2.530727, 2.617994, 2.705260, 2.792527, 2.879793, 2.967060, 3.054326, 3.141593, 3.228859, 3.316126, 3.403392, 3.490659, 3.577925, 3.665191, 3.752458, 3.839724, 3.926991, 4.014257, 4.101524, 4.188790, 4.276057, 4.363323, 4.450590, 4.537856, 4.625123, 4.712389};
  
  Int_t nCFSteps = 1;
  
  const Int_t nEvtVars = 2;  
  const Int_t iEvtBin[2]={11,15};
  
  fHistTriggerTrack= new AliTHn("fHistTriggerTrack", "fHistTriggerTrack", nCFSteps, nEvtVars, iEvtBin);
  fHistTriggerTrack->SetBinLimits(0,binning_pt_lead);
  fHistTriggerTrack->SetBinLimits(1,binning_cent);
  fHistTriggerTrack->SetVarTitle(0,"leading p_{T} GeV/c");
  fHistTriggerTrack->SetVarTitle(1,"centrality");
  
  fHistTriggerTrackMix= new AliTHn("fHistTriggerTrackMix", "fHistTriggerTrackMix", nCFSteps, nEvtVars, iEvtBin);
  fHistTriggerTrackMix->SetBinLimits(0,binning_pt_lead);
  fHistTriggerTrackMix->SetBinLimits(1,binning_cent);
  fHistTriggerTrackMix->SetVarTitle(0,"leading p_{T} GeV/c");
  fHistTriggerTrackMix->SetVarTitle(1,"centrality");
  
  fOutputList1->Add(fHistTriggerTrack);  
  fOutputList1->Add(fHistTriggerTrackMix);  
  
  const Int_t nTrackVars=5;
  const Int_t iTrackBin[5]={48,11,11,15,72};   
    
    if(fasso=="V0"||fasso=="Phi"||(fasso=="hadron"&&fPID))nCFSteps=7;
    if(fasso=="hadron"&&!fPID) nCFSteps=1;
    fHistReconstTrack= new AliTHn("fHistReconstTrack", "fHistReconstTrack", nCFSteps, nTrackVars,iTrackBin);
    fHistReconstTrack->SetBinLimits(0,binning_deta);
    fHistReconstTrack->SetBinLimits(1,binning_pt_assoc);
    fHistReconstTrack->SetBinLimits(2,binning_pt_lead);
    fHistReconstTrack->SetBinLimits(3,binning_cent);
    fHistReconstTrack->SetBinLimits(4,binning_dphi);
    fHistReconstTrack->SetVarTitle(0,"#Delta#eta");
    fHistReconstTrack->SetVarTitle(1,"p_{T} GeV/c");
    fHistReconstTrack->SetVarTitle(2,"leading p_{T} GeV/c");
    fHistReconstTrack->SetVarTitle(3,"centrality");
    fHistReconstTrack->SetVarTitle(4,"#Delta#phi");
    fHistReconstTrackMix= new AliTHn("fHistReconstTrackMix", "fHistReconstTrackMix", nCFSteps, nTrackVars,iTrackBin);
    fHistReconstTrackMix->SetBinLimits(0,binning_deta);
    fHistReconstTrackMix->SetBinLimits(1,binning_pt_assoc);
    fHistReconstTrackMix->SetBinLimits(2,binning_pt_lead);
    fHistReconstTrackMix->SetBinLimits(3,binning_cent);
    fHistReconstTrackMix->SetBinLimits(4,binning_dphi);
    fHistReconstTrackMix->SetVarTitle(0,"#Delta#eta");
    fHistReconstTrackMix->SetVarTitle(1,"p_{T} GeV/c");
    fHistReconstTrackMix->SetVarTitle(2,"leading p_{T} GeV/c");
    fHistReconstTrackMix->SetVarTitle(3,"centrality");
    fHistReconstTrackMix->SetVarTitle(4,"#Delta#phi");
    fOutputList1->Add(fHistReconstTrack);
    fOutputList1->Add(fHistReconstTrackMix);
    return;
}
void AliAnalysisTaskSEpPbCorrelationsYS::UserExec(Option_t *)
{
  
  fHist_Stat->Fill(0);
  
  AliAnalysisManager   *mgr      = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inEvMain = (AliInputEventHandler*)(mgr->GetInputEventHandler());
  if (!inEvMain) return;
  
  fPIDResponse = inEvMain->GetPIDResponse();
  if(!fPIDResponse) return;
  

  fHist_Stat->Fill(1);
  
  fEvent = dynamic_cast<AliAODEvent*>(inEvMain->GetEvent());
  //  fEvent = dynamic_cast<AliESDEvent*>(inEvMain->GetEvent());
  if(!fEvent) {
      AliWarning("ERROR: fEvent not available \n");
      return;
  }
  
  //Physics Selection
  UInt_t maskIsSelected = inEvMain->IsEventSelected();
  Bool_t isSelected=kFALSE;
  isSelected = ((maskIsSelected & AliVEvent::kINT7)== AliVEvent::kINT7); 
  if (!isSelected)    return;
  fHist_Stat->Fill(2);
  
  //Primary Vllllertex 
  lPrimaryBestVtx = fEvent->GetPrimaryVertex();
  if (!lPrimaryBestVtx) return;
  Int_t nTracksPrim= lPrimaryBestVtx->GetNContributors();
  if(nTracksPrim<1) return;
  if ( (TMath::Abs(lPrimaryBestVtx->GetZ())) >= fZVertex) return ;
  tPrimaryVtxPosition[0] = lPrimaryBestVtx->GetX();
  tPrimaryVtxPosition[1] = lPrimaryBestVtx->GetY();
  tPrimaryVtxPosition[2] = lPrimaryBestVtx->GetZ();  
  fHistzvertex->Fill(tPrimaryVtxPosition[2]);
  
  fPrimaryZVtx= lPrimaryBestVtx->GetZ();  
  
  fHist_Stat->Fill(3);     
  
  bSign = 0.;
  bSign = (InputEvent()->GetMagneticField()>0) ?  1 : -1;
  
  //Multiplicity Object
  
  AliMultSelection *multSelection = (AliMultSelection*)fEvent->FindListObject("MultSelection");
  if(!multSelection){
      AliWarning("AliMultSelection could not be found in the aod event list of objects");
  }
  
  
  if(fCentType=="V0A"){
    lCentrality=multSelection->GetMultiplicityPercentile("V0A");
    Int_t qual = multSelection->GetEvSelCode();
    if(qual == 199 ) lCentrality=-999;
  }
  
  //Centrality 
  /*
  AliCentrality *centobj = 0;
  lCentrality = 0;
  centobj = fEvent->GetCentrality();
  lCentrality = centobj->GetCentralityPercentile(fCentType);
  if(!centobj) lCentrality=-1.;
  */
  if(lCentrality<0.||lCentrality>100.-0.0000001) return;
  
  Double_t * CentBins = fCentBins;
  poolmin    = CentBins[0];
  poolmax    = CentBins[fNCentBins];
  fHist_Stat->Fill(4);

  //Pile up
  /*
  AliAnalysisUtils *fUtils=new AliAnalysisUtils;
  if(fUtils->IsFirstEventInChunk(fEvent)) return;
  if(!fUtils->IsVertexSelected2013pA(fEvent)) return;
  if(fUtils->IsPileUpEvent(fEvent))return;
  fHist_Stat->Fill(5);
  */
  //AnalysisPart 
  MakeAna();
 
  fHistCentrality->Fill(lCentrality);
  
  PostData(1,fOutputList);
  PostData(2,fOutputList1);
  PostData(3,fOutputList2);
}
  
void AliAnalysisTaskSEpPbCorrelationsYS::Terminate(Option_t *)
{
  //  AliInfo(Form("Number of Correlation Entries======================%f",fNEntries));
}

void AliAnalysisTaskSEpPbCorrelationsYS::MakeAna()
{
  TObjArray *selectedTracksLeading=new TObjArray;
  selectedTracksLeading->SetOwner(kTRUE);
  TObjArray *selectedTracksAssociated=new TObjArray;
  selectedTracksAssociated->SetOwner(kTRUE);
  
  // Leading Particle  
  selectedTracksLeading=GetAcceptedTracksLeading(fEvent);
  
  //associated Particle
  if(fasso=="Phi")selectedTracksAssociated=GetAcceptedTracksAssociated(fEvent);
  if(fasso=="V0")selectedTracksAssociated=GetAcceptedV0Tracks(fEvent);
  if(fasso=="hadron")selectedTracksAssociated=GetAcceptedTracksLeading(fEvent);

  /*
  //associated Tracks
  fvzero = fEvent->GetVZEROData();
  Double_t eta_min;
  Double_t eta_max;
  Double_t eta_ave;
  Double_t phi_vzero;
  Double_t mult_vzero;
  Double_t vzeroqa[3];
  Double_t mult_vzero_eq;
  
  for(Int_t imod=0;imod<64;imod++)
    {
      eta_min=fvzero->GetVZEROEtaMin(imod);
      eta_max=fvzero->GetVZEROEtaMax(imod);
      phi_vzero=fvzero->GetVZEROAvgPhi(imod);
      mult_vzero=fvzero->GetMultiplicity(imod);
      mult_vzero_eq=fEvent->GetVZEROEqMultiplicity(imod);
      eta_ave=(eta_min+eta_max)/2.;
      
      fHist_vzeromult->Fill(imod,mult_vzero);
      fHist_vzeromultEqweighted->Fill(imod,mult_vzero_eq);
      fHist2dmult->Fill(imod,mult_vzero_eq,mult_vzero);
      
      vzeroqa[0]=eta_ave;
      vzeroqa[1]=phi_vzero;
      vzeroqa[2]=lCentrality;
      //       cout<<vzeroqa[0]<<" "<<vzeroqa[1]<<" "<<mult_vzero_eq<<" "<<mult_vzero<<endl;
      if(fQA) fHistVZERO->Fill(vzeroqa,0,(Double_t)mult_vzero_eq);
  
      if(fAnaMode=="V0AV0C" && imod>31) selectedTracksLeading->Add(new AliAssociatedVZEROYS(mult_vzero_eq,eta_ave,phi_vzero,0.0,0,0));//V0A is used as Leading 
      
      if(fAnaMode=="TPCV0A" || fAnaMode=="SPDV0A") {
	if(imod>31) selectedTracksAssociated->Add(new AliAssociatedVZEROYS(mult_vzero_eq,eta_ave,phi_vzero,0.0,0,0));
      }else if(fAnaMode=="TPCV0C" || fAnaMode=="V0AV0C" || fAnaMode=="SPDV0C"){
	if(imod<32) selectedTracksAssociated->Add(new AliAssociatedVZEROYS(mult_vzero_eq,eta_ave,phi_vzero,0.0,0,0));
      }
  
    }  
  */
  
    FillCorrelationTracks(lCentrality,selectedTracksLeading,selectedTracksAssociated, fHistTriggerTrack,fHistReconstTrack,kFALSE,0.02,0.8,bSign,0); 
    FillCorrelationTracksMixing(lCentrality,lPrimaryBestVtx->GetZ(),poolmax,poolmin,selectedTracksLeading,selectedTracksAssociated, fHistTriggerTrackMix,fHistReconstTrackMix,kFALSE,0.02,0.8,bSign,0);
 
  selectedTracksLeading->Clear();
  delete selectedTracksLeading;
  selectedTracksAssociated->Clear();
  delete selectedTracksAssociated;
}

TObjArray *AliAnalysisTaskSEpPbCorrelationsYS::GetAcceptedTracksLeading(AliAODEvent*fAOD){
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  Int_t  nTracks= fAOD->GetNumberOfTracks();
  Double_t pidqa[4];  
  for (Int_t i = 0; i < nTracks; i++)
   {   AliAODTrack* aodTrack = dynamic_cast<AliAODTrack *>(fAOD->GetTrack(i));
      if (!aodTrack) continue;
      if (!IsAcceptedTrack(aodTrack)) continue;
      if(aodTrack->Charge()==0) continue;
      Double_t nSigmaKaonTPC   = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kKaon);
      Double_t nSigmaPionTPC   = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kPion);
      Double_t nSigmaProtonTPC   = fPIDResponse->NumberOfSigmasTPC(aodTrack, AliPID::kProton);
      Double_t nSigmaKaonTOF   = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kKaon);
      Double_t nSigmaPionTOF   = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kPion);
      Double_t nSigmaProtonTOF   = fPIDResponse->NumberOfSigmasTOF(aodTrack, AliPID::kProton);

      Double_t d2nsigmakaon=nSigmaKaonTPC*nSigmaKaonTPC+nSigmaKaonTOF*nSigmaKaonTOF;
      Double_t d2nsigmapion=nSigmaPionTPC*nSigmaPionTPC+nSigmaPionTOF*nSigmaPionTOF;
      Double_t d2nsigmaproton=nSigmaProtonTPC*nSigmaProtonTPC+nSigmaProtonTOF*nSigmaProtonTOF;

      Bool_t fPIDTOF=kTRUE;
      if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, fAOD->GetTrack(i))==0) fPIDTOF=kFALSE;
      else fPIDTOF=kTRUE;
      
      Double_t nSigmaKaonTOFTPC;
      Double_t nSigmaPionTOFTPC;
      Double_t nSigmaProtonTOFTPC;

      if(fPIDTOF &&  aodTrack->Pt()>0.5){
	nSigmaKaonTOFTPC=TMath::Sqrt(d2nsigmakaon);
	nSigmaPionTOFTPC=TMath::Sqrt(d2nsigmapion);
	nSigmaProtonTOFTPC=TMath::Sqrt(d2nsigmaproton);
      }else{
	nSigmaKaonTOFTPC=TMath::Abs(nSigmaKaonTPC);
	nSigmaPionTOFTPC=TMath::Abs(nSigmaPionTPC);
	nSigmaProtonTOFTPC=TMath::Abs(nSigmaProtonTPC);
      }
      
      Int_t SpPID=-999;
      if( (nSigmaKaonTOFTPC<fMaxnSigmaTPCTOF) && (nSigmaKaonTOFTPC<nSigmaPionTOFTPC) && (nSigmaKaonTOFTPC<nSigmaProtonTOFTPC))  SpPID=1;
      if( (nSigmaPionTOFTPC<fMaxnSigmaTPCTOF) && (nSigmaPionTOFTPC<nSigmaKaonTOFTPC) && (nSigmaPionTOFTPC<nSigmaProtonTOFTPC))  SpPID=2;
      if( (nSigmaProtonTOFTPC<fMaxnSigmaTPCTOF) && (nSigmaProtonTOFTPC<nSigmaKaonTOFTPC) && (nSigmaProtonTOFTPC<nSigmaPionTOFTPC)) SpPID=3;
      
      tracks->Add(new AliAssociatedTrackYS(aodTrack->Charge(),aodTrack->Eta(),aodTrack->Phi(),aodTrack->Pt(),aodTrack->GetID(),-999,-999,SpPID,1));
   }
  return tracks;
}

TObjArray *AliAnalysisTaskSEpPbCorrelationsYS::GetAcceptedTracksAssociated(AliAODEvent*fAODEvent){
  TObjArray *tracks=new TObjArray;
  tracks->SetOwner(kTRUE);
  TObjArray *dtrack=new TObjArray;
  dtrack->SetOwner(kTRUE);
  
  Int_t nTracks= fAODEvent->GetNumberOfTracks();
  for (Int_t i = 0; i < nTracks-1; i++){
    AliAODTrack* aodTrack1 = dynamic_cast<AliAODTrack *>(fAODEvent->GetTrack(i));                                                               
    if (!aodTrack1) continue;			
    Double_t nSigmaKaonTPC_phi1   = fPIDResponse->NumberOfSigmasTPC(aodTrack1, AliPID::kKaon); 
    Double_t nSigmaKaonTOF_phi1   = fPIDResponse->NumberOfSigmasTOF(aodTrack1, AliPID::kKaon); 
    Double_t nSigmaPionTPC_phi1   = fPIDResponse->NumberOfSigmasTPC(aodTrack1, AliPID::kPion);
    Double_t nSigmaPionTOF_phi1   = fPIDResponse->NumberOfSigmasTOF(aodTrack1, AliPID::kPion);
    Double_t nSigmaProtonTPC_phi1 = fPIDResponse->NumberOfSigmasTPC(aodTrack1, AliPID::kProton);
    Double_t nSigmaProtonTOF_phi1 = fPIDResponse->NumberOfSigmasTOF(aodTrack1, AliPID::kProton);
    Double_t d2sigmaphi1kaontpctof=nSigmaKaonTPC_phi1*nSigmaKaonTPC_phi1+nSigmaKaonTOF_phi1*nSigmaKaonTOF_phi1;
    Double_t d2sigmaphi1piontpctof=nSigmaPionTPC_phi1*nSigmaPionTPC_phi1+nSigmaPionTOF_phi1*nSigmaPionTOF_phi1;
    Double_t d2sigmaphi1protontpctof=nSigmaProtonTPC_phi1*nSigmaProtonTPC_phi1+nSigmaProtonTOF_phi1*nSigmaProtonTOF_phi1;
    Bool_t fPIDTOF_phi1;  
    
    if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, aodTrack1)==0)  fPIDTOF_phi1=kFALSE;
    else fPIDTOF_phi1=kTRUE;
    
    Double_t nSigmaKaonTOFTPC_phi1;
    Double_t nSigmaPionTOFTPC_phi1;
    Double_t nSigmaProtonTOFTPC_phi1;
    if(fPIDTOF_phi1 && aodTrack1->Pt()>0.5){
      nSigmaKaonTOFTPC_phi1=TMath::Sqrt(d2sigmaphi1kaontpctof);
      nSigmaPionTOFTPC_phi1=TMath::Sqrt(d2sigmaphi1piontpctof);
      nSigmaProtonTOFTPC_phi1=TMath::Sqrt(d2sigmaphi1protontpctof);
    }else{
      nSigmaKaonTOFTPC_phi1=TMath::Abs(nSigmaKaonTPC_phi1);
      nSigmaPionTOFTPC_phi1=TMath::Abs(nSigmaPionTPC_phi1);
      nSigmaProtonTOFTPC_phi1=TMath::Abs(nSigmaProtonTPC_phi1);
    }
    if(!IsAcceptedPhiDaughterTrack(aodTrack1)) continue;       
   
    
    fHistPhiDTPCNSig->Fill(aodTrack1->Pt(),nSigmaKaonTPC_phi1);
    fHistPhiDTOFNSig->Fill(aodTrack1->Pt(),nSigmaKaonTOF_phi1);
    fHistPhiDTPCTOFNSig->Fill(aodTrack1->Pt(),TMath::Sqrt(d2sigmaphi1kaontpctof));
   
    Bool_t isKaon1=kFALSE;
    //if(nSigmaKaonTOFTPC_phi1<5.0 && nSigmaKaonTOFTPC_phi1<nSigmaPionTOFTPC_phi1 && nSigmaKaonTOFTPC_phi1<nSigmaProtonTOFTPC_phi1) isKaon1=kTRUE;
    if(TMath::Abs(nSigmaKaonTPC_phi1)<5.0 && TMath::Abs(nSigmaKaonTPC_phi1)<TMath::Abs(nSigmaPionTPC_phi1) && TMath::Abs(nSigmaKaonTPC_phi1)<TMath::Abs(nSigmaProtonTPC_phi1)) isKaon1=kTRUE;
    if(!isKaon1) continue;
    
    dtrack->Add(new AliMixTrackYS(aodTrack1->Charge(),aodTrack1->Eta(),aodTrack1->Phi(),aodTrack1->Pt(),aodTrack1->Px(),aodTrack1->Py(),aodTrack1->Pz()));
    
    for (Int_t j = i+1; j < nTracks; j++)  {
      AliAODTrack* aodTrack2 = dynamic_cast<AliAODTrack *>(fAODEvent->GetTrack(j));                                                                   
      if (!aodTrack2) continue;     
      Double_t nSigmaKaonTPC_phi2   = fPIDResponse->NumberOfSigmasTPC(aodTrack2, AliPID::kKaon);                                                          
      Double_t nSigmaKaonTOF_phi2   = fPIDResponse->NumberOfSigmasTOF(aodTrack2, AliPID::kKaon);                                                          
      Double_t nSigmaPionTPC_phi2   = fPIDResponse->NumberOfSigmasTPC(aodTrack2, AliPID::kPion);
      Double_t nSigmaPionTOF_phi2   = fPIDResponse->NumberOfSigmasTOF(aodTrack2, AliPID::kPion);
      Double_t nSigmaProtonTPC_phi2 = fPIDResponse->NumberOfSigmasTPC(aodTrack2, AliPID::kProton);
      Double_t nSigmaProtonTOF_phi2 = fPIDResponse->NumberOfSigmasTOF(aodTrack2, AliPID::kProton);
      Double_t d2sigmaphi2kaontpctof=nSigmaKaonTPC_phi2*nSigmaKaonTPC_phi2+nSigmaKaonTOF_phi2*nSigmaKaonTOF_phi2;
      Double_t d2sigmaphi2piontpctof=nSigmaPionTPC_phi2*nSigmaPionTPC_phi2+nSigmaPionTOF_phi2*nSigmaPionTOF_phi2;
      Double_t d2sigmaphi2protontpctof=nSigmaProtonTPC_phi2*nSigmaProtonTPC_phi2+nSigmaProtonTOF_phi2*nSigmaProtonTOF_phi2;
      Bool_t fPIDTOF_phi2;

      if(fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, aodTrack2)==0)  fPIDTOF_phi2=kFALSE;
      else fPIDTOF_phi2=kTRUE;

      Double_t nSigmaKaonTOFTPC_phi2;
      Double_t nSigmaPionTOFTPC_phi2;
      Double_t nSigmaProtonTOFTPC_phi2;
      if(fPIDTOF_phi2 && aodTrack2->Pt()>0.5){
	nSigmaKaonTOFTPC_phi2=TMath::Sqrt(d2sigmaphi2kaontpctof);
	nSigmaPionTOFTPC_phi2=TMath::Sqrt(d2sigmaphi2piontpctof);
	nSigmaProtonTOFTPC_phi2=TMath::Sqrt(d2sigmaphi2protontpctof);
      }else{
	nSigmaKaonTOFTPC_phi2=TMath::Abs(nSigmaKaonTPC_phi2);
	nSigmaPionTOFTPC_phi2=TMath::Abs(nSigmaPionTPC_phi2);
	nSigmaProtonTOFTPC_phi2=TMath::Abs(nSigmaProtonTPC_phi2);
      }

      Bool_t isKaon2=kFALSE;
      //if(nSigmaKaonTOFTPC_phi2<5.0 && nSigmaKaonTOFTPC_phi2<nSigmaPionTOFTPC_phi2 && nSigmaKaonTOFTPC_phi2<nSigmaProtonTOFTPC_phi2) isKaon2=kTRUE;
      if(TMath::Abs(nSigmaKaonTPC_phi2)<5.0 && TMath::Abs(nSigmaKaonTPC_phi2)<TMath::Abs(nSigmaPionTPC_phi2) && TMath::Abs(nSigmaKaonTPC_phi2)<TMath::Abs(nSigmaProtonTPC_phi2)) isKaon2=kTRUE;
      if(!isKaon2) continue;                                                                                              
      
      if(!IsAcceptedPhiDaughterTrack(aodTrack2)) continue;
      if(aodTrack1->GetID()==aodTrack2->GetID()) continue;                                                                  
      if(aodTrack1->Charge()==aodTrack2->Charge()) continue;                                                                         
      
      Double_t mass_phi;                           
      Double_t px_phi = aodTrack1->Px()+aodTrack2->Px();                                                                                               
      Double_t py_phi = aodTrack1->Py()+aodTrack2->Py();                                                                                            
      Double_t pz_phi = aodTrack1->Pz()+aodTrack2->Pz();  
      Double_t p_phi = TMath::Sqrt(px_phi*px_phi+py_phi*py_phi+pz_phi*pz_phi);                                                                          
      Double_t phi_phi = atan2(py_phi,px_phi);
      if(phi_phi<0.) phi_phi += 2 * TMath::Pi();                                                                                                                  
      Double_t px1 = aodTrack1->Px();
      Double_t py1 = aodTrack1->Py();
      Double_t pz1 = aodTrack1->Pz();
      Double_t E1 = TMath::Sqrt(px1*px1+py1*py1+pz1*pz1+0.493677*0.493677);
      Double_t px2 = aodTrack2->Px();
      Double_t py2 = aodTrack2->Py();
      Double_t pz2 = aodTrack2->Pz();
      Double_t E2 = TMath::Sqrt(px2*px2+py2*py2+pz2*pz2+0.493677*0.493677);
      mass_phi = TMath::Sqrt((E1+E2)*(E1+E2)-(px1+px2)*(px1+px2)-(py1+py2)*(py1+py2)-(pz1+pz2)*(pz1+pz2));
      Double_t pt_phi = TMath::Sqrt((px1+px2)*(px1+px2)+(py1+py2)*(py1+py2));
      
      // if(pt_phi<0.5) continue;      
      
      Double_t eta_phi = 0.5 * log((p_phi+pz_phi)/(p_phi-pz_phi));                                                                                       
      
      Double_t PhiQA[3]={mass_phi, pt_phi, lCentrality};
      fHistMass_PhiMeson->Fill(PhiQA);
      if(TMath::Abs(eta_phi)>0.8) continue;
     
      Int_t SpPhi=-999.;
      
      const Double_t fPhiMass =  TDatabasePDG::Instance()->GetParticle(333)->Mass();
      if(fabs(mass_phi-fPhiMass)<0.006)  SpPhi=0; //same for step
      if(mass_phi>0.99 &&mass_phi<fPhiMass-0.006 ) SpPhi=1; //same as step
      if(mass_phi>fPhiMass+0.006 && mass_phi<fPhiMass+0.018) SpPhi=2; //same as step
      if(mass_phi>fPhiMass+0.018 && mass_phi<fPhiMass+0.030) SpPhi=3; //same as step
      if(mass_phi>fPhiMass+0.030 && mass_phi<fPhiMass+0.042) SpPhi=4; //same as step
      if(mass_phi>fPhiMass+0.042 && mass_phi<fPhiMass+0.054) SpPhi=5; //same as step
      if(mass_phi>fPhiMass+0.054 && mass_phi<fPhiMass+0.066) SpPhi=6; //same as step
      
      tracks->Add(new AliAssociatedTrackYS(0,eta_phi,phi_phi,pt_phi,-999,aodTrack1->GetID(),aodTrack2->GetID(),SpPhi,1));
      
      Double_t spPhiQA[4]={eta_phi, phi_phi, (Float_t)SpPhi,lCentrality};
      fHist_PhiQA->Fill(spPhiQA);
    
    }
  }

  //Mixed Event
  Double_t  pvxMix=fPrimaryZVtx;
  Double_t counterMix=0;
  AliEventPool* pool = fPoolMgr1->GetEventPool(lCentrality, pvxMix);
  if (!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", lCentrality, pvxMix));  
  if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks  || pool->GetCurrentNEvents() > fMinEventsToMix)
    {
      Int_t nMix = pool->GetCurrentNEvents();
      for (Int_t jMix=0; jMix<nMix; jMix++){
        TObjArray* mixEvents = pool->GetEvent(jMix);
	for(Int_t i=0;i<dtrack->GetEntriesFast();i++){
	  AliMixTrackYS* dtrack1=(AliMixTrackYS*)dtrack->At(i);
	  if(!dtrack1) continue;
	  Double_t pdx1 = dtrack1->Px();
	  Double_t pdy1 = dtrack1->Py();
	  Double_t pdz1 = dtrack1->Pz();
	  Double_t Ed1 = TMath::Sqrt(pdx1*pdx1+pdy1*pdy1+pdz1*pdz1+0.493677*0.493677);
	  counterMix++;
	  for(Int_t j=0;j<mixEvents->GetEntriesFast();j++){
	    AliMixTrackYS*dtrack2=(AliMixTrackYS*)mixEvents->At(j);
	    if(!dtrack2)continue;
	    if(dtrack1->Charge()==dtrack2->Charge()) continue;
	    Double_t pdx2 = dtrack2->Px();
	    Double_t pdy2 = dtrack2->Py();
	    Double_t pdz2 = dtrack2->Pz();
	    Double_t Ed2 = TMath::Sqrt(pdx2*pdx2+pdy2*pdy2+pdz2*pdz2+0.493677*0.493677);
	    Double_t mass_phi_mix = TMath::Sqrt((Ed1+Ed2)*(Ed1+Ed2)-(pdx1+pdx2)*(pdx1+pdx2)-(pdy1+pdy2)*(pdy1+pdy2)-(pdz1+pdz2)*(pdz1+pdz2));
	    Double_t pt_phi_mix = TMath::Sqrt((pdx1+pdx2)*(pdx1+pdx2)+(pdy1+pdy2)*(pdy1+pdy2));
	    // if(pt_phi_mix<0.5) continue;
	    Double_t px_phi_mix =  pdx1+pdx2;
	    Double_t py_phi_mix = pdy1+pdy2;
	    Double_t pz_phi_mix = pdz1+pdz2;
	    Double_t p_phi_mix = TMath::Sqrt(px_phi_mix*px_phi_mix+py_phi_mix*py_phi_mix+pz_phi_mix*pz_phi_mix);
	    Double_t phi_phi_mix = atan2(py_phi_mix,px_phi_mix);
	    if(phi_phi_mix<0.) phi_phi_mix += 2 * TMath::Pi();
	    Double_t eta_phi_mix = 0.5 * log((p_phi_mix+pz_phi_mix)/(p_phi_mix-pz_phi_mix));  
	    if(TMath::Abs(eta_phi_mix)>0.8) continue;   
	    Double_t PhiQAMix[3]={mass_phi_mix, pt_phi_mix, lCentrality};
	    fHistMass_PhiMeson_MIX->Fill(PhiQAMix,1./(Double_t)nMix);
	  }
	}
      }
    }

  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  
  for (Int_t i=0; i<dtrack->GetEntriesFast(); i++)
    {
      AliMixTrackYS*particle=(AliMixTrackYS*) dtrack->At(i);
      tracksClone->Add(new AliMixTrackYS(particle->Charge(),
				    particle->Eta(),
				    particle->Phi(),
				    particle->Pt(),
				    particle->Px(),
				    particle->Py(),
				    particle->Pz()
					  )
		       );
    }

  pool->UpdatePool(tracksClone);
  
  return tracks;
  
}

TObjArray *AliAnalysisTaskSEpPbCorrelationsYS::GetAcceptedV0Tracks(AliAODEvent*fAODEvent){
  TObjArray *tracks = new TObjArray;
  tracks->SetOwner(kTRUE);
  //V0Particles
  Int_t nv0s = fAODEvent->GetNumberOfV0s();
  for (Int_t iv0 = 0; iv0 < nv0s; iv0++)
    {
      AliAODv0* aodv0 = dynamic_cast<AliAODv0 *>(fAODEvent->GetV0(iv0));
      if (!aodv0) {
        AliError(Form("ERROR: Could not retrieve aodv0 %d", iv0));
        continue;
      }
      if(!IsAcceptedV0(aodv0)) continue; 
    
  //c*tau
  Double_t lDVx = aodv0->GetSecVtxX();
  Double_t lDVy = aodv0->GetSecVtxY();
  Double_t lDVz = aodv0->GetSecVtxZ();
  const Double_t kLambdaMass =  TDatabasePDG::Instance()->GetParticle(3122)->Mass();
  const Double_t kK0Mass =  TDatabasePDG::Instance()->GetParticle(311)->Mass();
  
  //Double_t cutcTauLam  = fcutctau*7.89;//c*mean life time=2.632*2.9979 cm                                                                           
  //Double_t cutcTauK0   = fcutctau*2.68; //c*mean life time=0.8954*2.9979 cm   
  Double_t lV0DecayLength = TMath::Sqrt(TMath::Power(lDVx -  tPrimaryVtxPosition[0],2) + TMath::Power(lDVy-  tPrimaryVtxPosition[1],2) + TMath::Power(lDVz -  tPrimaryVtxPosition[2],2 ));
  Double_t lPV0 = TMath::Sqrt((aodv0->Pt())*(aodv0->Pt())+(aodv0->Pz())*(aodv0->Pz()));
  Double_t lcTauLam = (lV0DecayLength*kLambdaMass)/lPV0;
  Double_t lcTauK0 = (lV0DecayLength*kK0Mass)/lPV0;
  //life time < 20(K0s), 30(lambda)
  Bool_t cutK0ctau = (lcTauK0<fcutcTauK0);
  Bool_t cutLambdactau = (lcTauLam<fcutcTauLam);
  Bool_t cutAntiLambdactau = (lcTauLam<fcutcTauLam);
  
  //cpa 
  Double_t cpa=aodv0->CosPointingAngle(fAODEvent->GetPrimaryVertex());
  Bool_t cpaK0s= cpa>fcosMinK0s;
  Bool_t cpaLambda=cpa>fcosMinLambda;
      
  AliAODTrack *myTrackPos=dynamic_cast<AliAODTrack *>(aodv0->GetDaughter(0)); //The first dauther track, which should be positive      
  AliAODTrack *myTrackNeg=dynamic_cast<AliAODTrack *>(aodv0->GetDaughter(1)); //The second dauther track, which should be negative
  
  if (!myTrackPos || !myTrackNeg) {
    Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
    continue;
  }
      
  if(!IsAcceptedDaughterTrack(myTrackPos)||!IsAcceptedDaughterTrack(myTrackNeg)) continue;

  if( (myTrackPos->Charge() == 1) && (myTrackNeg->Charge() == -1)){
    myTrackPos = myTrackPos;
    myTrackNeg = myTrackNeg;
  }
      
  if( (myTrackPos->Charge() == -1) && (myTrackNeg->Charge() == 1 )){
    myTrackPos = myTrackNeg;
    myTrackNeg = myTrackPos;
  }
      
  //PID Cut
  //Positive tracks
  Double_t nSigmaPosPionTPC = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kPion);
  Double_t nSigmaPosPionTOF = fPIDResponse->NumberOfSigmasTOF(myTrackPos, AliPID::kPion);
  Double_t nSigmaPosKaonTPC = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kKaon);
  Double_t nSigmaPosKaonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackPos, AliPID::kKaon);
  Double_t nSigmaPosProtonTPC = fPIDResponse->NumberOfSigmasTPC(myTrackPos, AliPID::kProton);
  Double_t nSigmaPosProtonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackPos, AliPID::kProton);
  //negative tracks
  Double_t nSigmaNegPionTPC = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kPion);
  Double_t nSigmaNegPionTOF = fPIDResponse->NumberOfSigmasTOF(myTrackNeg, AliPID::kPion);
  Double_t nSigmaNegKaonTPC = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kKaon);
  Double_t nSigmaNegKaonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackNeg, AliPID::kKaon);
  Double_t nSigmaNegProtonTPC = fPIDResponse->NumberOfSigmasTPC(myTrackNeg, AliPID::kProton);
  Double_t nSigmaNegProtonTOF = fPIDResponse->NumberOfSigmasTOF(myTrackNeg, AliPID::kProton);
  Bool_t bpPion = TMath::Abs(nSigmaPosPionTPC) <= fMaxnSigmaTPCV0;
  Bool_t bpProton = TMath::Abs(nSigmaPosProtonTPC) <= fMaxnSigmaTPCV0;
  Bool_t bnPion = TMath::Abs(nSigmaNegPionTPC) <= fMaxnSigmaTPCV0;
  Bool_t bnProton = TMath::Abs(nSigmaNegProtonTPC) <= fMaxnSigmaTPCV0;
  
  fHistPosNsig[0]->Fill(myTrackPos->Pt(),nSigmaPosPionTPC);
  fHistPosNsig[1]->Fill(myTrackPos->Pt(),nSigmaPosKaonTPC);
  fHistPosNsig[2]->Fill(myTrackPos->Pt(),nSigmaPosProtonTPC);
  fHistNegNsig[0]->Fill(myTrackNeg->Pt(),nSigmaNegPionTPC);
  fHistNegNsig[1]->Fill(myTrackNeg->Pt(),nSigmaNegKaonTPC);
  fHistNegNsig[2]->Fill(myTrackNeg->Pt(),nSigmaNegProtonTPC);
  
  Bool_t cutK0Pid = ( bpPion && bnPion);
  Bool_t cutLambdaPid = (bpProton && bnPion);
  Bool_t cutAntiLambdaPid = (bpPion && bnProton);
  
  //Armentros podoranski cut
  fHist_AP->Fill(aodv0->AlphaV0(),aodv0->PtArmV0());
  
  Bool_t k0APcut = (aodv0->PtArmV0()>(TMath::Abs(0.2*aodv0->AlphaV0())));
  //Mass cut
  Double_t lInvMassLambda = aodv0->MassLambda();
  Double_t lInvMassK0=aodv0->MassK0Short();
  Double_t lInvMassAntiLambda=aodv0->MassAntiLambda();
 
  Bool_t cutMassLambda = ((lInvMassLambda>1.05) && (lInvMassLambda<1.25));
  Bool_t cutMassAntiLambda = ((lInvMassAntiLambda>1.05) && (lInvMassAntiLambda<1.25));
  Bool_t cutMassK0 = (lInvMassK0>0.4) && (lInvMassK0<0.6);
  
  Bool_t IDK0s=cutK0Pid && cutK0ctau && cpaK0s && k0APcut && (!cutMassLambda) && (!cutMassAntiLambda);
  Bool_t IDLa= cutLambdaPid && cutLambdactau && cpaLambda && (!cutMassK0);
  Bool_t IDALa=  cutAntiLambdaPid && cutAntiLambdactau && cpaLambda && (!cutMassK0);
  
  Double_t K0sMassQA[3]={lInvMassK0,aodv0->Pt(),lCentrality};
  if(IDK0s)  fHistMass_K0s->Fill(K0sMassQA);
  Double_t LambdaMassQA[3]={lInvMassLambda,aodv0->Pt(),lCentrality};
  if(IDLa) fHistMass_Lambda->Fill(LambdaMassQA);
  Double_t ALambdaMassQA[3]={lInvMassAntiLambda,aodv0->Pt(),lCentrality};
  if(IDALa) fHistMass_ALambda->Fill(ALambdaMassQA);
  
  Int_t SpV0=-999;
  if(IDK0s && (fabs(lInvMassK0-kK0Mass)<0.01))  SpV0=0; //same for step
  if(IDK0s && (lInvMassK0>0.40)&&(lInvMassK0<0.44))  SpV0=1; //same as step
  if(IDK0s && (lInvMassK0>0.44)&&(lInvMassK0<kK0Mass-0.01))  SpV0=2; //same as step
  if(IDK0s && (lInvMassK0>kK0Mass+0.01)&&(lInvMassK0<0.56))  SpV0=3; //same as step
  if(IDK0s && (lInvMassK0>0.56)&&(lInvMassK0<0.60))  SpV0=4; //same as step
  if((IDLa || IDALa) && (fabs(lInvMassLambda-kLambdaMass)<0.005))  SpV0=5; //same as step
  if((IDLa || IDALa) && (lInvMassLambda>kLambdaMass+0.005)&&(lInvMassLambda<1.25))  SpV0=6; //same as step
  Double_t spV0QA[4]={aodv0->Eta(),aodv0->Phi(),(Float_t)SpV0,lCentrality};
  fHist_V0QA->Fill(spV0QA);
  
  tracks->Add(new AliAssociatedTrackYS(aodv0->Charge(),aodv0->Eta(),aodv0->Phi(),aodv0->Pt(),aodv0->GetID(),myTrackPos->GetID(),myTrackNeg->GetID(),SpV0,1));
    }
  
  return tracks;

}

Bool_t AliAnalysisTaskSEpPbCorrelationsYS::IsAcceptedPhiDaughterTrack(const AliAODTrack  * aodTrack){
  if(!aodTrack->TestFilterMask(BIT(ffilterbit))) return kFALSE;
  if(aodTrack->Pt()<0.15) return kFALSE;
  if(fabs(aodTrack->Eta())>0.8)  return kFALSE;
  if(aodTrack->Charge()==0) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsYS::IsAcceptedDaughterTrack(const AliAODTrack  * itrack)
{
  if(TMath::Abs(itrack->Eta())>fEtaMaxDaughter)return kFALSE;
  if (!itrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  Float_t nCrossedRowsTPC = itrack->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fclustermin) return kFALSE;
  Int_t findable=itrack->GetTPCNclsF();
  if (findable <= 0) return kFALSE;
  if (nCrossedRowsTPC/findable < fratiocluster) return kFALSE;
  if (itrack->GetKinkIndex(0)>0) return kFALSE;
  // if(itrack->Pt()<fPtMinDaughter) return kFALSE;                                                                                                                     
  return kTRUE;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsYS::IsAcceptedTrack(const AliAODTrack* aodTrack)
{
  if(!aodTrack) return kFALSE;
  //  if(!aodTrack->TestFilterMask(BIT(5))) return kFALSE; // standard cut with tight DCA cut                                                                     
  if(!aodTrack->TestFilterMask(BIT(5))) return kFALSE; // only tpc cut                                                                                    
  /*
  if (!aodTrack->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;
  Float_t nCrossedRowsTPC =aodTrack->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < 70) return kFALSE;
  Int_t findable=aodTrack->GetTPCNclsF();
  if (findable <= 0)return kFALSE;
  if (nCrossedRowsTPC/findable < 0.8) return kFALSE;
  */
  if(aodTrack->Pt()<fPtMin) return kFALSE;
  if(abs(aodTrack->Eta())>fEtaMax) return kFALSE;
  return kTRUE;
}

Bool_t AliAnalysisTaskSEpPbCorrelationsYS::IsAcceptedV0(const AliAODv0* aodV0){
  /*
    Pseudo rapidity V0 < 0.8
    daughter tracks DCA to P >0.06cm
    DCA between daughters <1 sigma
    V0 2D decay radius from 0.5cm to 200cm  
  */
  
  if(!aodV0) return kFALSE;
  if (aodV0->GetOnFlyStatus()) return kFALSE;// ofline reconstraction only
  
  Double_t leta_v0=aodV0->PseudoRapV0();
  if(TMath::Abs(leta_v0)>fEtaMinV0) return kFALSE;
  
  //DCA to primary vertex
  Float_t xyn=aodV0->DcaNegToPrimVertex();
  if (TMath::Abs(xyn)<fdcaDaughtersToPrimVtx) return kFALSE;
  Float_t xyp=aodV0->DcaPosToPrimVertex();
  if (TMath::Abs(xyp)<fdcaDaughtersToPrimVtx) return kFALSE;
  
  //DCA between dauther tracks
  Double_t dca=aodV0->DcaV0Daughters();
  if (dca>fdcaBetweenDaughters) return kFALSE;
  
  //Fiducial volume 
  Double_t xyz[3];
  aodV0->GetSecondaryVtx(xyz);
  Double_t r2=xyz[0]*xyz[0] + xyz[1]*xyz[1];
  if (r2<fRadiMin*fRadiMin) return kFALSE;
  if (r2>fRadiMax*fRadiMax) return kFALSE;

  return kTRUE;
}

void AliAnalysisTaskSEpPbCorrelationsYS::FillCorrelationTracks(Double_t centrality,
                                                                       TObjArray *triggerArray,
                                                                       TObjArray*selectedTrackArray,
                                                                       AliTHn *triggerHist,
                                                                       AliTHn *associateHist,
                                                                       Bool_t twoTrackEfficiencyCut,
                                                                       Float_t twoTrackEfficiencyCutValue,
                                                                       Float_t fTwoTrackCutMinRadius,
                                                                       Float_t bSign,
                                                                       Int_t step
                                                                       )
{
  if(!triggerHist || !associateHist)  return;
  
    Double_t  binscontTrig[2];
    Double_t  binscont[5];
    Double_t pidqa[4];
    for(Int_t i=0;i<triggerArray->GetEntriesFast();i++)  {
      AliAssociatedTrackYS * trigger = (AliAssociatedTrackYS*) triggerArray->At(i);
      if(!trigger)continue;
      Double_t triggerPt  = trigger->Pt();
      Double_t triggerEta  = trigger->Eta();
      Double_t triggerPhi  = trigger->Phi();
      Int_t trigFirstID=trigger->GetIDFirstDaughter();
      Int_t trigSecondID=trigger->GetIDSecondDaughter();
      Int_t trigID=trigger->GetID();
      binscontTrig[0]=triggerPt;
      binscontTrig[1]=centrality;
      
      pidqa[0]=trigger->Pt();
      pidqa[1]=trigger->Eta();
      pidqa[2]=trigger->Phi();
      pidqa[3]=lCentrality;
      if(fQA) fHistPIDQA->Fill(pidqa,0);
      
      triggerHist->Fill(binscontTrig,step); 
      for (Int_t j=0; j<selectedTrackArray->GetEntriesFast(); j++){
	AliAssociatedTrackYS* associate = (AliAssociatedTrackYS*) selectedTrackArray->At(j);
	if(!associate)continue;
	Int_t AssoFirstID=associate->GetIDFirstDaughter();
	Int_t AssoSecondID=associate->GetIDSecondDaughter();
	//	if(AssoFirstID==trigFirstID || AssoFirstID==trigSecondID ||AssoSecondID==trigFirstID || AssoSecondID==trigSecondID) continue;
	if(trigID==AssoFirstID || trigID==AssoSecondID) continue;
	binscont[0]=triggerEta-associate->Eta();
	binscont[1]=associate->Pt();
	binscont[2]=triggerPt;
	binscont[3]=centrality;
	binscont[4]=RangePhi(triggerPhi-associate->Phi());
	Int_t SpAsso=associate->WhichCandidate();
	if(fasso=="V0"||fasso=="Phi"||(fasso=="hadron"&&fPID)){
	  if(SpAsso<0) continue;
	  associateHist->Fill(binscont,SpAsso);
	}else if(fasso=="hadron"&&!fPID){
	  associateHist->Fill(binscont,0);
	}
      }
    }
  
}

void AliAnalysisTaskSEpPbCorrelationsYS::FillCorrelationTracksMixing(Double_t centrality,
								    Double_t pvxMix,
								    Double_t poolmax,
								    Double_t poolmin,
								    TObjArray *triggerArray,
								    TObjArray *selectedTrackArray,
								    AliTHn *triggerHist,
								    AliTHn *associateHist,
								    Bool_t twoTrackEfficiencyCut,
								    Float_t twoTrackEfficiencyCutValue,
								    Float_t fTwoTrackCutMinRadius,
								    Float_t bSign,
								     Int_t step
								     )
{

  if(!triggerHist || !associateHist)  return;
  
  Double_t  counterMix=0;
  AliEventPool* pool = fPoolMgr->GetEventPool(centrality, pvxMix);
  if (!pool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centrality, pvxMix));
  if (pool->IsReady() || pool->NTracksInPool() > fPoolMinNTracks  || pool->GetCurrentNEvents() > fMinEventsToMix)
    {
      Int_t nMix = pool->GetCurrentNEvents();
      for (Int_t jMix=0; jMix<nMix; jMix++){  
	TObjArray* mixEvents = pool->GetEvent(jMix);
	
	Double_t  binscontTrig[2];
	  Double_t  binscont[5];
	  for(Int_t i=0;i<triggerArray->GetEntriesFast();i++){
	    AliAssociatedTrackYS*trig=(AliAssociatedTrackYS*)triggerArray->At(i);
	    if(!trig) continue;
	    Double_t triggerPhi=trig->Phi();
	    Double_t triggerEta=trig->Eta();
	    Double_t triggerPt=trig->Pt();
	    counterMix++;
	    binscontTrig[0]=triggerPt;
	    binscontTrig[1]=centrality;
	    triggerHist->Fill(binscontTrig,step);
	    for(Int_t j=0;j<mixEvents->GetEntriesFast();j++){
	      AliAssociatedTrackYS*associate=(AliAssociatedTrackYS*)mixEvents->At(j);
	      if(!associate) continue;
	      binscont[0]=triggerEta-associate->Eta();
	      binscont[1]=associate->Pt();
	      binscont[2]=triggerPt;
	      binscont[3]=centrality;
	      binscont[4]=RangePhi(triggerPhi-associate->Phi());
	      Int_t SpAsso=associate->WhichCandidate();
	      if(fasso=="V0"||fasso=="Phi"||(fasso=="hadron"&&fPID)){
		if(SpAsso<0) continue;
		associateHist->Fill(binscont,SpAsso,1./(Double_t)nMix);
	      }else if(fasso=="hadron"&&!fPID){
		associateHist->Fill(binscont,0,1./(Double_t)nMix);
	      }
		
	    }
	  }
      }
    }
  
  TObjArray* tracksClone = new TObjArray;
  tracksClone->SetOwner(kTRUE);
  for (Int_t i=0; i<selectedTrackArray->GetEntriesFast(); i++)
    {
      AliAssociatedTrackYS* particle=(AliAssociatedTrackYS*) selectedTrackArray->At(i);
      tracksClone->Add(new AliAssociatedTrackYS(particle->Charge(),
					   particle->Eta(),
					   particle->Phi(),
					   particle->Pt(),
					   particle->GetID(),
					   particle->GetIDFirstDaughter(),
					   particle->GetIDSecondDaughter(),
					   particle->WhichCandidate(),
					   particle->Multiplicity()));
      
    }
  
  pool->UpdatePool(tracksClone);
  
}


Double_t AliAnalysisTaskSEpPbCorrelationsYS::RangePhi(Double_t DPhi)
{
  if (DPhi < -TMath::Pi()/2)  DPhi += 2*TMath::Pi();
  if (DPhi > 3*TMath::Pi()/2) DPhi -= 2*TMath::Pi();
  return DPhi;
}


Double_t AliAnalysisTaskSEpPbCorrelationsYS::RangePhi2(Double_t DPhi)
{
  DPhi=TMath::ATan2(TMath::Sin(DPhi),TMath::Cos(DPhi));
  if(DPhi<-1.178097) DPhi+=2*TMath::Pi();
  return DPhi;
}
