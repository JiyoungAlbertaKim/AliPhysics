/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnalysisTaskEmcalJetHadCorrAA.cxx
Author : Jiyoung Kim, jiyoung.kim@cern.ch */

#include "AliAnalysisTaskEmcalJetHadCorrAA.h"

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <THnSparse.h>
#include <THashList.h>

#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliLog.h>

#include "AliTLorentzVector.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAODHeader.h"
#include "AliAODEvent.h"
#include "AliAODVertex.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliMultSelection.h"
#include "AliPicoTrack.h"

/// PID includes
#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliEventCuts.h"

using namespace std;

/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalJetHadCorrAA);
/// \endcond

/**
 * Default constructor. Needed by ROOT I/O
 */
AliAnalysisTaskEmcalJetHadCorrAA::AliAnalysisTaskEmcalJetHadCorrAA() :
  AliAnalysisTaskEmcalJet("EmcalJetHadCorrAATask", kTRUE),
  fHistManager(),
  fPoolMgr(0x0),
  MxTracksArray(0x0),
  JetsArray(0x0),
  fNzvtxBins(11),
  fPIDResponse(0x0),
  fMultSelection(0x0),
  fHistCentrality_jy(0),
  fHistQACentrality_jy(0),
  fHistMixNumOfEventsInPool(0),
  fHistMixNumOfTracksInPool(0),
  doMakeQAplot(kFALSE),
  fEventCut()
{
  SetMakeGeneralHistograms(kTRUE);
  SetUseAliAnaUtils(kTRUE, kTRUE); //vertex selection for pileup rejection
  //SetUseAliAnaUtils(Bool_t b, Bool_t bRejPilup = kTRUE) { fUseAliAnaUtils = b ; fRejectPileup = bRejPilup  ; }

}

/**
 * Standard constructor. Should be used by the user.
 *
 * @param[in] name Name of the task
 */
AliAnalysisTaskEmcalJetHadCorrAA::AliAnalysisTaskEmcalJetHadCorrAA(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fPoolMgr(0x0),
  MxTracksArray(0x0),
  JetsArray(0x0),
  fHistManager(name),
  fNzvtxBins(11),
  fPIDResponse(0x0),
  fMultSelection(0x0),
  fHistCentrality_jy(0),
  fHistQACentrality_jy(0), 
  fHistMixNumOfEventsInPool(0),
  fHistMixNumOfTracksInPool(0),
  doMakeQAplot(kFALSE),
  fEventCut()
{
  SetMakeGeneralHistograms(kTRUE);
  SetUseAliAnaUtils(kTRUE, kTRUE);
}

/**
 * Destructor
 */
AliAnalysisTaskEmcalJetHadCorrAA::~AliAnalysisTaskEmcalJetHadCorrAA()
{
}

/**
 * Performing run-independent initialization.
 * Here the histograms should be instantiated.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::UserCreateOutputObjects()
{

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  AllocateTrackHistograms();
  AllocateJetHistograms();
  AllocateMixHistograms();

  TIter next(fHistManager.GetListOfHistograms());
  TObject* obj = 0;
  while ((obj = next())) {
    fOutput->Add(obj);
  }
 
  // pileup rejection qa figures
  /*fList = new TList();
  fList->SetOwner(true);  
  fEventCuts.AddQAplotsToList(fList);
  fOutput->Add(fList);*/    

  fEventCuts.AddQAplotsToList(fOutput);

  fHistMixNumOfEventsInPool = new TH3F("fHistMixNumOfEventsInPool", "histMixNumOfEventsInPool;Centrality [%];Z vertex [cm];N_{events}", 100, 0, 100, 20, -10, 10, 1500, 0, 1500); 
  fOutput->Add(fHistMixNumOfEventsInPool);

  fHistMixNumOfTracksInPool = new TH3F("fHistMixNumOfTracksInPool","histMixNumOfTracksInPool;Centrality [%];Z vertex [cm];N_{tracks}", 100, 0, 100, 20, -10, 10, 700, 0, 70000); 
  fOutput->Add(fHistMixNumOfTracksInPool);

  // Call for AliEventPoolManager (setup for mixing)  
  SetupForMixing();

  PostData(1, fOutput); // Post data for ALL output slots > 0 here.
}


/**
 * Setup for mixing events
 *
 */
void AliAnalysisTaskEmcalJetHadCorrAA::SetupForMixing()
{
   
  const Int_t iPoolsize = 1000;
  const Int_t iTrackDepth = 50000;
  const Int_t nCentBins = 5+1;
  const Int_t nZvtxBins = 10+1;
//  const Int_t nCentBins = 1;
//  const Int_t nZvtxBins = 1;


//  Double_t CentBinsArray[] = {0, 99};
//  Double_t ZvtxBinsArray[] = {-10, 10};
  Double_t CentBinsArray[] = {0, 10, 30, 50, 90, 100};
  Double_t ZvtxBinsArray[] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10};

  fPoolMgr = new AliEventPoolManager(iPoolsize, iTrackDepth, nCentBins, CentBinsArray, nZvtxBins, ZvtxBinsArray);
  fPoolMgr->SetTargetValues(iTrackDepth, 0.1, 5); // (trackDepth, fraction, event) for IsReady()
  //fPoolMgr->SetDebug(0);

}

/*
 * This function is for event-mixing 
 * Event taken from TrackLoop gives Cent, Zvtx info.
 * Based on those, EventPool can be called
 * To Fill mixedHistograms with events in the pool
 */
void AliAnalysisTaskEmcalJetHadCorrAA::DoMixing(double cent, double zvertex, TObjArray* trackArray)
{
  AliEventPool* pool = fPoolMgr->GetEventPool(cent, zvertex);
  if (!pool) {
    AliInfo(Form("No pool found for centrality = %f, zVtx = %f", cent, zvertex));
  } else {

      if (pool->IsReady()) { // check whether eventpool is ready or not
 
        //pool->PrintInfo();
        Int_t nMix = pool->GetCurrentNEvents();
        Int_t nTrk = pool->NTracksInPool();
        //AliInfo(Form("Get current number of events in the pool: %i", nMix));

        fHistMixNumOfEventsInPool->Fill(cent, zvertex, nMix); 
        fHistMixNumOfTracksInPool->Fill(cent, zvertex, nTrk); 

        // Get cBin and zBin to pass those to FillMixedHistograms
        Int_t nCentBin = 5+1;
        float cBinArray[] = {0, 10, 30, 50, 90, 100}; // 6 bins for LHC15o

        // Get cBin from percentage
        int cBin = -1;
        for(int icbin = 0; icbin<nCentBin; icbin++) {
          if(0 > cent) continue;
          float cBinMin = cBinArray[icbin];
          float cBinMax = cBinArray[icbin+1];
          if((cBinMin <= cent) && (cent < cBinMax)) cBin = icbin;
        }

        // Fill zVertex and zBin
        Int_t nZvtxBin = 10+1;
        float zBinArray[] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}; // 11 bins

        int zBin = -1;
        for(int izbin = 0; izbin<nZvtxBin; izbin++) {
          if(zvertex == -999.) continue; // exclude what have initial value
          float zBinMin = zBinArray[izbin];
          float zBinMax = zBinArray[izbin+1];
          if((zBinMin <= zvertex) && (zvertex < zBinMax)) zBin = izbin; // to fill from the 1st bin
        }
      
        for (Int_t jMix=0; jMix<nMix; jMix++) {
          MxTracksArray = pool->GetEvent(jMix);
          if (!MxTracksArray) continue;
          ///FillMixedHistos( trigJets, fTracksMixing, cBin, zBin, EPangle, bSignDoMixing, 1./nMix ); 
          FillMixedHistograms(trackArray, MxTracksArray, cBin, zBin, nMix);
        }
      } 
      pool->UpdatePool(trackArray);
      //Printf("event pool update : %i", pool->WasUpdated());
  }
}  

/*
 * This function allocates the histograms for basic tracking QA.
 * A set of histograms (pT, eta, phi, difference between kinematic properties
 * at the vertex and at the EMCal surface, number of tracks) is allocated
 * per each particle container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::AllocateTrackHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray); // loop 1
  while ((partCont = static_cast<AliParticleContainer*>(next()))) { // loop 2 // these loop are needed in calling the name
    groupname = partCont->GetName();
    fHistManager.CreateHistoGroup(groupname);
    Printf("%s", partCont->GetName());

    for (Int_t cent = 0; cent < fNcentBins; cent++) {

      // To check z vertex dependency in delta eta distribution
      for (Int_t zv = 0; zv < fNzvtxBins; zv++) {

        //histname = TString::Format("%s/histTrackPhiEtaPt_Zvtx_%d_%d", groupname.Data(), cent, zv);
        //histtitle = TString::Format("%s;#it{#phi}_{track};counts", histname.Data());
        //fHistManager.CreateTH3(histname, histtitle, 64, 0, TMath::TwoPi(), 40, -1, 1, fNbins, fMinBinPt, fMaxBinPt / 2);

        histname = TString::Format("%s/histTrackPhiEta_Zvtx_%d_%d", groupname.Data(), cent, zv);
        histtitle = TString::Format("%s;#it{#phi}_{track};#it{#eta}_{track};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 64, 0, TMath::TwoPi(), 40, -1, 1);
 
      } 

      //Plots for Vz vs. eta
      histname = TString::Format("%s/histTrackZvEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};Z vertex", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 40, -1, 1, 20, -10, 10);
     
      // combination plot for projections
      histname = TString::Format("%s/histTrackPhiEtaPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi}_{track};#it{#eta}_{track};#it{p}_{T,track} (GeV/#it{c})", histname.Data());
      fHistManager.CreateTH3(histname, histtitle, 64, 0, TMath::TwoPi(), 40, -1, 1, fNbins, fMinBinPt, fMaxBinPt /2); 
    
      // Particle Identification (PID) plots 
      histname = TString::Format("%s/histPID_TPCsignal_%d", groupname.Data(), cent);
      //histtitle = TString::Format("%s;#it{p}_{T,track} (GeV/#it{c});TPC signal (arb. units)", histname.Data()); 
      histtitle = TString::Format("%s;p (GeV);TPC signal (arb. units)", histname.Data()); 
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 500, 0, 1000);

      histname = TString::Format("%s/histPID_ITSsignal_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;p (GeV);ITS signal (arb. units)", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 300, 0, 300);

      histname = TString::Format("%s/histPID_TOFsignal_%d", groupname.Data(), cent);       
      histtitle = TString::Format("%s;p (GeV);TOF signal [ns]", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 300, 0, 30);

      histname = TString::Format("%s/histPID_TPCNSigmaProton_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;p (GeV);n#sigma", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 200, -10, 10);

      histname = TString::Format("%s/histPID_TOFNSigmaProton_%d", groupname.Data(), cent); 
      histtitle = TString::Format("%s;p (GeV);n#sigma", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 200, -10, 10);

      histname = TString::Format("%s/histPID_TOFNSigma_TPCNsigma_Proton_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;n#sigma_{TPC};n#sigma_{TOF};p (GeV)", histname.Data());
      fHistManager.CreateTH3(histname, histtitle, 200, -10, 10, 200, -10, 10, 200, 0.1, 30);  

      /// nsigma after TPC or TOF cut    
      histname = TString::Format("%s/histPID_TPC_3sigmaTOF_Proton_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;p (GeV);n#sigma", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 200, -10, 10);
 
      histname = TString::Format("%s/histPID_TPC_2sigmaTOF_Proton_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;p (GeV);n#sigma", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 200, -10, 10);
 
      histname = TString::Format("%s/histPID_TOF_3sigmaTPC_Proton_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;p (GeV);n#sigma", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 200, -10, 10);
 
      histname = TString::Format("%s/histPID_TOF_2sigmaTPC_Proton_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;p (GeV);n#sigma", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 200, 0.1, 30, 200, -10, 10);
 

      histname = TString::Format("%s/histNTracks_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;number of tracks;events", histname.Data());
      if (fForceBeamType != kpp) {
        fHistManager.CreateTH1(histname, histtitle, 5000, 0, 5000);
      }
      else {
        fHistManager.CreateTH1(histname, histtitle, 200, 0, 200);
      }
    }//centrality loop
  }//container loop
}

/*
 * This function allocates the histograms for basic jet QA.
 * A set of histograms (pT, eta, phi, area, number of jets, corrected pT) is allocated
 * per each jet container and per each centrality bin.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::AllocateJetHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;
  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    fHistManager.CreateHistoGroup(groupname);
    Printf("%s", jetCont->GetName());

    for (Int_t cent = 0; cent < fNcentBins; cent++) {

      // To check z vertex dependency in delta eta distribution
      for (Int_t ztx = 0; ztx < fNzvtxBins; ztx++) {
        /// for jets
        ///histname = TString::Format("%s/histJetPhiEtaPt_%d_%d", groupname.Data(), cent, ztx);
        ///histtitle = TString::Format("%s;#it{#phi};#it{#eta};#it{p}_{T,jet} (GeV/#it{c})", histname.Data());
        ///fHistManager.CreateTH3(histname, histtitle, 64, 0, TMath::TwoPi(), 40, -1, 1, fNbins, fMinBinPt, fMaxBinPt / 2);
        histname = TString::Format("%s/histJetPhiEta_Zvtx_%d_%d", groupname.Data(), cent, ztx);
        histtitle = TString::Format("%s;#it{#phi};#it{#eta};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 64, 0, TMath::TwoPi(), 40, -1, 1);
        
        /// for correlations
        ///histname = TString::Format("%s/histCorrPhiEtaPt_%d_%d", groupname.Data(), cent, ztx);
        ///histtitle = TString::Format("%s;#it{#Delta#phi};#it{#Delta#eta};#it{p}_{T,track} (GeV/#it{c})", histname.Data());
        ///fHistManager.CreateTH3(histname, histtitle, 64, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 64, -1.6, 1.6, fNbins, fMinBinPt, fMaxBinPt / 2);
	histname = TString::Format("%s/histCorrPhiEta_Zvtx_%d_%d", groupname.Data(), cent, ztx);
        histtitle = TString::Format("%s;#it{#Delta#phi};#it{#Delta#eta};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 64, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 64, -1.6, 1.6);
      }

      //Plots for Vz vs. eta
      histname = TString::Format("%s/histJetZvEta_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#eta}_{track};Z vertex", histname.Data());
      fHistManager.CreateTH2(histname, histtitle, 40, -1, 1, 20, -10, 10);
 
      /// for jets
      histname = TString::Format("%s/histJetPhiEtaPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#phi};#it{#eta};#it{p}_{T,jet} (GeV/#it{c})", histname.Data());
      fHistManager.CreateTH3(histname, histtitle, 64, 0, TMath::TwoPi(), 40, -1, 1, fNbins, fMinBinPt, fMaxBinPt / 2);
 
      // Correlation plots
      histname = TString::Format("%s/histCorrPhiEtaPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#Delta#phi};#it{#Delta#eta};#it{p}_{T,track} (GeV/#it{c})", histname.Data());
      fHistManager.CreateTH3(histname, histtitle, 64, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 64, -1.6, 1.6, fNbins, fMinBinPt, fMaxBinPt / 2); 
      
      histname = TString::Format("%s/histCorrJetPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
    }
  }
}

/*
 * This function allocates the histograms for mixing plots.
 * need to make explanations
 * 
 */
void AliAnalysisTaskEmcalJetHadCorrAA::AllocateMixHistograms()
{
  TString histname;
  TString histtitle;
  TString groupname;

  AliJetContainer* jetCont = 0;
  TIter next(&fJetCollArray);
  while ((jetCont = static_cast<AliJetContainer*>(next()))) {
    groupname = jetCont->GetName();
    groupname += "_";
    groupname += "Mix";
    fHistManager.CreateHistoGroup(groupname);

    /*histname = TString::Format("%s/histMixNumOfEventsInPool", groupname.Data());
    histtitle = TString::Format("%s;Centrality [%];Z vertex [cm];N_{events}", histname.Data());
    fHistManager.CreateTH3(histname, histtitle, 100, 0, 100, 20, -10, 10, 1500, 0, 1500); 

    histname = TString::Format("%s/histMixNumOfTracksInPool", groupname.Data());
    histtitle = TString::Format("%s;Centrality [%];Z vertex [cm];N_{tracks}", histname.Data());
    fHistManager.CreateTH3(histname, histtitle, 100, 0, 100, 20, -10, 10, 7000, 7000, 14000); */
          
    for (Int_t cent = 0; cent < fNcentBins; cent++) {
      // To check z vertex dependency in delta eta distribution
      for (Int_t zvtx = 0; zvtx < fNzvtxBins; zvtx++) {
        ///histname = TString::Format("%s/histMixCorrPhiEtaPt_%d_%d", groupname.Data(), cent, zvtx);
        ///histtitle = TString::Format("%s;#it{#Delta#phi};#it{#Delta#eta};#it{p}_{T,track} (GeV/#it{c})", histname.Data());
        ///fHistManager.CreateTH3(histname, histtitle, 64, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 64, -1.6, 1.6, fNbins, fMinBinPt, fMaxBinPt / 2);
        histname = TString::Format("%s/histMixCorrPhiEta_Zvtx_%d_%d", groupname.Data(), cent, zvtx);
        histtitle = TString::Format("%s;#it{#Delta#phi};#it{#Delta#eta};counts", histname.Data());
        fHistManager.CreateTH2(histname, histtitle, 64, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 64, -1.6, 1.6);
      }
      histname = TString::Format("%s/histMixCorrPhiEtaPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{#Delta#phi};#it{#Delta#eta};#it{p}_{T,track} (GeV/#it{c})", histname.Data());
      fHistManager.CreateTH3(histname, histtitle, 64, -0.5*TMath::Pi(), 1.5*TMath::Pi(), 64, -1.6, 1.6, fNbins, fMinBinPt, fMaxBinPt / 2); 
    
      // For others
      histname = TString::Format("%s/histMixJetCorrPt_%d", groupname.Data(), cent);
      histtitle = TString::Format("%s;#it{p}_{T,jet}^{corr} (GeV/#it{c});counts", histname.Data());
      fHistManager.CreateTH1(histname, histtitle, fNbins, -fMaxBinPt / 2, fMaxBinPt / 2);
    }
  }
}

/**
 * The body of this function should contain instructions to fill the output histograms.
 * This function is called inside the event loop, after the function Run() has been
 * executed successfully (i.e. it returned kTRUE).
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetHadCorrAA::FillHistograms()
{
  DoJetLoop();
  DoTrackLoop();

  return kTRUE;
}

/**
 * This function is called inside the Run() funtions
 * To fill MixingHistograms, call DoMixLoop() inside this function
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetHadCorrAA::FillMixedHistograms(TObjArray* TrkArray, TObjArray* MxTrkArray, int CentBin, int ZvtxBin, int nMix)
{ 
  DoMixLoop(TrkArray, MxTrkArray, CentBin, ZvtxBin, nMix);
  return kTRUE;
}

/**
 * This function performs a loop over the reconstructed jets & mixed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::DoMixLoop(TObjArray* part, TObjArray* partMix, int cBin, int zBin, int nMix)
{

  int nMxtrks = partMix->GetEntriesFast();
  int ntrks = part->GetEntriesFast();

  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  Double_t kJetCorrPtCut = 0; 
  Float_t JetRadius;
  Float_t jetCorrPt = 0;

  Double_t zvertex = -999;
  zvertex = InputEvent()->GetPrimaryVertex()->GetZ();

  TIter next(&fJetCollArray); // loop 1 : definition for next sentence While ~
  while ((jetCont = static_cast<AliJetContainer*>(next()))) { // loop 2 : call two jet containers R020, R040. if you need you can add more containers
    groupname = jetCont->GetName();
    groupname += "_";
    groupname += "Mix";
 
    UInt_t count = 0;

    JetRadius = jetCont->GetJetRadius();    
    if (TMath::Abs(JetRadius - 0.2) < 0.01 ) {
      kJetCorrPtCut = 30.0;
    } else if (TMath::Abs(JetRadius - 0.4) < 0.01) {
      kJetCorrPtCut = 40.0;
    }

    /*histname = TString::Format("%s/histMixNumOfEventsInPool", groupname.Data());
    fHistManager.FillTH3(histname, fCent, zvertex, nMix); 
 
    histname = TString::Format("%s/histMixNumOfTracksInPool", groupname.Data());
    fHistManager.FillTH3(histname, fCent, zvertex, nMxtrks);*/ 
 
    for(auto jet : jetCont->accepted()) { // loop 3 : loop for jet by jet
      if (!jet) continue;
      count++;

      if (jetCont->GetRhoParameter()) { 
        jetCorrPt = jet->Pt() - jetCont->GetRhoVal() * jet->Area();

        if(jetCorrPt >= kJetCorrPtCut) {

          histname = TString::Format("%s/histMixJetCorrPt_%d", groupname.Data(), cBin);
          fHistManager.FillTH1(histname, jetCorrPt);

          for(int i = 0; i < ntrks; i++) { // trackloop with same events 
            AliPicoTrack *particle = static_cast<AliPicoTrack*>(part->At(i));
            if(!(particle)) return;

            // Calculation of delta phi between a jet and a track
            Float_t deltaPhiSame = 0;
            Float_t deltaEtaSame = 0;
            deltaPhiSame = jet->Phi() - particle->Phi();
            deltaEtaSame = jet->Eta() - particle->Eta();
 
            // For axis adjustment, from -pi/2 to +3*pi/2
            if (deltaPhiSame > TMath::Pi()*3.0/2.0) {
              deltaPhiSame -= TMath::TwoPi();
            }else if (deltaPhiSame < -TMath::Pi()/2.0) {
              deltaPhiSame += TMath::TwoPi();
            }
          } 
 
          for(int i = 0; i < nMxtrks; i++) { // trackloop with mixed events
            AliPicoTrack *Mixparticle = static_cast<AliPicoTrack*>(partMix->At(i));
            if(!(Mixparticle)) return;

            // Calculation of delta phi between a jet and a track
            Float_t deltaPhi = 0;
            Float_t deltaEta = 0;
            deltaPhi = jet->Phi() - Mixparticle->Phi();
            deltaEta = jet->Eta() - Mixparticle->Eta();
 
            // For axis adjustment, from -pi/2 to +3*pi/2
            if (deltaPhi > TMath::Pi()*3.0/2.0) {
              deltaPhi -= TMath::TwoPi();
            }else if (deltaPhi < -TMath::Pi()/2.0) {
              deltaPhi += TMath::TwoPi();
            }

            // filling histogram - delta phi, eta ..  
            // To check z vertex dependency in delta eta distribution
            ///histname = TString::Format("%s/histMixCorrPhiEtaPt_%d_%d", groupname.Data(), cBin, zBin);
            ///fHistManager.FillTH3(histname, deltaPhi, deltaEta, Mixparticle->Pt()); 
            histname = TString::Format("%s/histMixCorrPhiEta_Zvtx_%d_%d", groupname.Data(), cBin, zBin);
            fHistManager.FillTH2(histname, deltaPhi, deltaEta); 
 
            // Correlation plots
            histname = TString::Format("%s/histMixCorrPhiEtaPt_%d", groupname.Data(), cBin);
            fHistManager.FillTH3(histname, deltaPhi, deltaEta, Mixparticle->Pt());
          }  
        }// corrected JetPt cut 
      }// Check rho values
    }// jet loop
  }// jet Container 
}


/**
 * This function performs a loop over the reconstructed jets
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::DoJetLoop()
{
  TString histname;
  TString groupname;
  AliJetContainer* jetCont = 0;
  Double_t kJetCorrPtCut = 0; 
  Float_t JetRadius;
  Float_t jetCorrPt = 0;

  TIter next(&fJetCollArray); // loop 1 : definition for next sentence While ~
  while ((jetCont = static_cast<AliJetContainer*>(next()))) { // loop 2 : call two jet containers R020, R040. if you need you can add more containers
    groupname = jetCont->GetName();
    UInt_t count = 0;

    Double_t zvertex = -999;
    zvertex = InputEvent()->GetPrimaryVertex()->GetZ();

    // Fill zVertex and zBin
    Int_t nZvtxBin = 10+1;
    float zBinArray[] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}; // 11 bins

    int zBin = -1;
    for(int izbin = 0; izbin<nZvtxBin; izbin++) {
      if(zvertex == -999.) continue; // exclude what have initial value
      float zBinMin = zBinArray[izbin];
      float zBinMax = zBinArray[izbin+1];
      if((zBinMin <= zvertex) && (zvertex < zBinMax)) zBin = izbin; // to fill from the 1st bin
    }
      
    JetRadius = jetCont->GetJetRadius();    
    if (TMath::Abs(JetRadius - 0.2) < 0.01 ) {
      //Printf("Wow, 02, %f", JetRadius);
      kJetCorrPtCut = 30.0;
      //Printf("Applied corr pt cut: %f", kJetCorrPtCut);
    } else if (TMath::Abs(JetRadius - 0.4) < 0.01) {
      //Printf("Wow, 04, %f", JetRadius);
      kJetCorrPtCut = 40.0;
      //Printf("Applied corr pt cut: %f", kJetCorrPtCut);
    }

    for(auto jet : jetCont->accepted()) { // loop 3 : loop for jet by jet
      if (!jet) continue;
      count++;
      if (jetCont->GetRhoParameter()) { 

        jetCorrPt = jet->Pt() - jetCont->GetRhoVal() * jet->Area();

        if(jetCorrPt >= kJetCorrPtCut) {
      
          // combination plot for projections
          histname = TString::Format("%s/histJetPhiEtaPt_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH3(histname, jet->Phi(), jet->Eta(), jetCorrPt);

          // To check z vertex dependency in delta eta distribution
          ///histname = TString::Format("%s/histJetPhiEtaPt_%d_%d", groupname.Data(), fCentBin, zBin);
          ///fHistManager.FillTH3(histname, jet->Phi(), jet->Eta(), jetCorrPt);
 	        histname = TString::Format("%s/histJetPhiEta_Zvtx_%d_%d", groupname.Data(), fCentBin, zBin);
          fHistManager.FillTH2(histname, jet->Phi(), jet->Eta());
 
          // Vz vs. eta
          histname = TString::Format("%s/histJetZvEta_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH2(histname, jet->Eta(), zvertex);
 
          // Particle (track) loop
          AliParticleContainer* partCont = 0;
          TIter next(&fParticleCollArray);
          while ((partCont = static_cast<AliParticleContainer*>(next()))) {
            UInt_t trkcount = 0;
            for(auto part : partCont->accepted()) {
              if (!part) continue;
              trkcount++; // number of tracks in a event (in jet loop)
              
              // Calculation of delta phi between a jet and a track
              Float_t deltaPhi = 0;
              Float_t deltaEta = 0;
              deltaPhi = jet->Phi() - part->Phi();
              deltaEta = jet->Eta() - part->Eta();
 
              // For axis adjustment, from -pi/2 to +3*pi/2
              if (deltaPhi > TMath::Pi()*3.0/2.0) {
                deltaPhi -= TMath::TwoPi();
              }else if (deltaPhi < -TMath::Pi()/2.0) {
                deltaPhi += TMath::TwoPi();
              } 

              // filling histogram - delta phi, eta ..  
              histname = TString::Format("%s/histCorrPhiEtaPt_%d", groupname.Data(), fCentBin);
              fHistManager.FillTH3(histname, deltaPhi, deltaEta, part->Pt()); 

              // To check z vertex dependency in delta eta distribution
              ///histname = TString::Format("%s/histCorrPhiEtaPt_%d_%d", groupname.Data(), fCentBin, zBin);
              ///fHistManager.FillTH3(histname, deltaPhi, deltaEta, part->Pt());
              histname = TString::Format("%s/histCorrPhiEta_Zvtx_%d_%d", groupname.Data(), fCentBin, zBin);
              fHistManager.FillTH2(histname, deltaPhi, deltaEta); 

            } // track loop
          } // track Container
          
          // To check correlation entries: entries =  total number of jets
          histname = TString::Format("%s/histCorrJetPt_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH1(histname, jetCorrPt); 

        } // jet pt
      } // check rho value      
    } // jet loop
  }// jet Container 
}

/**
 * This function performs a loop over the reconstructed tracks
 * in the current event and fills the relevant histograms.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::DoTrackLoop()
{
  AliClusterContainer* clusCont = GetClusterContainer(0);

  TObjArray *trkArray = new TObjArray; //Track array for event-mixing 

  Double_t zvertex = -999;
  zvertex = InputEvent()->GetPrimaryVertex()->GetZ();

  // Fill zVertex and zBin
  Int_t nZvtxBin = 10+1;
  float zBinArray[] = {-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10}; // 11 bins

  int zBin = -1;
  for(int izbin = 0; izbin<nZvtxBin; izbin++) {
    if(zvertex == -999.) continue; // exclude what have initial value
    float zBinMin = zBinArray[izbin];
    float zBinMax = zBinArray[izbin+1];
    if((zBinMin <= zvertex) && (zvertex < zBinMax)) zBin = izbin; // to fill from the 1st bin
  }
  
  TString histname;
  TString groupname;
  AliParticleContainer* partCont = 0;
  TIter next(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    groupname = partCont->GetName();
    UInt_t count = 0;

    for(auto part : partCont->accepted()) {
      if (!part) continue;
      count++;

      // Add tracks in the array for event-mixing 
      trkArray->Add(new AliPicoTrack(part->Pt(), part->Eta(), part->Phi(), part->Charge(), 0, 0, 0, 0));

      // combination plot for projections
      histname = TString::Format("%s/histTrackPhiEtaPt_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH3(histname, part->Phi(), part->Eta(), part->Pt()); 

      histname = TString::Format("%s/histTrackZvEta_%d", groupname.Data(), fCentBin);
      fHistManager.FillTH2(histname, part->Eta(), zvertex);
 
      // To check z vertex dependency in delta eta distribution
      ///histname = TString::Format("%s/histTrackPhiEtaPt_Zvtx_%d_%d", groupname.Data(), fCentBin, zBin);
      ///fHistManager.FillTH3(histname, part->Phi(), part->Eta(), part->Pt());
      histname = TString::Format("%s/histTrackPhiEta_Zvtx_%d_%d", groupname.Data(), fCentBin, zBin);
      fHistManager.FillTH2(histname, part->Phi(), part->Eta());


      // some variables for PID
      Double_t pt = -999, dEdx = -999, ITSsig = -999, TOFsig = -999, charge = -999;
      Double_t TPCmom = -999, P = -999;

      // nSigma of particles in TPC, TOF, and ITS
      Double_t nSigmaPion_TPC, nSigmaProton_TPC, nSigmaKaon_TPC;
      Double_t nSigmaPion_TOF, nSigmaProton_TOF, nSigmaKaon_TOF;
      Double_t nSigmaPion_ITS, nSigmaProton_ITS, nSigmaKaon_ITS;

      if(fPIDResponse) {
        // get parameters of track
        charge = part->Charge();    // charge of track
        pt     = part->Pt();        // pT of track
        P      = part->P();         // absolute momentum? = sqrt(pt^2+pz^2)

        if (partCont->GetLoadedClass()->InheritsFrom("AliVTrack")) {
          const AliVTrack* track = static_cast<const AliVTrack*>(part);
       
          // get momentums
          TPCmom = track->GetTPCmomentum(); //momemtum at the inner wall of tpc
           
          // get detector signals
          //dEdx = fPIDResponse->GetTPCResponse().GetTrackdEdx(track); //TPCsig
          dEdx = track->GetTPCsignal(); //TPCsig
          ITSsig = track->GetITSsignal();
          TOFsig = track->GetTOFsignal()/1000.;

        }

        // TPC nSigma's // maybe do not need to call nSigma here, just need to call once
        nSigmaPion_TPC = fPIDResponse->NumberOfSigmasTPC(part,AliPID::kPion);
        nSigmaKaon_TPC = fPIDResponse->NumberOfSigmasTPC(part,AliPID::kKaon);
        nSigmaProton_TPC = fPIDResponse->NumberOfSigmasTPC(part,AliPID::kProton);

        // TOF nSigma's
        nSigmaPion_TOF = fPIDResponse->NumberOfSigmasTOF(part,AliPID::kPion);
        nSigmaKaon_TOF = fPIDResponse->NumberOfSigmasTOF(part,AliPID::kKaon);
        nSigmaProton_TOF = fPIDResponse->NumberOfSigmasTOF(part,AliPID::kProton);

        // ITS nSigma's
        nSigmaPion_ITS = fPIDResponse->NumberOfSigmasITS(part,AliPID::kPion);
        nSigmaKaon_ITS = fPIDResponse->NumberOfSigmasITS(part,AliPID::kKaon);
        nSigmaProton_ITS = fPIDResponse->NumberOfSigmasITS(part,AliPID::kProton);


        // To fill histograms
        histname = TString::Format("%s/histPID_TPCsignal_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH2(histname, P, dEdx);
        histname = TString::Format("%s/histPID_ITSsignal_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH2(histname, P, ITSsig);
        histname = TString::Format("%s/histPID_TOFsignal_%d", groupname.Data(), fCentBin);
        ///fHistManager.FillTH2(histname, part->Pt(), TOFsig);
        fHistManager.FillTH2(histname, P, TOFsig);
       
        histname = TString::Format("%s/histPID_TPCNSigmaProton_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH2(histname, P, nSigmaProton_TPC);
        histname = TString::Format("%s/histPID_TOFNSigmaProton_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH2(histname, P, nSigmaProton_TOF);

        histname = TString::Format("%s/histPID_TOFNSigma_TPCNsigma_Proton_%d", groupname.Data(), fCentBin);
        fHistManager.FillTH3(histname, nSigmaProton_TPC, nSigmaProton_TOF, P);
        
        //TPC after TOF cut
        if (TMath::Abs(nSigmaProton_TOF)<3.) {
          histname = TString::Format("%s/histPID_TPC_3sigmaTOF_Proton_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH2(histname, P, nSigmaProton_TPC);
 
          if (TMath::Abs(nSigmaProton_TOF)<2.) {
            histname = TString::Format("%s/histPID_TPC_2sigmaTOF_Proton_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH2(histname, P, nSigmaProton_TPC);
          }
        }
        //TOF after TPC cut
        if (TMath::Abs(nSigmaProton_TPC)<3.) {
          histname = TString::Format("%s/histPID_TOF_3sigmaTPC_Proton_%d", groupname.Data(), fCentBin);
          fHistManager.FillTH2(histname, P, nSigmaProton_TOF);
 
          if (TMath::Abs(nSigmaProton_TPC)<2.) {
            histname = TString::Format("%s/histPID_TOF_2sigmaTPC_Proton_%d", groupname.Data(), fCentBin);
            fHistManager.FillTH2(histname, P, nSigmaProton_TOF);
          }  
        }

        
/*
///////// trial

 // Fill PID qa histograms for the TOF
  //   Here also the TPC histograms after TOF selection are filled
  //

  AliVEvent *event=InputEvent();

  Int_t ntracks=event->GetNumberOfTracks();
  for(Int_t itrack = 0; itrack < ntracks; itrack++){
    AliVTrack *track=(AliVTrack*)event->GetTrack(itrack);

    //
    //basic track cuts
    //
    ULong_t status=track->GetStatus();
    // not that nice. status bits not in virtual interface
    // TPC refit + ITS refit +
    // TOF out + TOFpid +
    // kTIME
    if (!((status & AliVTrack::kTPCrefit) == AliVTrack::kTPCrefit) ||
        !((status & AliVTrack::kITSrefit) == AliVTrack::kITSrefit) ||
//         !( (status & AliVTrack::kTPCpid  ) == AliVTrack::kTPCpid ) || //removes light nuclei, so it is out for the moment
        !((status & AliVTrack::kTOFout  ) == AliVTrack::kTOFout  ) ||
        //!((status & AliVTrack::kTOFpid  ) == AliVTrack::kTOFpid  ) || // not valid any longer with new TOF structure
        !((status & AliVTrack::kTIME    ) == AliVTrack::kTIME    ) ) continue;

    Float_t nCrossedRowsTPC = track->GetTPCClusterInfo(2,1);
    Float_t  ratioCrossedRowsOverFindableClustersTPC = 1.0;
    if (track->GetTPCNclsF()>0) {
      ratioCrossedRowsOverFindableClustersTPC = nCrossedRowsTPC/track->GetTPCNclsF();
    }

    if ( nCrossedRowsTPC<70 || ratioCrossedRowsOverFindableClustersTPC<.8 ) continue;


    Double_t mom=track->P();
    Double_t momTPC=track->GetTPCmomentum();

    for (Int_t ispecie=0; ispecie<AliPID::kSPECIESC; ++ispecie){
      //TOF nSigma
      Double_t nSigmaTOF=fPIDResponse->NumberOfSigmasTOF(track, (AliPID::EParticleType)ispecie);
      Double_t nSigmaTPC=fPIDResponse->NumberOfSigmasTPC(track, (AliPID::EParticleType)ispecie);

      //TPC after TOF cut
      TH2 *h=(TH2*)fListQAtpctof->At(ispecie);
      if (h && TMath::Abs(nSigmaTOF)<3.) h->Fill(momTPC,nSigmaTPC);

      //TOF after TPC cut
      h=(TH2*)fListQAtpctof->At(ispecie+AliPID::kSPECIESC);
      if (h && TMath::Abs(nSigmaTPC)<3.) h->Fill(mom,nSigmaTOF);

      //EMCAL after TOF and TPC cut
      h=(TH2*)fListQAtpctof->At(ispecie+2*AliPID::kSPECIESC);
      if (h && TMath::Abs(nSigmaTOF)<3. && TMath::Abs(nSigmaTPC)<3. ){

  Int_t nMatchClus = track->GetEMCALcluster();
  Double_t pt      = track->Pt();
  Double_t eop     = -1.;

  if(nMatchClus > -1){

    AliVCluster *matchedClus = (AliVCluster*)event->GetCaloCluster(nMatchClus);

    if(matchedClus){

      // matched cluster is EMCAL
      if(matchedClus->IsEMCAL()){

        Double_t fClsE       = matchedClus->E();
        eop                  = fClsE/mom;

        h->Fill(pt,eop);


      }
    }
  }
      }
    }
  }
////////////// here !!!
*/

      }// PID
    }// partCont->accepted track loop
    histname = TString::Format("%s/histNTracks_%d", groupname.Data(), fCentBin);
    fHistManager.FillTH1(histname, count);
  }// End of ParticleContainer loop 

  // For event-mixing 
  bool doMixing; doMixing = kTRUE;
  if(doMixing) DoMixing(fCent, zvertex, trkArray);
}

/**
 * This function is executed automatically for the first event.
 * Some extra initialization can be performed here.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::ExecOnce()
{
  AliAnalysisTaskEmcalJet::ExecOnce();
}

/**
 * Run analysis code here, if needed.
 * It will be executed before FillHistograms().
 * If this function return kFALSE, FillHistograms() will *not*
 * be executed for the current event
 * @return Always kTRUE
 */
Bool_t AliAnalysisTaskEmcalJetHadCorrAA::Run()
{
  
  AliVEvent *ev = InputEvent();
  if (!fEventCut.AcceptEvent(ev)) {
    return kFALSE;
  }

/*  AliMultSelection *fMultSelection = static_cast<AliMultSelection*>(ev->FindListObject("MultSelection"));
  if (fMultSelection) {
    //fCent = fMultSelection->GetMultiplicityPercentile(fCentEst.Data());
    fCent = fMultSelection->GetMultiplicityPercentile(fCentEst.Data(), kTRUE);
    Int_t qual = fMultSelection->GetEvSelCode(); // modified by JY
    //Printf("qual: %i", qual);
    fHistCentrality_jy->Fill(fCent);
    fHistQACentrality_jy->Fill(qual); // modified by JY
  } else {
    AliWarning(Form("%s: Could not retrieve centrality information! Assuming 99", GetName()));
  }*/

  // pileup rejection (in your UserExec) - JY
  UInt_t fSelectMask= fInputHandler->IsEventSelected();
  Bool_t isINT7selected = fSelectMask& AliVEvent::kINT7;

  // additional event cut for pileup rejection - JY
  Int_t multEsd = ((AliAODHeader*)ev->GetHeader())->GetNumberOfESDTracks();
  const Int_t nTracks = ev->GetNumberOfTracks();
  Int_t multTPC=0;
  for (Int_t it1 = 0; it1 < nTracks; it1++) {
    AliAODTrack* aodTrk1 = (AliAODTrack*)ev->GetTrack(it1);
    if (!aodTrk1) continue;
    if (aodTrk1->TestFilterBit(128)) multTPC++;
  }
  if(multEsd -3.38*multTPC<15000){
    // keep the event 
    return kTRUE;
  }else{
   // reject the event 
    return kFALSE;
  }

  //return kTRUE;
}

/**
 * This function is called once at the end of the analysis.
 */
void AliAnalysisTaskEmcalJetHadCorrAA::Terminate(Option_t *) 
{
}
