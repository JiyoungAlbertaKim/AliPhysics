#ifndef ALIANALYSISTASKEMCALJETHADCORRAA_H
#define ALIANALYSISTASKEMCALJETHADCORRAA_H
/**
 * \file AliAnalysisTaskEmcalJetHadCorrAA.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetHadCorrAA
 *
 * In this header file the class AliAnalysisTaskEmcalJetHadCorrAA is declared.
 * This is a sample task that shows how to write a simple user analysis task
 * using the EMCal jet framework. It is also used to do automatic benchmark
 * tests of the software.
 *
 * \author Jiyoung Kim <jiyoung.kim@cern.ch>, Heidelberg University
 * \date Nov 8, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

#include "AliEventPoolManager.h"
#include "TObjArray.h"
#include "AliPIDResponse.h"
#include "AliMultSelection.h"
#include "AliEventCuts.h"

/**
 * \class AliAnalysisTaskEmcalJetHadCorrAA
 * \brief Implementation of a sample jet correlation analysis task with event mixing.
 *
 * This class in an implementation of a sample task for EMCal jet analysis.
 * It derives from AliAnalysisTaskEmcalJet.
 * It performs a simple analysis, producing track, cluster and jet spectra.
 * It also performs a QA of the cluster-track matching.
 * Note: if jets are not used this class can be simplified by deriving
 * from AliAnalysisTaskEmcal and removing the functions DoJetLoop()
 * and AllocateJetHistograms().
 */
class AliAnalysisTaskEmcalJetHadCorrAA : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetHadCorrAA()                                               ;
  AliAnalysisTaskEmcalJetHadCorrAA(const char *name)                               ;
  virtual ~AliAnalysisTaskEmcalJetHadCorrAA()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      FillMixedHistograms(TObjArray* TrkArray, TObjArray* MxTrkArray, int CentBin, int ZvtxBin, int nMix);/// To call DoMixLoop
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        AllocateTrackHistograms()                         ;
  void                        AllocateMixHistograms()                           ;///< Set of histograms for mixing jets

  void                        DoJetLoop()                                       ;
  void                        DoTrackLoop()                                     ;
  void                        DoMixLoop(TObjArray* trkArray, TObjArray* MxtrkArray, int cBin, int zBin, int nMix);/// To fill mixed histograms

  void                        SetupForMixing()                                  ;//Call EventPoolManager - JY
  void                        DoMixing(double cent, double zvertex, TObjArray* trigJets);// doing mixing 


  THistManager                fHistManager                                      ;///< Histogram manager
  AliEventPoolManager*        fPoolMgr                                          ;///< EventPool manager
  TObjArray*                  MxTracksArray                                     ;///< jets array for mixing
  TObjArray*                  JetsArray                                         ;///< tracks array for mixing
  Int_t                       fNzvtxBins                                        ;///< Number of z vertex Bins for AllocateMixHistograms 
  AliPIDResponse*             fPIDResponse                                      ;//! PID response object
  AliMultSelection*           fMultSelection                                    ;//! MultSelection object for centrality
  TH1                        *fHistCentrality_jy;             //!<!event centrality distribution
  TH1                        *fHistQACentrality_jy;           //!<!Centrality QA plot about fEvSelCode - modified by JY
  TH3                        *fHistMixNumOfEventsInPool;             //!<!event centrality distribution
  TH3                        *fHistMixNumOfTracksInPool;           //!<!Centrality QA plot about fEvSelCode - modified by JY
  Bool_t                      doMakeQAplot                                      ;///switch for QA plots 
  AliEventCuts                fEventCut                                         ;///< EventCuts for pileup rejection 

 private:
  AliAnalysisTaskEmcalJetHadCorrAA(const AliAnalysisTaskEmcalJetHadCorrAA&)           ; // not implemented
  AliAnalysisTaskEmcalJetHadCorrAA &operator=(const AliAnalysisTaskEmcalJetHadCorrAA&); // not implemented
  TList *flist; //output  

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetHadCorrAA, 8);
  /// \endcond
};
#endif
