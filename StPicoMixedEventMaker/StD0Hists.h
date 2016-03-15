#ifndef StD0Hists__h
#define StD0Hists__h

/* **************************************************
 *  A class to create and save my D0 analysis histograms.
 *
 *  Authors: Mustafa Mustafa (mmustafa@lbl.gov)
 *
 * **************************************************
 */

#include <string>
#include "TObject.h"

class TH1D;
class TH1F;
class TH2F;
class TH3F;
class THn;
class TProfile;
class TFile;
class TString;

#include "StMixerPair.h"
#include "StMixerCuts.h"

class StD0Hists: public TObject
{
public:
   StD0Hists(std::string fileBaseName);
   virtual ~StD0Hists();
   void closeFile();
   void fillSameEvt_US_QADist(StMixerPair const&, int const);
   void fillSameEvt_LS_QADist(StMixerPair const&, int const);
   void fillMixedEvtQADist(StMixerPair const&, int const);
   //event level qa
   TH1D*      hTotalNumberOfEvents;
   TH2F*      hVzVpdVz;
   TH1F*      hVzDiff;
   TH2F*      hVxy;
   TH1F*      hRefMult;
   TH1F*      hGRefMult;
   TH1F*      hTrigger;
   TH1F*      hCentrality;
   TH1F*      hCentralityWeighted;

   TH2F*      hOneOverBetaDiffKaonP;
   TH2F*      hOneOverBetaDiffPionP;

   //mix event counts
   TH3F* hCentVzPsi;
   TH3F* hCentVzPsiSameEventNoWeight;
   TH3F* hCentVzPsiMixedNoWeight;
   TH3F* hCentVzPsiSameEvent;
   TH3F* hCentVzPsiMixed;

   //d0 v2 histograms
   THn* hD0CentPtEtaMDphi;
   THn* hD0CentPtEtaMDphiLikeSign;
   THn* hD0CentPtEtaMDphiMixed;
   THn* hD0CentPtEtaMDphiLikeSignMixed;

   //add daughter pT
   THn* hD0CentPtEtaMDphiDaug;
   THn* hD0CentPtEtaMDphiDaugLikeSign;
   THn* hD0CentPtEtaMDphiDaugMixed;
   THn* hD0CentPtEtaMDphiDaugLikeSignMixed;

   THn* hD0CentPtMDphiEtaGap;
   THn* hD0CentPtMDphiEtaGapLikeSign;
   THn* hD0CentPtMDphiEtaGapMixed;
   THn* hD0CentPtMDphiEtaGapLikeSignMixed;

   //Topology distribution
   TH3F* mSE_Vtx;
   TH3F* mME_Vtx;
   TH3F* mSE_LS;
   TH3F* mSE_US;
   TH3F* mME_LS;
   TH3F* mME_US;
   //QA hists
   THn*  mSE_US_DecayTopology;
   TH3F* mSE_US_PointingAngle;
   TH3F* mSE_US_DecayL;
   TH3F* mSE_US_Dca12;
   TH3F* mSE_US_PionDca2Vtx;
   TH3F* mSE_US_KaonDca2Vtx;
   TH3F* mSE_US_D0Dca2Vtx;
   //
   TH3F* mSE_LS_PointingAngle;
   TH3F* mSE_LS_DecayL;
   TH3F* mSE_LS_Dca12;
   TH3F* mSE_LS_PionDca2Vtx;
   TH3F* mSE_LS_KaonDca2Vtx;
   TH3F* mSE_LS_D0Dca2Vtx;
   //
   THn*  mME_US_DecayTopology;
   TH3F* mME_US_PointingAngle;
   TH3F* mME_US_DecayL;
   TH3F* mME_US_Dca12;
   TH3F* mME_US_PionDca2Vtx;
   TH3F* mME_US_KaonDca2Vtx;
   TH3F* mME_US_D0Dca2Vtx;
private:

   ClassDef(StD0Hists, 1);
};

#endif
