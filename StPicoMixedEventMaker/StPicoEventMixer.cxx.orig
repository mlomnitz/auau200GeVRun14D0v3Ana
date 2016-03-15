#include <limits>

#include "TH3F.h"
#include "THn.h"
#include "StPicoEventMixer.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StEventPlane/StEventPlane.h"
#include "StPicoMixedEventMaker.h"
#include "StMixerEvent.h"
#include "StMixerPair.h"
#include "StMixerTriplet.h"
#include "StD0Hists.h"
#include "StBTofUtil/tofPathLength.hh"

StPicoEventMixer::StPicoEventMixer(int centBin, int vzBin, int psiBin, StEventPlane* eventPlaneMaker, StD0Hists* d0Hists):
   mEvents(NULL), mD0Hists(d0Hists)
{
   mCentBin = centBin;
   mVzBin = vzBin;
   mPsiBin = psiBin;
   mEventPlaneMaker = eventPlaneMaker;
   setEventsBufferSize(1);

}
StPicoEventMixer::~StPicoEventMixer()
{
   for (int i = 0 ; i < mEvents.size() ; i++)
   {
      delete mEvents.at(i);
   }
}
void StPicoEventMixer::finish()
{
   if (!mFirstEvents.size())
   {
      cout << "warning: not enough events to mix!   centBin: " << mCentBin << "  vzBin: " << mVzBin << "  psiBin: " << mPsiBin << "  nEvents: " << mEvents.size() << endl;
      mEventsBufferSize = mEvents.size();
      mixEvents();
   }
   for (int i = 0; i < mEventsBufferSize - 1; i++)
   {
      mEvents.push_back(mFirstEvents.at(0));
      mFirstEvents.erase(mFirstEvents.begin());
      mixEvents();
   }

}
bool StPicoEventMixer::addPicoEvent(StPicoDst const* const picoDst, StThreeVectorF pVertex, float weight)
{
   if (!isGoodEvent(picoDst, pVertex))
      return false;
   int nTracks = picoDst->numberOfTracks();
   StMixerEvent* event = new StMixerEvent(pVertex, picoDst->event()->bField(), mEventPlaneMaker, weight);
   //Event.setNoTracks( nTracks );
   for (int iTrk = 0; iTrk < nTracks; ++iTrk)
   {
      StPicoTrack const* trk = picoDst->track(iTrk);
      bool saveTrack = false;
      bool isPion_ = false;
      bool isKaon_ = false;

      if (!isGoodTrack(trk)  || isCloseTrack(*trk, pVertex)) continue;//good track and Not close trak
      if (isPion(trk, picoDst, pVertex))
      {
         event->addPion(event->getNoTracks());
         isPion_ = true;
         saveTrack = true;
      }
      if (isKaon(trk, picoDst, pVertex))
      {
         event->addKaon(event->getNoTracks());
         isKaon_ = true;
         saveTrack = true;
      }
      if (saveTrack == true)
      {
         StMixerTrack mTrack(pVertex, picoDst->event()->bField(), *trk, isPion_, isKaon_, mEventPlaneMaker->q(iTrk));
         event->addTrack(mTrack);
      }
   }
   //   if (event->getNoPions() > 0 ||  event->getNoKaons() > 0)
   //   {
   mEvents.push_back(event);
   //   }
   /*
   else
   {
   delete event;
   return false;
   }
   */
   //Returns true if need to do mixing, false if buffer has space still
   if (mEvents.size() == mEventsBufferSize)
      return true;
   return false;
}
void StPicoEventMixer::mixEvents()
{
   //-------
   size_t const nEvents = mEvents.size();
   if (!nEvents) return;
   //Template for D0 studies
   for (size_t iEvt2 = 0; iEvt2 < nEvents; iEvt2++)
   {
      int const nTracksEvt1 = mEvents.at(0)->getNoPions();
      int const nTracksEvt2 = mEvents.at(iEvt2)->getNoKaons();

      if (iEvt2 == 0)
      {
         mD0Hists->hCentVzPsiSameEventNoWeight->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10);
         mD0Hists->hCentVzPsiSameEvent->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10, mEvents.at(0)->weight());
      }
      else
      {
         mD0Hists->hCentVzPsiMixedNoWeight->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10);
         mD0Hists->hCentVzPsiMixed->Fill(mCentBin + 0.5, mVzBin * 1.2 - 5.5, (mPsiBin + 0.5)*TMath::Pi() / 10, mEvents.at(0)->weight());
      }
      for (int iTrk2 = 0; iTrk2 < nTracksEvt2; iTrk2++)
      {

         for (int iTrk1 = 0; iTrk1 < nTracksEvt1; iTrk1++)
         {
            if (iEvt2 == 0)
            {
               if (mEvents.at(0)->pionId(iTrk1) == mEvents.at(iEvt2)->kaonId(iTrk2))
                  continue;
            }
            StMixerPair pair(mEvents.at(0)->pionAt(iTrk1), mEvents.at(iEvt2)->kaonAt(iTrk2),
                             mxeCuts::pidMass[mxeCuts::kPion], mxeCuts::pidMass[mxeCuts::kKaon],
                             mEvents.at(0)->vertex(), mEvents.at(iEvt2)->vertex(),
                             mEvents.at(0)->field());
            int charge2 = mEvents.at(0)->pionAt(iTrk1).charge() * mEvents.at(iEvt2)->kaonAt(iTrk2).charge();

            //Topology histos, fill before checking cuts
            if (iEvt2 == 0)
            {
               if (charge2 < 0) mD0Hists->fillSameEvt_US_QADist(pair, mCentBin);
               else mD0Hists->fillSameEvt_LS_QADist(pair, mCentBin);
            }
            else
            {
               if (charge2 < 0) mD0Hists->fillMixedEvtQADist(pair, mCentBin);
            }

            if (!isGoodPair(&pair)) continue;
	    TVector2 QSub = mEvents.at(0)->Q()-mEvents.at(0)->pionAt(iTrk1).q();
	    if(iEvt2 == 0) QSub -= mEvents.at(iEvt2)->kaonAt(iTrk2).q();
	    float dPhi = pair.phi()-QSub.Phi()/2;
	    while(dPhi<0) dPhi += TMath::Pi();
	    while(dPhi>=TMath::Pi()) dPhi -= TMath::Pi();

            double toFill[5] = {mCentBin + 0.5, pair.pt(), pair.eta(), pair.m(), dPhi};

            if (iEvt2 == 0)
            {
               //       if(pair.m()>1.6 && pair.m()<2.1)
               //    cout<<"pair: "<<pair.m()<<" "<<pair.pt()<<" "<<pair.eta()<<" "<<pair.particle1Dca()<<" "<<pair.particle2Dca()<<" "<<pair.dcaDaughters()<<" "<<pair.decayLength()<<" "<<cos(pair.pointingAngle())<<" "<<pair.decayLength() * sin(pair.pointingAngle())<<" "<<dPhi<<endl;
               if (charge2 < 0) mD0Hists->hD0CentPtEtaMDphi->Fill(toFill, mEvents.at(0)->weight());
               else mD0Hists->hD0CentPtEtaMDphiLikeSign->Fill(toFill, mEvents.at(0)->weight());
            }
            else
            {
               if (charge2 < 0) mD0Hists->hD0CentPtEtaMDphiMixed->Fill(toFill, mEvents.at(0)->weight());
               else mD0Hists->hD0CentPtEtaMDphiLikeSignMixed->Fill(toFill, mEvents.at(0)->weight());
            }

	    int iEta = (int)(pair.eta()*10+10);
	    //	    if(iEta<0 || iEta>=20) { cout<<"pair.eta(): "<<pair.eta()<<"   iEta: "<<iEta<<"   etaPion: "<<mEvents.at(0)->pionAt(iTrk1).gMom().pseudoRapidity()<<"  etaKaon: "<<mEvents.at(iEvt2)->kaonAt(iTrk2).gMom().pseudoRapidity()<<endl; continue;}
	    for(int nEtaGaps=0; nEtaGaps<8; nEtaGaps++)
	      {
		TVector2 QSubEtaGap = mEvents.at(0)->QEtaGap(iEta, nEtaGaps);
		int iEta_ = iEta;
		if(iEta_ < nEtaGaps) iEta_ = nEtaGaps-1;
		if(iEta_ > 20-nEtaGaps) iEta_= 20-nEtaGaps;
		int iEtaPion = (int)(mEvents.at(0)->pionAt(iTrk1).gMom().pseudoRapidity()*10+10);
		if(fabs(iEtaPion-iEta_) >= nEtaGaps)
		  QSubEtaGap -= mEvents.at(0)->pionAt(iTrk1).q();
		if(iEvt2 == 0)
		  {
		    int iEtaKaon = (int)(mEvents.at(iEvt2)->kaonAt(iTrk2).gMom().pseudoRapidity()*10+10);
		    if(fabs(iEtaKaon-iEta_) >= nEtaGaps)
		      QSubEtaGap -= mEvents.at(iEvt2)->kaonAt(iTrk2).q();
		  }
		if(QSubEtaGap.Mod()==0) {cout<<"QSubEtaGap.Mod()==0  nEtaGaps: "<<nEtaGaps<<endl; continue;}
		float dPhiEtaGap = pair.phi()-QSubEtaGap.Phi()/2;
		while(dPhiEtaGap<0) dPhiEtaGap += TMath::Pi();
		while(dPhiEtaGap>=TMath::Pi()) dPhiEtaGap -= TMath::Pi();
		double toFill[5] = {mCentBin + 0.5, pair.pt(), pair.m(), dPhiEtaGap, 0.1*nEtaGaps+0.05};
		
		if (iEvt2 == 0)
		  {
		    if (charge2 < 0) mD0Hists->hD0CentPtMDphiEtaGap->Fill(toFill, mEvents.at(0)->weight());
		    else mD0Hists->hD0CentPtMDphiEtaGapLikeSign->Fill(toFill, mEvents.at(0)->weight());
		  }
		else
		  {
		    if (charge2 < 0) mD0Hists->hD0CentPtMDphiEtaGapMixed->Fill(toFill, mEvents.at(0)->weight());
		    else mD0Hists->hD0CentPtMDphiEtaGapLikeSignMixed->Fill(toFill, mEvents.at(0)->weight());
		  }

	      }

         } //second event track loop
      } //first event track loop
   } //loop over second events

   if (mFirstEvents.size() == mEventsBufferSize - 1)
      delete mEvents.at(0);
   else
      mFirstEvents.push_back(mEvents.at(0));
   mEvents.erase(mEvents.begin());

}
bool StPicoEventMixer::isGoodEvent(StPicoDst const * const picoDst, StThreeVectorF pVertex)
{
   StPicoEvent* picoEvent = picoDst->event();
   return ((picoEvent->triggerWord() & mxeCuts::triggerWord) &&
           fabs(pVertex.z()) < mxeCuts::maxVz &&
           fabs(pVertex.z() - picoEvent->vzVpd()) < mxeCuts::vzVpdVz &&
           sqrt(TMath::Power(pVertex.x(), 2) + TMath::Power(pVertex.y(), 2)) <= mxeCuts:: Vrcut);
}
bool StPicoEventMixer::isKaon(const StPicoTrack* trk, const StPicoDst* picoDst, StThreeVectorF pVertex)
{
   if (!isTpcKaon(trk)) return false;
   float beta = getTofBeta(trk, picoDst, pVertex);
   if (beta < 0) return true;
   float p = trk->gMom(pVertex, picoDst->event()->bField()).mag();
   float mKaon = 0.493677;
   float oneOverBetaExpected = sqrt(mKaon * mKaon / p / p + 1);
   //  hOneOverBetaDiffPionP->Fill(pPion, 1./betaPion-oneOverBetaExpectedPion);
   if (fabs(1. / beta - oneOverBetaExpected) > mxeCuts::tofOneOverBetaDiffPion) return false;
   return true;
}
bool StPicoEventMixer::isPion(const StPicoTrack* trk, const StPicoDst* picoDst, StThreeVectorF pVertex)
{
   if (!isTpcPion(trk)) return false;
   float beta = getTofBeta(trk, picoDst, pVertex);
   if (beta < 0) return true;
   float p = trk->gMom(pVertex, picoDst->event()->bField()).mag();
   float mPion = 0.13957;
   float oneOverBetaExpected = sqrt(mPion * mPion / p / p + 1);
   //  hOneOverBetaDiffPionP->Fill(pPion, 1./betaPion-oneOverBetaExpectedPion);
   if (fabs(1. / beta - oneOverBetaExpected) > mxeCuts::tofOneOverBetaDiffKaon) return false;
   return true;
}
bool StPicoEventMixer::isTpcPion(StPicoTrack const * const trk)
{
   return (fabs(trk->nSigmaPion()) < mxeCuts::nSigmaPion);
}
bool StPicoEventMixer::isTpcKaon(StPicoTrack const * const trk)
{
   return (fabs(trk->nSigmaKaon()) < mxeCuts::nSigmaKaon);
}
bool StPicoEventMixer::isGoodTrack(StPicoTrack const * const trk)
{
   return ((!mxeCuts::mRequireHft || trk->isHFTTrack()) &&
           trk->nHitsFit() >= mxeCuts::nHitsFit && trk->gPt() > mxeCuts::minPt);
}
bool StPicoEventMixer::isCloseTrack(StPicoTrack const& trk, StThreeVectorF const& pVtx)
{
   StPhysicalHelixD helix = trk.dcaGeometry().helix();
   helix.moveOrigin(helix.pathLength(pVtx));
   if ((helix.origin() - pVtx).mag() > mxeCuts::dca2pVtx) return false;
   return true;
}
bool StPicoEventMixer::isGoodPair(StMixerPair const& pair)
{
   int ptIndex = getD0PtIndex(pair);
   return (pair.m() > mxeCuts::massMin && pair.m() < mxeCuts::massMax &&
           pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
           pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
           pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
           std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
           ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex]);
}
//-----------------------------------------------------------------------------
int StPicoEventMixer::getD0PtIndex(StMixerPair const& pair) const
{
   for (int i = 0; i < mxeCuts::nPtBins; i++)
   {
      if ((pair.pt() >= mxeCuts::PtEdge[i]) && (pair.pt() < mxeCuts::PtEdge[i + 1]))
         return i;
   }
   return mxeCuts::nPtBins - 1;
}
float StPicoEventMixer::getTofBeta(const StPicoTrack* trk, const StPicoDst* picoDst, StThreeVectorF pVertex) const
{
   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0)
   {
      StPicoBTofPidTraits *tofPid = picoDst->btofPidTraits(index2tof);

      if (tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            StPhysicalHelixD helix = trk->helix();
            float L = tofPathLength(&pVertex, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (2.99792458e10 / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
