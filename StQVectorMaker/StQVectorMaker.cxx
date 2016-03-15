#include "StQVectorMaker.h"
#include "StRoot/StPicoDstMaker/StPicoDst.h"
#include "StRoot/StPicoDstMaker/StPicoEvent.h"
#include "StRoot/StPicoDstMaker/StPicoTrack.h"
#include "StRoot/StPicoDstMaker/StPicoDstMaker.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"
#include "StRoot/StRefMultCorr/CentralityMaker.h"
#include "PhysicalConstants.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "TProfile.h"

ClassImp(StQVectorMaker)

//-----------------------------------------------------------------------------
StQVectorMaker::StQVectorMaker(const char* name, StPicoDstMaker *picoMaker)
: StMaker(name)
{
    mPicoDstMaker = picoMaker;
    mPicoDst = 0;
}

//-----------------------------------------------------------------------------
StQVectorMaker::~StQVectorMaker()
{ /*  */ }

//-----------------------------------------------------------------------------
Int_t StQVectorMaker::Init()
{

    mAcceptEvent = false;

    mRefMultCorr = new StRefMultCorr("grefmult");

    //Event Cuts 
    mVzMax = 6.0;
    mDeltaVzMax = 3.0;

    //Track Cuts
    mNHitsFitMin = 15;
    mNHitsFitRatioMin = 0.52;
    mEtaMax = 1.0;
    mPtMin = 0.15;
    mPtMax = 2.;
    mDcaMax = 3.0;

    mFileOut = new TFile(mOutputName, "recreate");

    // event level QA
    hVzVpdVz = new TH2F("hVzVpdVz","hVzVpdVz",200,-100,100,200,-100,100);
    hVzDiff  = new TH1F("hVzDiff","hVzDiff",500,-100,100);
    hVxy = new TH2F("hVxy","hVxy",500,-1,1,500,-1,1);
    hRefMult = new TH1I("hRefMult","hRefMult",1000,0,1000);
    hGRefMult = new TH1I("hGRefMult","hGRefMult",1000,0,1000);
    hTrigger = new TH1I("hTrigger","hTrigger",32,0,32);
    hCentrality = new TH1I("hCentrality", "hCentrality", 9, 0, 9);

    // track level QA
    hNHitsFit = new TH1I("hNHitsFit", "hNHitsFit", 50, 0, 50);
    hDca = new TH1F("hDca", "hDca", 20000, -10, 10);
    hEta = new TH1F("hEta", "hEta", 30, -1.5, 1.5);
    hPt = new TH1F("hPt", "hPt", 200, 0, 10);
    
    // event plane and Q vector
    float PI = TMath::Pi();

    hPhiCentEtaPlusZPlus = new TH2F("hPhiCentEtaPlusZPlus","hPhiCentEtaPlusZPlus",9,0,9,120,-PI,PI);
    hPhiCentEtaPlusZMinus = new TH2F("hPhiCentEtaPlusZMinus","hPhiCentEtaPlusZMinus",9,0,9,120,-PI,PI);
    hPhiCentEtaMinusZPlus = new TH2F("hPhiCentEtaMinusZPlus","hPhiCentEtaMinusZPlus",9,0,9,120,-PI,PI);
    hPhiCentEtaMinusZMinus = new TH2F("hPhiCentEtaMinusZMinus","hPhiCentEtaMinusZMinus",9,0,9,120,-PI,PI);

    prfQxCentEtaPlus = new TProfile("prfQxCentEtaPlus","prfQxCentEtaPlus",9,0,9);
    prfQyCentEtaPlus = new TProfile("prfQyCentEtaPlus","prfQyCentEtaPlus",9,0,9);
    prfQxCentEtaMinus = new TProfile("prfQxCentEtaMinus","prfQxCentEtaMinus",9,0,9);
    prfQyCentEtaMinus = new TProfile("prfQyCentEtaMinus","prfQyCentEtaMinus",9,0,9);

    hEventPlaneCent = new TH2F("hEventPlaneCent","hEventPlaneCent",9,0,9,60,0,PI);
    hQyQxCent = new TH3F("hQyQxCent", "hQyQxCent", 9,0,9,1000,-50,50,1000,-50,50);

    return kStOK;
}

//-----------------------------------------------------------------------------
Int_t StQVectorMaker::Finish() {
    mFileOut->Write();
    return kStOK;
}

//-----------------------------------------------------------------------------
void StQVectorMaker::Clear(Option_t *opt) {
}

//-----------------------------------------------------------------------------
Int_t StQVectorMaker::Make() {
    if(!mPicoDstMaker) {
        LOG_WARN << " No PicoDstMaker! Skip! " << endm;
        return kStWarn;
    }

    mPicoDst = mPicoDstMaker->picoDst();
    if(!mPicoDst) {
        LOG_WARN << " No PicoDst! Skip! " << endm;
        return kStWarn;
    }

    getEventInfo();//get event info
    if(mAcceptEvent){
        getTrackInfo();
    }

    return kStOK;
}

/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StQVectorMaker::setOutputName(Char_t* dir, Char_t* file)
{
    TString dirName(dir);
    TString fileName(file);
    mOutputName = dirName+"/"+fileName+".qVector.root";	
    
    return;
}
/*--------------------------------------------------------------------------------------------------------------------------------------------------------*/
void StQVectorMaker::getEventInfo()
{
    mAcceptEvent = false;
    
    if(!mPicoDst) return;
    
    //Load event
    mPicoEvent = (StPicoEvent*)mPicoDst->event();
    if(!mPicoEvent){
        cerr<<"Error opening picoDst Event, skip!"<<endl;
        return;
    }
    
    for(int i=0; i<32; i++)
        if(mPicoEvent->triggerWord()>>i & 0x1)
            hTrigger->Fill(i);

    bool isMinBias=kFALSE;
    for(int i=0;i<11;i++) { if(mPicoEvent->triggerWord() & (1<<i)) isMinBias=kTRUE ;}//Select MB trigger   
    //if (!(isMinBias)) {cout<<"not a mb trigger"<<endl;return 0;}
    bool isVPDMB5=kFALSE;
    for(int i=0;i<9;i++) { if(mPicoEvent->triggerWord() & (1<<i)) isVPDMB5=kTRUE ;}//Select MB trigger   
    if (!(isVPDMB5)) {
	//cout<<"not a VPDmb trigger"<<endl;
	return;
    }
    
    //Remove bad vertices
    mVertexPos = mPicoEvent->primaryVertex();
    
    hVzVpdVz->Fill(mVertexPos.z(), mPicoEvent->vzVpd());
    hVzDiff->Fill(mPicoEvent->vzVpd() - mVertexPos.z());
    hVxy->Fill(mVertexPos.y(), mVertexPos.x());

    if(TMath::Abs(mVertexPos.z()) > mVzMax) return;
    if(TMath::Abs(mVertexPos.z() - mPicoEvent->vzVpd()) > mDeltaVzMax) return;

    hRefMult->Fill(mPicoEvent->refMult());
    hGRefMult->Fill(mPicoEvent->grefMult());

    if(!mRefMultCorr) {  
      LOG_WARN << " No mRefMultCorr! Skip! " << endl;
      return;
    } 
    mRefMultCorr->init(mPicoDst->event()->runId());
    mRefMultCorr->initEvent(mPicoDst->event()->grefMult(),mVertexPos.z(),mPicoDst->event()->ZDCx()) ;
    mCent  = mRefMultCorr->getCentralityBin9();

    hCentrality->Fill(mCent);

    mAcceptEvent = true;

    mBField = mPicoEvent->bField();
}

/*----------------------------------------------------------------------------------------------------------------------*/
void StQVectorMaker::getTrackInfo()
{
    float vertexZ = mPicoEvent->primaryVertex().z();
    float Qx=0., Qy=0.;

    //Load tracks for default vertex index
    for(int iTrack=0; iTrack<mPicoDst->numberOfTracks(); iTrack++)
        {
            StPicoTrack* picoTrack =(StPicoTrack*) mPicoDst->track(iTrack);
            if(!picoTrack){
                break;
            }
            
            hNHitsFit->Fill(picoTrack->nHitsFit());
            if(picoTrack->nHitsFit() < mNHitsFitMin) continue;
	    if(1.*picoTrack->nHitsFit()/picoTrack->nHitsMax() < mNHitsFitRatioMin) continue;

	    StPhysicalHelix* helix = &picoTrack->dcaGeometry().helix();
	    float dca = helix->geometricSignedDistance(mVertexPos);
            hDca->Fill(dca);
            if(TMath::Abs(dca) > mDcaMax) continue;
            
	    float pathLengthToPrimaryVertex =helix->pathLength(mVertexPos.x(), mVertexPos.y());
	    StThreeVectorF momentum = helix->momentumAt(pathLengthToPrimaryVertex, mBField*kilogauss);
            float pt = momentum.perp();
            float eta = momentum.pseudoRapidity();
            float phi = momentum.phi();
            hEta->Fill(eta);
            hPt->Fill(pt);
            if(fabs(eta) > mEtaMax) continue;
            if(pt<mPtMin || pt>mPtMax) continue;

            if(eta>0 && vertexZ>0) hPhiCentEtaPlusZPlus->Fill(mCent, phi);
            if(eta>0 && vertexZ<0) hPhiCentEtaPlusZMinus->Fill(mCent, phi);
            if(eta<0 && vertexZ>0) hPhiCentEtaMinusZPlus->Fill(mCent, phi);
            if(eta<0 && vertexZ<0) hPhiCentEtaMinusZMinus->Fill(mCent, phi);

            float qx = cos(3*phi)*pt;
            float qy = sin(3*phi)*pt;

            if(eta>0)
                {
                    prfQxCentEtaPlus->Fill(mCent, qx);
                    prfQyCentEtaPlus->Fill(mCent, qy);
                }
            else
                {
                    prfQxCentEtaMinus->Fill(mCent, qx);
                    prfQyCentEtaMinus->Fill(mCent, qy);
                }

            Qx += qx;
            Qy += qy;
        }//loop thru picoTracks
    TVector2 Q(Qx, Qy);
    float eventPlane = Q.Phi()/3.0;
    hQyQxCent->Fill(mCent, Qx, Qy);
    if(Q.Mod()>0)
      hEventPlaneCent->Fill(mCent, eventPlane);
}

int StQVectorMaker::Centrality(int gRefMult)
{
    int centrality;
    int centFull[9] = {8,18,35,62,103,161,240,347,415};
    if      (gRefMult>=centFull[8]) centrality=0;
    else if (gRefMult>=centFull[7]) centrality=1;
    else if (gRefMult>=centFull[6]) centrality=2;
    else if (gRefMult>=centFull[5]) centrality=3;
    else if (gRefMult>=centFull[4]) centrality=4;
    else if (gRefMult>=centFull[3]) centrality=5;
    else if (gRefMult>=centFull[2]) centrality=6;
    else if (gRefMult>=centFull[1]) centrality=7;
    else if (gRefMult>=centFull[0]) centrality=8;
    else centrality = 9;
    
    return centrality;
}

