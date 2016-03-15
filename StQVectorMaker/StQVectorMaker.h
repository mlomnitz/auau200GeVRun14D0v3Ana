#ifndef STAR_StQVectorMaker
#define STAR_StQVectorMaker
#include "StMaker.h"
#include "StThreeVectorF.hh"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoDstMaker;
class StRefMultCorr;
class TH1I;
class TH1F;
class TH2F;
class TH3F;
class TProfile;

class StQVectorMaker : public StMaker {
  public:
    StQVectorMaker(const char *name, StPicoDstMaker *picoMaker);
    virtual ~StQVectorMaker();
    
    virtual Int_t Init();
    virtual Int_t Make();
    virtual void  Clear(Option_t *opt="");
    virtual Int_t Finish();
    
    void getEventInfo();
    void getTrackInfo();
    void setOutputName(Char_t* dir=".", Char_t* name="test");
    int Centrality(int gRefMult);
  private:
    StPicoDstMaker *mPicoDstMaker;
    StPicoDst      *mPicoDst;
    StPicoEvent	   *mPicoEvent;
    StRefMultCorr  *mRefMultCorr;

    bool	mAcceptEvent;
    int         mCent;
    float       mBField;
    StThreeVectorF mVertexPos;

    //cuts
    Float_t mVzMax;
    Float_t mRefMultMin;
    Float_t mDeltaVzMax;
    Float_t mPtMin;
    Float_t mPtMax;
    UShort_t mNHitsFitMin;
    Float_t mNHitsFitRatioMin;
    Float_t mDcaMax;
    Float_t mEtaMax;

    TString mOutputName;
    TFile* mFileOut;
        
    //event level qa
    TH2F*      hVzVpdVz;
    TH1F*      hVzDiff;
    TH2F*      hVxy;
    TH1I*      hRefMult;
    TH1I*      hGRefMult;
    TH1I*      hTrigger;
    TH1I*      hCentrality;

    //track level qa
    TH1I*      hNHitsFit;
    TH1F*      hDca;
    TH1F*      hEta;
    TH1F*      hPt;

    //event plane and Q vector
    TH2F*      hPhiCentEtaPlusZPlus;
    TH2F*      hPhiCentEtaPlusZMinus;
    TH2F*      hPhiCentEtaMinusZPlus;
    TH2F*      hPhiCentEtaMinusZMinus;

    TProfile*  prfQxCentEtaPlus;
    TProfile*  prfQyCentEtaPlus;
    TProfile*  prfQxCentEtaMinus;
    TProfile*  prfQyCentEtaMinus;

    TH2F*      hEventPlaneCent;
    TH3F*      hQyQxCent;

    ClassDef(StQVectorMaker, 1)
};

#endif
