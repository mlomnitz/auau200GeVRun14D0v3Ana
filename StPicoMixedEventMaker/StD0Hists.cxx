#include <cmath>

#include "StD0Hists.h"
#include "StMixerPair.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THn.h"
#include "TProfile.h"

ClassImp(StD0Hists);

StD0Hists::StD0Hists(std::string fileBaseName = "")
{
   // event level QA
   hTotalNumberOfEvents = new TH1D("hTotalNumberOfEvents","hTotalNumberOfEvents",1,0,1);
   hVzVpdVz = new TH2F("hVzVpdVz", "hVzVpdVz", 200, -100, 100, 200, -100, 100);
   hVzDiff  = new TH1F("hVzDiff", "hVzDiff", 500, -100, 100);
   hVxy = new TH2F("hVxy", "hVxy", 500, -1, 1, 500, -1, 1);
   hRefMult = new TH1F("hRefMult", "hRefMult", 1000, 0, 1000);
   hGRefMult = new TH1F("hGRefMult", "hGRefMult", 1000, 0, 1000);
   hTrigger = new TH1F("hTrigger", "hTrigger", 32, 0, 32);
   hCentrality = new TH1F("hCentrality", "hCentrality", 9, 0, 9);
   hCentralityWeighted = new TH1F("hCentralityWeighted", "hCentralityWeighted", 9, 0, 9);

   // track level QA
   hOneOverBetaDiffKaonP = new TH2F("hOneOverBetaDiffKaonP", "hOneOverBetaDiffKaonP", 50, 0, 5, 100, -0.5, 0.5);
   hOneOverBetaDiffPionP = new TH2F("hOneOverBetaDiffPionP", "hOneOverBetaDiffPionP", 50, 0, 5, 100, -0.5, 0.5);

   //mix event counts
   hCentVzPsi = new TH3F("hCentVzPsi", "hCentVzPsi", 9, 0, 9, 10, -6, 6, 10, 0, TMath::Pi());
   hCentVzPsiSameEventNoWeight = new TH3F("hCentVzPsiSameEventNoWeight", "hCentVzPsiSameEventNoWeight", 9, 0, 9, 10, -6, 6, 10, 0, TMath::Pi());
   hCentVzPsiMixedNoWeight = new TH3F("hCentVzPsiMixedNoWeight", "hCentVzPsiMixedNoWeight", 9, 0, 9, 10, -6, 6, 10, 0, TMath::Pi());
   hCentVzPsiSameEvent = new TH3F("hCentVzPsiSameEvent", "hCentVzPsiSameEvent", 9, 0, 9, 10, -6, 6, 10, 0, TMath::Pi());
   hCentVzPsiMixed = new TH3F("hCentVzPsiMixed", "hCentVzPsiMixed", 9, 0, 9, 10, -6, 6, 10, 0, TMath::Pi());

   //D0 histograms
   const int nDim = 5;
   int nBins[nDim] = {9, 100, 20, 250, 10};//cent, pt, eta, m, dPhi
   double xMin[nDim] = {0, 0, -1, 0, 0};
   double xMax[nDim] = {9, 10, 1, 2.5, TMath::Pi()};
   hD0CentPtEtaMDphi = new THnF("hD0CentPtEtaMDphi", "hD0CentPtEtaMDphi", nDim, nBins, xMin, xMax);
   hD0CentPtEtaMDphiLikeSign = new THnF("hD0CentPtEtaMDphiLikeSign", "hD0CentPtEtaMDphiLikeSign", nDim, nBins, xMin, xMax);
   hD0CentPtEtaMDphiMixed = new THnF("hD0CentPtEtaMDphiMixed", "hD0CentPtEtaMDphiMixed", nDim, nBins, xMin, xMax);
   hD0CentPtEtaMDphiLikeSignMixed = new THnF("hD0CentPtEtaMDphiLikeSignMixed", "hD0CentPtEtaMDphiLikeSignMixed", nDim, nBins, xMin, xMax);

   //D0 histograms
   const int nDimDaug = 5;
   int nBinsDaug[nDimDaug] = {9, 100, 10, 250, 10};//cent, pt, daughterpt1, m, daughterpt2
   double xMinDaug[nDimDaug] = {0, 0, 0.6, 0, 0.6};
   double xMaxDaug[nDimDaug] = {9, 10, 1.6, 2.5, 1.6};
   hD0CentPtEtaMDphiDaug = new THnF("hD0CentPtEtaMDphiDaug", "hD0CentPtEtaMDphiDaug", nDimDaug, nBinsDaug, xMinDaug, xMaxDaug);
   hD0CentPtEtaMDphiDaugLikeSign = new THnF("hD0CentPtEtaMDphiDaugLikeSign", "hD0CentPtEtaMDphiDaugLikeSign", nDimDaug, nBinsDaug, xMinDaug, xMaxDaug);
   hD0CentPtEtaMDphiDaugMixed = new THnF("hD0CentPtEtaMDphiDaugMixed", "hD0CentPtEtaMDphiDaugMixed", nDimDaug, nBinsDaug, xMinDaug, xMaxDaug);
   hD0CentPtEtaMDphiDaugLikeSignMixed = new THnF("hD0CentPtEtaMDphiDaugLikeSignMixed", "hD0CentPtEtaMDphiDaugLikeSignMixed", nDimDaug, nBinsDaug, xMinDaug, xMaxDaug);

   const int nDim2 = 5;
   int nBins2[nDim2] = {9, 100, 50, 10, 8};//cent, pt, m, dPhi, etaGap
   double xMin2[nDim2] = {0, 0, 1.6, 0, 0};
   double xMax2[nDim2] = {9, 10, 2.1, TMath::Pi(), 0.8};
   hD0CentPtMDphiEtaGap = new THnF("hD0CentPtMDphiEtaGap", "hD0CentPtMDphiEtaGap", nDim2, nBins2, xMin2, xMax2);
   hD0CentPtMDphiEtaGapLikeSign = new THnF("hD0CentPtMDphiEtaGapLikeSign", "hD0CentPtMDphiEtaGapLikeSign", nDim2, nBins2, xMin2, xMax2);
   hD0CentPtMDphiEtaGapMixed = new THnF("hD0CentPtMDphiEtaGapMixed", "hD0CentPtMDphiEtaGapMixed", nDim2, nBins2, xMin2, xMax2);
   hD0CentPtMDphiEtaGapLikeSignMixed = new THnF("hD0CentPtMDphiEtaGapLikeSignMixed", "hD0CentPtMDphiEtaGapLikeSignMixed", nDim2, nBins2, xMin2, xMax2);

   /*
   const int nDimMixed=7;
   //  int nBinsMixed[nDimMixed] = {9, 10, 10, 10, 20, 50, 10};//cent, vZ, psi, pt, eta, m, dPhi
   int nBinsMixed[nDimMixed] = {9, 10, 10, 10, 20, 5, 10};//cent, vZ, psi, pt, eta, m, dPhi
   double xMinMixed[nDimMixed] = {0, -6, 0, 0, -1, 1.6, 0};
   double xMaxMixed[nDimMixed] = {9, 6, TMath::Pi(), 10, 1, 2.1, TMath::Pi()};
   hD0CentVzPsiPtEtaMDphiMixed = new THnF("hD0CentVzPsiPtEtaMDphiMixed", "hD0CentVzPsiPtEtaMDphiMixed", nDimMixed, nBinsMixed, xMinMixed, xMaxMixed);
   hD0CentVzPsiPtEtaMDphiLikeSignMixed = new THnF("hD0CentVzPsiPtEtaMDphiLikeSignMixed", "hD0CentVzPsiPtEtaMDphiLikeSignMixed", nDimMixed, nBinsMixed, xMinMixed, xMaxMixed);
   */

   mSE_LS = new TH3F(Form("%s_se_ls_mass", fileBaseName.c_str()), "Same Event LS pair Invariant mass(K#pi);p_{T}(K#pi)(GeV/c);centrality;Mass_{K#pi}(GeV/c^{2})", 150, 0, 15, 9, 0, 9, 250, 0, 2.5);
   mSE_US = new TH3F(Form("%s_se_us_mass", fileBaseName.c_str()), "Same Event US pair Invariant mass(K#pi);p_{T}(K#pi)(GeV/c);centrality;Mass_{K#pi}(GeV/c^{2})", 150, 0, 15, 9, 0, 9, 250, 0, 2.5);
   mME_LS = new TH3F(Form("%s_me_ls_mass", fileBaseName.c_str()), "Mixed Event LS pair Invariant mass(K#pi);p_{T}(K#pi)(GeV/c);centrality;Mass_{K#pi}(GeV/c^{2})", 150, 0, 15, 9, 0, 9, 250, 0, 2.5);
   mME_US = new TH3F(Form("%s_me_us_mass", fileBaseName.c_str()), "Mixed Event US pair Invariant mass(K#pi);p_{T}(K#pi)(GeV/c);centrality;Mass_{K#pi}(GeV/c^{2})", 150, 0, 15, 9, 0, 9, 250, 0, 2.5);

   float const maxDca = 0.3;
   int const nDcaBins = 300;

   float const maxDca12 = 0.1;
   int const nDca12Bins = 200;

   int const decayTopologyNDim = 5;
   int const decayTopologyBins[decayTopologyNDim] = {16,14,6,6,95}; // pT, d0Dca2Vtx, pDca, kDca, decayLanegth
   double const decayTopologyXMin[decayTopologyNDim] = {0, 0.0050, 0.0050, 0.0050, 0.0050}; 
   double const decayTopologyXMax[decayTopologyNDim] = {8, 0.0120, 0.0110, 0.0110, 0.1000}; 

   //QA Foreground
   mSE_US_DecayTopology = new THnF("mSE_US_DecayTopology","mSE_US_DecayTopology;p_{T}(GeV/c);d0Dca2Vtx(cm);pDca(cm);kDca(cm);decayL(cm)",
                                                           decayTopologyNDim,decayTopologyBins,decayTopologyXMin,decayTopologyXMax);
   mSE_US_PointingAngle = new TH3F(Form("%s_se_us_pointingangle", fileBaseName.c_str()), "Same Event US pointing angle; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
   mSE_US_DecayL = new TH3F(Form("%s_se_us_decayL", fileBaseName.c_str()), "Same Event US Decay Length; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mSE_US_Dca12 = new TH3F(Form("%s_se_us_dcaDaughters", fileBaseName.c_str()), "Same Event US dca daughters; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDca12Bins, 0, maxDca12);
   mSE_US_PionDca2Vtx = new TH3F(Form("%s_se_us_pionDca", fileBaseName.c_str()), "Same Event #pi dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mSE_US_KaonDca2Vtx = new TH3F(Form("%s_se_us_kaonDca", fileBaseName.c_str()), "Same Event US K dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mSE_US_D0Dca2Vtx = new TH3F(Form("%s_se_us_D0Dca2Vtx", fileBaseName.c_str()), "SameEvent US D0 dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);
   //
   mSE_LS_PointingAngle = new TH3F(Form("%s_se_ls_pointingangle", fileBaseName.c_str()), "Same Event LS pointing angle; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
   mSE_LS_DecayL = new TH3F(Form("%s_se_ls_decayL", fileBaseName.c_str()), "Same Event LS Decay Length; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mSE_LS_Dca12 = new TH3F(Form("%s_se_ls_dcaDaughters", fileBaseName.c_str()), "Same Event LS dca daughters; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDca12Bins, 0, maxDca12);
   mSE_LS_PionDca2Vtx = new TH3F(Form("%s_se_ls_pionDca", fileBaseName.c_str()), "Same Event #pi dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mSE_LS_KaonDca2Vtx = new TH3F(Form("%s_se_ls_kaonDca", fileBaseName.c_str()), "Same Event LS K dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mSE_LS_D0Dca2Vtx = new TH3F(Form("%s_se_ls_D0Dca2Vtx", fileBaseName.c_str()), "SameEvent LS D0 dca 2 vertex; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);
   //
   mME_US_DecayTopology = new THnF("mME_US_DecayTopology","mME_US_DecayTopology;p_{T}(GeV/c);d0Dca2Vtx(cm);pDca(cm);kDca(cm);decayL(cm)",
                                                           decayTopologyNDim,decayTopologyBins,decayTopologyXMin,decayTopologyXMax);
   mME_US_PointingAngle = new TH3F(Form("%s_me_us_pointingangle", fileBaseName.c_str()), "Same Event US pointing angle ; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
   mME_US_DecayL = new TH3F(Form("%s_me_us_decayL", fileBaseName.c_str()), "Same Event US Decay Length ; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mME_US_Dca12 = new TH3F(Form("%s_me_us_dcaDaughters", fileBaseName.c_str()), "Same Event US dca daughters ; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDca12Bins, 0, maxDca12);
   mME_US_PionDca2Vtx = new TH3F(Form("%s_me_us_pionDca", fileBaseName.c_str()), "Same Event #pi dca 2 vertex ; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mME_US_KaonDca2Vtx = new TH3F(Form("%s_me_us_kaonDca", fileBaseName.c_str()), "Same Event US K dca 2 vertex ; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
   mME_US_D0Dca2Vtx = new TH3F(Form("%s_me_us_D0Dca2Vtx", fileBaseName.c_str()), "SameEvent US D0 dca 2 vertex ; p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);

}
StD0Hists::~StD0Hists()
{
}
void StD0Hists::closeFile()
{

   return;
}
// QA histogram filling
// --------------------------------------
void StD0Hists::fillMixedEvtQADist(StMixerPair const&  pair, int const centrality)
{
   int ptIndex;
   for (int i = 0; i < mxeCuts::nPtBins; i++)
   {
      if ((pair.pt() >= mxeCuts::PtEdge[i]) && (pair.pt() < mxeCuts::PtEdge[i + 1]))
      {
         ptIndex = i;
         break;
      }
   }
   if (pair.m() <  mxeCuts::QAmassMin || pair.m() > mxeCuts::QAmassMax) return;

   // Decay Topology
   if (centrality >= mxeCuts::decayTopologyCentrality &&
       pair.particle1Dca() > mxeCuts::decayTopologyMinDca && pair.particle2Dca() > mxeCuts::decayTopologyMinDca &&
       pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
       pair.decayLength() > mxeCuts::decayTopologyMinDca &&
       std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
       ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::decayTopologyMaxD0Dca2Vtx)
   {
     double toFill[] = {pair.pt(), pair.decayLength() * std::sin(pair.pointingAngle()), 
                        pair.particle1Dca(), pair.particle2Dca(), pair.decayLength()};
     mME_US_DecayTopology->Fill(toFill);
   }

   //Cos theta
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         //std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mME_US_PointingAngle->Fill(pair.pt(), centrality,  std::cos(pair.pointingAngle()));

   //DecayL
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         //pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mME_US_DecayL->Fill(pair.pt(), centrality, pair.decayLength());

   //DcaDaughter
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         //pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mME_US_Dca12->Fill(pair.pt(), centrality, pair.dcaDaughters());

   //PionDca
   if (//pair.particle1Dca() > mxeCuts::pDca[ptIndex] &&
      pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
      pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
      pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
      std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
      ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mME_US_PionDca2Vtx->Fill(pair.pt(), centrality, pair.particle1Dca());

   //Kaon Dca
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] &&
         //pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mME_US_KaonDca2Vtx->Fill(pair.pt(), centrality,  pair.particle2Dca());

   //D0 dca
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex]
      )
      mME_US_D0Dca2Vtx->Fill(pair.pt(), centrality, (pair.decayLength()) * sin(pair.pointingAngle()));
}
// --------------------------------------
void StD0Hists::fillSameEvt_US_QADist(StMixerPair const&  pair, int const centrality)
{
   int ptIndex;
   for (int i = 0; i < mxeCuts::nPtBins; i++)
   {
      if ((pair.pt() >= mxeCuts::PtEdge[i]) && (pair.pt() < mxeCuts::PtEdge[i + 1]))
      {
         ptIndex = i;
         break;
      }
   }
   if (pair.m() <  mxeCuts::QAmassMin || pair.m() > mxeCuts::QAmassMax) return;

   // Decay Topology
   if (centrality >= mxeCuts::decayTopologyCentrality &&
       pair.particle1Dca() > mxeCuts::decayTopologyMinDca && pair.particle2Dca() > mxeCuts::decayTopologyMinDca &&
       pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
       pair.decayLength() > mxeCuts::decayTopologyMinDca &&
       std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
       ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::decayTopologyMaxD0Dca2Vtx)
   {
     double toFill[] = {pair.pt(), pair.decayLength() * std::sin(pair.pointingAngle()), 
                        pair.particle1Dca(), pair.particle2Dca(), pair.decayLength()};
     mSE_US_DecayTopology->Fill(toFill);
   }

   //Cos theta
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         //std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_US_PointingAngle->Fill(pair.pt(), centrality, std::cos(pair.pointingAngle()));

   //DecayL
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         //pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_US_DecayL->Fill(pair.pt(),  centrality, pair.decayLength());

   //DcaDaughter
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         //pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_US_Dca12->Fill(pair.pt(), centrality, pair.dcaDaughters());

   //PionDca
   if (//pair.particle1Dca() > mxeCuts::pDca[ptIndex] &&
      pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
      pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
      pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
      std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
      ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_US_PionDca2Vtx->Fill(pair.pt(), centrality, pair.particle1Dca());

   //Kaon Dca
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] &&
         //pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_US_KaonDca2Vtx->Fill(pair.pt(), centrality, pair.particle2Dca());

   //D0 dca
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex]
      )
      mSE_US_D0Dca2Vtx->Fill(pair.pt(), centrality, (pair.decayLength()) * sin(pair.pointingAngle()));
}
// --------------------------------------
void StD0Hists::fillSameEvt_LS_QADist(StMixerPair const&  pair, int const centrality)
{
   int ptIndex;
   for (int i = 0; i < mxeCuts::nPtBins; i++)
   {
      if ((pair.pt() >= mxeCuts::PtEdge[i]) && (pair.pt() < mxeCuts::PtEdge[i + 1]))
      {
         ptIndex = i;
         break;
      }
   }
   if (pair.m() <  mxeCuts::QAmassMin || pair.m() > mxeCuts::QAmassMax) return;

   //Cos theta
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         //std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_LS_PointingAngle->Fill(pair.pt(), centrality, std::cos(pair.pointingAngle()));

   //DecayL
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         //pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_LS_DecayL->Fill(pair.pt(), centrality, pair.decayLength());

   //DcaDaughter
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         //pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_LS_Dca12->Fill(pair.pt(), centrality, pair.dcaDaughters());

   //PionDca
   if (//pair.particle1Dca() > mxeCuts::pDca[ptIndex] &&
      pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
      pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
      pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
      std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
      ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_LS_PionDca2Vtx->Fill(pair.pt(), centrality, pair.particle1Dca());

   //Kaon Dca
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] &&
         //pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex] &&
         ((pair.decayLength()) * sin(pair.pointingAngle())) < mxeCuts::dcaV0ToPv[ptIndex])
      mSE_LS_KaonDca2Vtx->Fill(pair.pt(), centrality, pair.particle2Dca());

   //D0 dca
   if (pair.particle1Dca() > mxeCuts::pDca[ptIndex] && pair.particle2Dca() > mxeCuts::kDca[ptIndex] &&
         pair.dcaDaughters() < mxeCuts::dcaDaughters[ptIndex] &&
         pair.decayLength() > mxeCuts::decayLength[ptIndex] &&
         std::cos(pair.pointingAngle()) > mxeCuts::cosTheta[ptIndex]
      )
      mSE_LS_D0Dca2Vtx->Fill(pair.pt(), centrality, (pair.decayLength()) * sin(pair.pointingAngle()));
}
