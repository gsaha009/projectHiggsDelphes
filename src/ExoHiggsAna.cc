#include <iostream>
#include <iomanip>
#include <algorithm>
#include <iterator>
#include <functional>
#include <numeric>
#include <string>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include <fstream>
#include <utility> 
#include <typeinfo>
#include <memory>

#include "TROOT.h"
#include "TSystem.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TRandom.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TH1K.h"

#include "AnaUtil.h"
#include "ZCandidate.h"
#include "ExoHiggsAna.h"


using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::ostringstream;
using std::vector;
using std::map;
using std::pair;
using std::setprecision;
using std::setw;
using std::ios;
using std::setiosflags;
using std::resetiosflags;

using namespace ZSpace;
// -----------
// Constructor
// -----------
ExoHiggsAna::ExoHiggsAna()
{}

// ----------
// Destructor
// ----------
ExoHiggsAna::~ExoHiggsAna() 
{}

// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
//bool ExoHiggsAna::beginJob(const char* inputFile, char* outputFile)  
bool ExoHiggsAna::beginJob(const string& inputFile, string& outputFile)  
{ 
  //string out (outputFile);
  histf = new TFile (outputFile.c_str(), "RECREATE");
  if (!histf) return false;
  histf->cd();
  histf->mkdir("ObjectSelection");
  histf->mkdir("LepPairSelection");
  histf->mkdir("Analysis");

  chain = new TChain("Delphes");
  //string in (inputFile);
  //string out (outputFile);
  string file;
  int line=0;

  std::ifstream myf (inputFile.c_str(), ios::in);
  while (!myf.eof()) {
    line++;
    myf >> file;
    std::cout<<file<<std::endl;
    if (myf) chain->Add(file.c_str());
  }
  treeReader = new ExRootTreeReader(chain);
  if(!treeReader) return false;

  
  bookHistograms();
  eventLoop(treeReader);

  return true;
}

// ---------------
// Book histograms
// ---------------
void ExoHiggsAna::bookHistograms()
{
  histf->cd();
  histf->cd("ObjectSelection");
  //Muon Hists
  new TH1D("MuonCutFlow", "Muon CutFlow", 4, -0.5, 3.5);
  new TH1D("MuonPt", "pT of all Muons",300, 0.0, 300.0);
  new TH1D("MuonEta","eta of all Muons",100, 0.0, 5.0);
  new TH1D("MuonIsoCorr", "Muon Isolation (corr)",100, 0.0, 5.0);
  new TH1D("MuonSumPt_Iso", "Muon SumPt (isolation Variable)",200, 0.0, 200.0);
  new TH1D("MuonSumPtNeutral_Iso", "Muon SumPtNeutral (isolation Variable)",200, 0.0, 200.0);
  new TH1D("MuonSumPtCharged_Iso", "Muon SumPtCharged (isolation Variable)",200, 0.0, 200.0);

  //Electron Hists
  new TH1D("ElectronCutFlow", "Electron CutFlow", 4, -0.5, 3.5);
  new TH1D("ElectronPt", "pT of all Electrons",300, 0.0, 300.0);
  new TH1D("ElectronEta","eta of all Electrons",100, 0.0, 5.0);
  new TH1D("ElectronIsoCorr", "Electron Isolation (corr)",100, 0.0, 5.0);
  new TH1D("ElectronSumPt_Iso", "Electron SumPt (isolation Variable)",200, 0.0, 200.0);
  new TH1D("ElectronSumPtNeutral_Iso", "Electron SumPtNeutral (isolation Variable)",200, 0.0, 200.0);
  new TH1D("ElectronSumPtCharged_Iso", "ELectronon SumPtCharged (isolation Variable)",200, 0.0, 200.0);

  //Jet Hists
  new TH1D ("JetCutFlow", "Jet CutFlow", 4, -0.5, 3.5);
  new TH1D ("JetPt", "pT of all Jets",300, 0.0, 300.0);
  new TH1D ("JetEta","eta of all Jets",100, 0.0, 5.0);
  new TH1D ("JetPhi", "", 100, 0., 4.0);
  new TH1D ("JetDeltaEta", "", 100, 0.0, 2.0);
  new TH1D ("JetDeltaPhi", "", 100, 0.0, 2.0);
  new TH1D ("JetRadius", "", 100, 0.0, 2.0);
  new TH1D ("JetFlavor", "", 10, -0.5, 9.5);
  new TH1D ("JetFlavorAlgo","", 10, -0.5, 9.5);
  new TH1D ("JetFlavorPhys","", 10, -0.5, 9.5);
  new TH1D ("JetBTag", "", 2, -0.5, 1.5);
  new TH1D ("JetBTagAlgo", "", 2, -0.5, 1.5);
  new TH1D ("JetBTagPhys", "", 2, -0.5, 1.5);
  new TH1D ("JetTauTag", "", 2, -0.5, 1.5);
  new TH1D ("JetTauCharge", "TauCharge", 30, -1.5, 1.5);
  new TH1D ("JetEhadOverEem", "", 200, 0.0, 50.0);
  new TH1D ("JetNCharged", "", 50, -0.5, 49.5);
  new TH1D ("JetNNeutrals", "", 50, -0.5, 49.5);
  new TH1D ("JetBeta", "", 100, 0., 2.);
  new TH1D ("JetBetaStar", "", 100, 0., 2.);
  new TH1D ("JetMeanSqDeltaR", "",100, 0., 5.);
  new TH1D ("JetPTD", "",100, 0., 5.0);
  
  new TH1D ("JetNSubJetsTrimmed", "", 50, -0.5, 49.5);
  new TH1D ("JetNSubJetsPruned", "", 50, -0.5, 49.5);
  new TH1D ("JetNSubJetsSoftDropped", "", 50, -0.5, 49.5);
  
  
  histf->cd();
  histf->cd("LepPairSelection");
  new TH1D("dR_l1Zl1", "", 100, 0., 5.0);
  new TH1D("dR_l1Zl2", "", 100, 0., 5.0);
  new TH1D("dR_l1Z", "", 100, 0., 5.0);
  new TH1D("dR_l2Zl1", "", 100, 0., 5.0);
  new TH1D("dR_l2Zl2", "", 100, 0., 5.0);
  new TH1D("dR_l2Z", "", 100, 0., 5.0);
  new TH1D("dR_l1l2", "", 100, 0., 5.0);
  new TH1D("leptonPairCutFlow", "", 6, -0.5, 5.5);
  new TH1D("hCandCutFlow", "", 6, -0.5, 5.5);





  histf->cd();
  histf->cd("Analysis");

  //Event Hists
  new TH1D("evtCutFlow", "Event CutFlow", 8, -0.5, 7.5);
  new TH1D("eventWt", "HepMCEvent Weight",5, -2.5, 2.5);
  new TH1D("nvtx", "Number of Primary vertices", 60, -0.5, 59.5);
  new TH1D("met", "Missing Transver Energy",200, 0, 200);
  new TH1D("nRawEle", "Number of Raw electrons", 10, -0.5, 9.5);
  new TH1D("nRawMuon", "Number of Raw mons", 10, -0.5, 9.5);
  new TH1D("nRawJet", "Number of Raw jets", 10, -0.5, 9.5);
  new TH1D("nGoodMuon", "Number of Good mons", 10, -0.5, 9.5);
  new TH1D("Muon1Pt", "pT of 1st Muon",400, 0.0, 800.0);
  new TH1D("Muon2Pt", "pT of 2nd Muon",400, 0.0, 800.0);
  new TH1D("Muon3Pt", "pT of 3rd Muon",400, 0.0, 800.0);
  new TH1D("Muon4Pt", "pT of 4th Muon",400, 0.0, 800.0);
  new TH1D("nGoodElectron", "Number of Good electrons", 10, -0.5, 9.5);
  new TH1D("Electron1Pt", "pT of 1st Electron",400, 0.0, 800.0);
  new TH1D("Electron2Pt", "pT of 2nd Electron",400, 0.0, 800.0);
  new TH1D("Electron3Pt", "pT of 3rd Electron",400, 0.0, 800.0);
  new TH1D("nGoodJet", "Number of Good jets", 10, -0.5, 9.5);
  new TH1D("Jet1Pt", "pT of 1st Jet",400, 0.0, 800.0);
  new TH1D("Jet2Pt", "pT of 2nd Jet",400, 0.0, 800.0);
  new TH1D("Jet3Pt", "pT of 3rd Jet",400, 0.0, 800.0);
  new TH1D("Jet4Pt", "pT of 4th Jet",400, 0.0, 800.0);
  new TH1D("Jet5Pt", "pT of 5th Jet",400, 0.0, 800.0);
  new TH1D("Jet6Pt", "pT of 6th Jet",400, 0.0, 800.0);
  new TH1D("hLepPt", "LeptonPt", 400, 0.0, 400.0);
  new TH1D("nZCand", "nZ Candidates", 5, -0.5, 4.5);
  new TH1D("massZcand", "ZMass", 500, 0.0, 500.0);
  new TH1D("nGoodLeptons", "", 10, -0.5, 9.5);
  new TH1D("Z2mass", "hMass", 500, 0.0, 500.0);

  histf->cd();
  histf->ls();
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------

void ExoHiggsAna::clearLists() {
  MuonColl.clear();
  ElectronColl.clear();
  JetColl.clear();
}

// -------------------
// The main event loop
// -------------------

void ExoHiggsAna::eventLoop(ExRootTreeReader *treeReader)
{
  /////////////////////////Objects for analysis////////////////////////////

  TClonesArray *BrElec    = treeReader->UseBranch("Electron");
  TClonesArray *BrMuon    = treeReader->UseBranch("Muon");
  TClonesArray *BrJet     = treeReader->UseBranch("Jet");
  TClonesArray *BrMet     = treeReader->UseBranch("MissingET");
  TClonesArray *BrEvent   = treeReader->UseBranch("Event");


  //Long64_t all = treeReader->GetEntries();
  size_t nEntries = treeReader->GetEntries();
  cout << "** Chain contains " << nEntries << " events" << endl;

  ////////////////////********Event Loop*********///////////////////// 
  for (size_t iEntry = 0; iEntry < nEntries; ++iEntry) {
    
    bool verbose {false};
    if ((iEntry/1000) > 0 && iEntry%1000 == 0) verbose = true;
    if (verbose) std::cout<<"Events processed :"<<"\t"<<iEntry<<std::endl;

    clearLists(); // reset analysis related lists for each event
    treeReader->ReadEntry(iEntry);// Load selected branches with data from specified event

    //////////////////EventWeight and EventWeightSum//////////////////
    int evWt = 0;
    HepMCEvent *ev = (HepMCEvent*)BrEvent -> At(0);
    if (ev->Weight < 0.0) evWt = -1;
    else if (ev->Weight > 0.0) evWt = 1;
    else if (ev->Weight == 0.0) evWt = 0;
    AnaUtil::fillHist1D("eventWt", evWt, 1.0);
    /////////////////////////////////////////////////////////////////


    histf->cd();
    histf->cd("Analysis");


    AnaUtil::fillHist1D("evtCutFlow", 0, evWt); //events processed
    int nelec = BrElec->GetEntriesFast();
    AnaUtil::fillHist1D("nRawEle", nelec, evWt);
    int nmuon = BrMuon->GetEntriesFast();
    AnaUtil::fillHist1D("nRawMuon", nmuon, evWt);
    int njet  = BrJet->GetEntriesFast();
    AnaUtil::fillHist1D("nRawJet", njet, evWt);
    MissingET* met_ = (MissingET*) BrMet->At(0);
    AnaUtil::fillHist1D("met", (met_-> MET), evWt);


    histf->cd();
    histf->cd("ObjectSelection");

    //////////////////////////M U O N  S E L E C T I O N////////////////////////////
    for (TObject *_mu: *BrMuon){ 
      const Muon& mu = dynamic_cast<const Muon&> (*_mu);
      AnaUtil::fillHist1D ("MuonCutFlow", 0, evWt);
      AnaUtil::fillHist1D ("MuonPt", mu.PT, evWt);
      AnaUtil::fillHist1D ("MuonEta", std::abs(mu.Eta), evWt);
      AnaUtil::fillHist1D ("MuonPhi", std::abs(mu.Phi), evWt);
      AnaUtil::fillHist1D ("MuonIsoCorr", mu.IsolationVarRhoCorr, evWt);
      AnaUtil::fillHist1D ("MuonSumPt_Iso", mu.SumPt, evWt);
      AnaUtil::fillHist1D ("MuonSumPtNeutral_Iso", mu.SumPtCharged, evWt);
      AnaUtil::fillHist1D ("MuonSumPtCharged_Iso", mu.SumPtNeutral, evWt);

      if (mu.PT < 10) continue;
      AnaUtil::fillHist1D ("MuonCutFlow", 1, evWt);
      if (std::abs(mu.Eta) > 2.4) continue;
      AnaUtil::fillHist1D ("MuonCutFlow", 2, evWt);
      //if (mu.IsolationVarRhoCorr > 0.15) continue;
      AnaUtil::fillHist1D ("MuonCutFlow", 3, evWt);
      TLorentzVector mup4 = mu.P4();
      MuonColl.push_back(mu); 
    }
    std::sort(std::begin(MuonColl), std::end(MuonColl), PtComparator<Muon>()); //sorting Muons


    //////////////E L E C T R O N  S E L E C T I O N////////////////
    for (TObject *_el: *BrElec){ 
      const Electron& ele = dynamic_cast<const Electron&> (*_el);
      AnaUtil::fillHist1D ("ElectronCutFlow", 0, evWt);
      AnaUtil::fillHist1D ("ElectronPt", ele.PT, evWt);
      AnaUtil::fillHist1D ("ElectronEta", std::abs(ele.Eta), evWt);
      AnaUtil::fillHist1D ("ElectronPhi", std::abs(ele.Phi), evWt);
      AnaUtil::fillHist1D ("ElectronIsoCorr", ele.IsolationVarRhoCorr, evWt);
      AnaUtil::fillHist1D ("ElectronSumPt_Iso", ele.SumPt, evWt);
      AnaUtil::fillHist1D ("ElectronSumPtNeutral_Iso", ele.SumPtCharged, evWt);
      AnaUtil::fillHist1D ("ElectronSumPtCharged_Iso", ele.SumPtNeutral, evWt);

      if (ele.PT < 10) continue;
      AnaUtil::fillHist1D ("ElectronCutFlow", 1, evWt);
      if (std::abs(ele.Eta) > 2.5) continue;
      AnaUtil::fillHist1D ("ElectronCutFlow", 2, evWt);
      //if (ele.IsolationVarRhoCorr > 0.1) continue;
      AnaUtil::fillHist1D ("ElectronCutFlow", 3, evWt);
      ElectronColl.push_back(ele); 
    }
    std::sort(std::begin(ElectronColl), std::end(ElectronColl), PtComparator<Electron>()); //sorting Electrons

    ///////////////////J E T  S E L E C T I O N//////////////////////
    for (TObject *_j: *BrJet){ 
      const Jet& jet = dynamic_cast<const Jet&> (*_j);
      AnaUtil::fillHist1D ("JetCutFlow", 0, evWt);
      AnaUtil::fillHist1D ("JetPt", jet.PT, evWt);
      AnaUtil::fillHist1D ("JetEta", std::abs(jet.Eta), evWt);
      AnaUtil::fillHist1D ("JetPhi", std::abs(jet.Phi), evWt);
      AnaUtil::fillHist1D ("JetDeltaEta", jet.DeltaEta, evWt);
      AnaUtil::fillHist1D ("JetDeltaPhi", jet.DeltaPhi, evWt);
      AnaUtil::fillHist1D ("JetRadius", TMath::Sqrt(std::pow(jet.DeltaPhi, 2) + std::pow(jet.DeltaEta, 2)), evWt);
      AnaUtil::fillHist1D ("JetFlavor", jet.Flavor, evWt);
      AnaUtil::fillHist1D ("JetFlavorAlgo", jet.FlavorAlgo, evWt);
      AnaUtil::fillHist1D ("JetFlavorPhys", jet.FlavorPhys, evWt);
      AnaUtil::fillHist1D ("JetBTag", jet.BTag, evWt);
      AnaUtil::fillHist1D ("JetBTagAlgo", jet.BTagAlgo, evWt);
      AnaUtil::fillHist1D ("JetBTagPhys", jet.BTagPhys, evWt);
      AnaUtil::fillHist1D ("JetTauTag", jet.TauTag, evWt);
      AnaUtil::fillHist1D ("JetTauCharge", jet.Charge, evWt);
      AnaUtil::fillHist1D ("JetEhadOverEem", jet.EhadOverEem, evWt);
      AnaUtil::fillHist1D ("JetNCharged", jet.NCharged, evWt);
      AnaUtil::fillHist1D ("JetNNeutrals", jet.NNeutrals, evWt);
      AnaUtil::fillHist1D ("JetBeta", jet.Beta, evWt);
      AnaUtil::fillHist1D ("JetBetaStar", jet.BetaStar, evWt);
      AnaUtil::fillHist1D ("JetMeanSqDeltaR", jet.MeanSqDeltaR, evWt);
      AnaUtil::fillHist1D ("JetPTD", jet.PTD, evWt);

      AnaUtil::fillHist1D ("JetNSubJetsTrimmed", jet.NSubJetsTrimmed, evWt);
      AnaUtil::fillHist1D ("JetNSubJetsPruned", jet.NSubJetsPruned, evWt);
      AnaUtil::fillHist1D ("JetNSubJetsSoftDropped", jet.NSubJetsSoftDropped, evWt);


      if (jet.PT < 30) continue;
      AnaUtil::fillHist1D ("JetCutFlow", 1, evWt);
      if (std::abs(jet.Eta) > 4.7) continue;
      AnaUtil::fillHist1D ("JetCutFlow", 2, evWt);
      //if (jet.IsolationVarRhoCorr > 0.15) continue;
      AnaUtil::fillHist1D ("JetCutFlow", 3, evWt);
      JetColl.push_back(jet); 
    }
    std::sort(std::begin(JetColl), std::end(JetColl), PtComparator<Jet>()); //sorting Jets


    histf->cd();
    histf->cd("Analysis");

    int nGoodMuon = MuonColl.size();
    int nGoodEle  = ElectronColl.size();
    int nGoodJet  = JetColl.size();
    AnaUtil::fillHist1D("nGoodMuon", nGoodMuon, evWt);
    if (MuonColl.size() > 0) AnaUtil::fillHist1D("Muon1Pt", MuonColl[0].PT, evWt);
    if (MuonColl.size() > 1) AnaUtil::fillHist1D("Muon2Pt", MuonColl[1].PT, evWt);
    if (MuonColl.size() > 2) AnaUtil::fillHist1D("Muon3Pt", MuonColl[2].PT, evWt);
    if (MuonColl.size() > 3) AnaUtil::fillHist1D("Muon4Pt", MuonColl[3].PT, evWt);

    AnaUtil::fillHist1D("nGoodElectron", nGoodEle, evWt);
    if (ElectronColl.size() > 0) AnaUtil::fillHist1D("Electron1Pt", ElectronColl[0].PT, evWt);
    if (ElectronColl.size() > 1) AnaUtil::fillHist1D("Electron2Pt", ElectronColl[1].PT, evWt);
    if (ElectronColl.size() > 2) AnaUtil::fillHist1D("Electron3Pt", ElectronColl[2].PT, evWt);

    AnaUtil::fillHist1D("nGoodJet", nGoodJet, evWt);
    if (JetColl.size() > 0) AnaUtil::fillHist1D("Jet1Pt", JetColl[0].PT, evWt);
    if (JetColl.size() > 1) AnaUtil::fillHist1D("Jet2Pt", JetColl[1].PT, evWt);
    if (JetColl.size() > 2) AnaUtil::fillHist1D("Jet3Pt", JetColl[2].PT, evWt);
    if (JetColl.size() > 3) AnaUtil::fillHist1D("Jet4Pt", JetColl[3].PT, evWt);
    if (JetColl.size() > 4) AnaUtil::fillHist1D("Jet5Pt", JetColl[4].PT, evWt);
    if (JetColl.size() > 5) AnaUtil::fillHist1D("Jet6Pt", JetColl[5].PT, evWt);
    


    if (nGoodMuon < 2 && nGoodEle < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 1, evWt);

    // leptons coming from the h1 should have higher pT
    double hLepPt = (nGoodMuon > 0) ? MuonColl[0].PT : 0.0;
    if (nGoodEle > 0 && ElectronColl[0].PT > hLepPt)
      hLepPt = ElectronColl[0].PT;
    AnaUtil::fillHist1D("hLepPt", hLepPt, evWt);

    if (hLepPt < 50.0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2, evWt);


    // Find Z candidates
    std::vector<ZCandidate> ZCandList;
    if (nGoodMuon >= 2) ZSelector(MuonColl, ZCandList);
    if (nGoodEle  >= 2) ZSelector(ElectronColl, ZCandList);
    AnaUtil::fillHist1D("nZcand", ZCandList.size(), evWt);

    if (ZCandList.empty()) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3, evWt);
    
    // Z candidates                                        
    if (ZCandList.size() > 1) 
      std::sort(std::begin(ZCandList), std::end(ZCandList), dmComparator);
    for (const auto& z: ZCandList) 
      AnaUtil::fillHist1D("massZcand", z.mass, evWt);
    
    // The first Z candidate
    const ZCandidate& ZCand = ZCandList[0];
    if (ZCand.mass < 70 || ZCand.mass > 120) continue;
    AnaUtil::fillHist1D("evtCutFlow", 4, evWt);

    // The second Z candidate                                                   
    if (ZCandList.size() > 1) AnaUtil::fillHist1D("Z2mass", ZCandList[1].mass, evWt);

    // Now require at least four tight isolated leptons
    AnaUtil::fillHist1D("nGoodLeptons", nGoodMuon + nGoodEle, evWt);
    if (nGoodMuon + nGoodEle < 4) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5, evWt);

    // Go closer to final topology         
    if (nGoodMuon < 1 || nGoodEle < 1) continue;
    if (nGoodMuon < 3 && nGoodEle < 3) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6, evWt);

    // Now find the ele-mu candidate          
    std::vector<ZCandidate> leptonPairCandList;
    histf->cd();
    histf->cd("LepPairSelection"); //To fill the histograms mentioned in leptonPairSelector function

    auto result = leptonPairSelector(MuonColl, ElectronColl, ZCand, leptonPairCandList, 170.0, evWt, 0.02);
    if (leptonPairCandList.empty()) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7, evWt);

    if (leptonPairCandList.size() > 1)
      std::sort(std::begin(leptonPairCandList), std::end(leptonPairCandList), massComparator);
    const ZCandidate& h1Cand = leptonPairCandList[0];
    AnaUtil::fillHist1D("Z2mass", h1Cand.mass, evWt);


    if (verbose) {
      cout << ">>> "
	   << ", <nMuon>: "     << MuonColl.size()
	   << ", <nElectron>: " << ElectronColl.size()
	   << ", <nJet> "       << JetColl.size()
	   << endl;
    }
  }
  // Analysis over
}

void ExoHiggsAna::endJob() {
  
  histf->cd();
  histf->cd("Analysis");
  vector<string> evLabels {
        "0) Events processed  : ",
	  "1) nEle & nMu >=2  : ",
	  "2) hLepPt > 50     : ",
	  "3) no ZCand        : ",
	  "4) 70<ZCandMass<120: ",
	  "5) nMu+nEle >= 4   : ",
	  "6) sigTopologyCond : ",
	  "7) have hCand      : "
	  
	};
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection");  

#if 0
  double lumiFac = 1.0;
  lumiFac = lumiWt(evtWeightSum_);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWeightSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
       << endl;
  
  
  LL4JMETUtil::scaleHistogram("evtCutFlowWt", lumiFac);
  LL4JMETUtil::showEfficiency("evtCutFlowWt", evLabels, "Event Selection (Weighted)", "Events");
#endif  
}

void ExoHiggsAna::closeHistFiles(){
  histf->cd();
  histf->Write();
  histf->Close();
}
void ExoHiggsAna::closeFiles(){
  closeHistFiles();
  delete treeReader;
  delete chain;
}

#if 0
void ExoHiggsAna::printJob(ostream& os) const
{
  AnaBase::printJob(os);
  os << endl;
  os << "   useEventList: " << std::boolalpha << useEventList_ << endl
     << "  skipDuplicate: " << std::boolalpha << skipDuplicate_ << endl
     << " dumpEventCount: " << dumpEventCount_ << endl
     << "   syncDumpFile: " << dumpFilename_ << endl
     << "   dumpEventMax: " << dumpEventCount_ << endl
     << "  selectPartons: " << std::boolalpha << selectPM_ << endl
     << "     nMEPartons: " << nMEPartons_ << endl;
}
#endif
