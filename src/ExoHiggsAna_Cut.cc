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
{
  fileList_.clear();
  hmap_.clear();
}

// ----------
// Destructor
// ----------
ExoHiggsAna::~ExoHiggsAna() 
{}

// -------------------------------------------------------
// Prepare for the run, do necessary initialisation etc.
// -------------------------------------------------------
//bool ExoHiggsAna::beginJob(const string& inputFile, string& outputFile)  
//bool ExoHiggsAna::beginJob(const std::string& jobFile, std::string& outputFile)  
//bool ExoHiggsAna::beginJob(const std::string& jobFile)  
bool ExoHiggsAna::beginJob()  
{ 
  histf = new TFile (histFile_.c_str(), "RECREATE");
  if (!histf) return false;
  histf->cd();
  histf->mkdir("ObjectSelection");
  histf->mkdir("LepPairSelection");
  histf->mkdir("Analysis");


  /*
  chain = new TChain("Delphes");
  //string in (inputFile);
  //string out (outputFile);
  string file;
  int line=0;

  std::ifstream myf (inputFile.c_str(), ios::in);
  std::cout<<"*****Making Chain*****"<<std::endl;
  while (!myf.eof()) {
    line++;
    myf >> file;
    std::cout<<"Adding File: "<<file<<std::endl;
    if (myf) chain->Add(file.c_str());
  }
  treeReader = new ExRootTreeReader(chain);
  if(!treeReader) return false;

  std::cout<<"*****Booking Histograms*****"<<std::endl;
  bookHistograms();
  std::cout<<"*****Starting Event Loop*****"<<std::endl;
  eventLoop(treeReader);
  */

  std::cout<<"*****Booking Histograms*****"<<std::endl;
  bookHistograms();

  if (_createMVATree) {
    #ifdef TRUE_CPP14
    skimObj_ = std::make_unique <MVASkim> (_mvaInputFile);
    #else
    skimObj_ = std::unique_ptr <MVASkim>(new MVASkim (_mvaInputFile));
    #endif
    if (!skimObj_) return false;
  }

  else if (_readMVA) {
    #ifdef TRUE_CPP14
    _mvaObj = std::make_unique<MVAnalysis>(_MVAnetwork, _MVAxmlFile);
    #else
    _mvaObj = std::unique_ptr <MVAnalysis>(new MVAnalysis (_MVAnetwork, _MVAxmlFile));
    #endif
    if (!_mvaObj) return false;
  }
  

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
  //Muon Histos
  new TH1D("MuonCutFlow", "Muon CutFlow", 4, -0.5, 3.5);
  new TH1D("MuonPt", "pT of all Muons",300, 0.0, 300.0);
  new TH1D("MuonEta","eta of all Muons",100, 0.0, 5.0);
  new TH1D("MuonIsoCorr", "Muon Isolation (corr)",100, 0.0, 5.0);
  new TH1D("MuonSumPt_Iso", "Muon SumPt (isolation Variable)",200, 0.0, 200.0);
  new TH1D("MuonSumPtNeutral_Iso", "Muon SumPtNeutral (isolation Variable)",200, 0.0, 200.0);
  new TH1D("MuonSumPtCharged_Iso", "Muon SumPtCharged (isolation Variable)",200, 0.0, 200.0);

  //Electron Histos
  new TH1D("ElectronCutFlow", "Electron CutFlow", 4, -0.5, 3.5);
  new TH1D("ElectronPt", "pT of all Electrons",300, 0.0, 300.0);
  new TH1D("ElectronEta","eta of all Electrons",100, 0.0, 5.0);
  new TH1D("ElectronIsoCorr", "Electron Isolation (corr)",100, 0.0, 5.0);
  new TH1D("ElectronSumPt_Iso", "Electron SumPt (isolation Variable)",200, 0.0, 200.0);
  new TH1D("ElectronSumPtNeutral_Iso", "Electron SumPtNeutral (isolation Variable)",200, 0.0, 200.0);
  new TH1D("ElectronSumPtCharged_Iso", "ELectronon SumPtCharged (isolation Variable)",200, 0.0, 200.0);

  //Tau Histos
  new TH1D ("TauCutFlow", "Tau CutFlow", 4, -0.5, 3.5);
  new TH1D ("TauPt", "pT of all Taus",300, 0.0, 300.0);
  new TH1D ("TauEta","eta of all Taus",100, -5.0, 5.0);
  new TH1D ("TauPhi", "", 100, 0., 4.0);

  
  //Jet Histos
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

  new TH1D ("jetMuDR", "", 100, 0.0, 5.);
  new TH1D ("jetEleDR", "", 100, 0.0, 5.);
  
  
  histf->cd();
  histf->cd("LepPairSelection");
  new TH1D("dR_l1Zl1", "", 100, 0., 5.0);
  new TH1D("dR_l1Zl2", "", 100, 0., 5.0);
  new TH1D("dR_l1Z", "", 100, 0., 5.0);
  new TH1D("lepCh", "", 3, -1.5, 1.5);
  new TH1D("tauCh", "", 3, -1.5, 1.5);
  new TH1D("sumCh", "", 5, -2.5, 2.5);
  new TH1D("dR_l2Zl1", "", 100, 0., 5.0);
  new TH1D("dR_l2Zl2", "", 100, 0., 5.0);
  new TH1D("dR_l2Z", "", 100, 0., 5.0);
  new TH1D("dR_l1l2", "", 100, 0., 5.0);
  new TH1D("leptonPairCutFlow", "", 6, -0.5, 5.5);
  new TH1D("hCandCutFlow", "", 6, -0.5, 5.5);
  new TH1D("SumCharge", "", 4, -0.5, 3.5);
  new TH1D("lepDR", "", 100, 0.0, 5.0);





  histf->cd();
  histf->cd("Analysis");

  //Event Hists
  new TH1D("evtCutFlow", "Event CutFlow", 12, -0.5, 11.5);
  //  new TH1D("eventWt", "HepMCEvent Weight",100000000000, -5000000000., 5000000000.0);
  new TH1D("eventWt_", "HepMCEvent Weight",3, -1.5, 1.5);
  new TH1D("nvtx", "Number of Primary vertices", 60, -0.5, 59.5);
  new TH1D("met", "Missing Transver Energy",200, 0, 200);
  new TH1D("nRawEle", "Number of Raw electrons", 10, -0.5, 9.5);
  new TH1D("nRawMuon", "Number of Raw mons", 10, -0.5, 9.5);
  new TH1D("nRawJet", "Number of Raw jets", 10, -0.5, 9.5);
  new TH1D("nRawTau", "Number of Raw taus", 10, -0.5, 9.5);
  new TH1D("nGoodMuon", "Number of Good mons", 10, -0.5, 9.5);
  new TH1D("Muon1Pt", "pT of 1st Muon",400, 0.0, 800.0);
  new TH1D("Muon2Pt", "pT of 2nd Muon",400, 0.0, 800.0);
  new TH1D("Muon3Pt", "pT of 3rd Muon",400, 0.0, 800.0);
  new TH1D("Muon4Pt", "pT of 4th Muon",400, 0.0, 800.0);
  new TH1D("nGoodElectron", "Number of Good electrons", 10, -0.5, 9.5);
  new TH1D("Electron1Pt", "pT of 1st Electron",400, 0.0, 800.0);
  new TH1D("Electron2Pt", "pT of 2nd Electron",400, 0.0, 800.0);
  new TH1D("Electron3Pt", "pT of 3rd Electron",400, 0.0, 800.0);
  new TH1D("nGoodTau", "Number of Good taus", 10, -0.5, 9.5);
  new TH1D("Tau1Pt", "pT of 1st Tau",400, 0.0, 800.0);
  new TH1D("Tau2Pt", "pT of 2nd Tau",400, 0.0, 800.0);
  new TH1D("Tau3Pt", "pT of 3rd Tau",400, 0.0, 800.0);
  new TH1D("nGoodJet", "Number of Good jets", 10, -0.5, 9.5);
  new TH1D("Jet1Pt", "pT of 1st Jet",400, 0.0, 800.0);
  new TH1D("Jet2Pt", "pT of 2nd Jet",400, 0.0, 800.0);
  new TH1D("Jet3Pt", "pT of 3rd Jet",400, 0.0, 800.0);
  new TH1D("Jet4Pt", "pT of 4th Jet",400, 0.0, 800.0);
  new TH1D("Jet5Pt", "pT of 5th Jet",400, 0.0, 800.0);
  new TH1D("Jet6Pt", "pT of 6th Jet",400, 0.0, 800.0);

  new TH1D("nNonTauJet",  "nNonTauJet",10, -0.5, 9.5);
  new TH1D("nonTauJetPt",  "pT of all qcd Jets",400, 0.0, 800.0);
  new TH1D("nonTauJet1Pt", "pT of 1st Jet",400, 0.0, 800.0);
  new TH1D("nonTauJet2Pt", "pT of 2nd Jet",400, 0.0, 800.0);
  new TH1D("nonTauJet3Pt", "pT of 3rd Jet",400, 0.0, 800.0);
  new TH1D("nonTauJet4Pt", "pT of 4th Jet",400, 0.0, 800.0);

  new TH1D("hMass_almost", "hMass", 500, 0.0, 1000.0);
  new TH1D("muTauDR", "DeltaR between mu and tau", 100, 0.0, 5.0);
  new TH1D("hLepPt", "LeptonPt", 800, 0.0, 800.0);
  new TH1D("hMuPt", "Leading Muon Pt", 800, 0.0, 800.0);
  new TH1D("nZCand", "nZ Candidates", 5, -0.5, 4.5);
  new TH1D("Z1CandMass", "ZMass", 500, 0.0, 500.0);
  new TH1D("Z2CandMass", "Z2Mass", 500, 0.0, 500.0);
  new TH1D("nGoodLeptons", "", 10, -0.5, 9.5);
  new TH1D("Tau1pT_5", "Tau1Pt", 800, 0.0, 800.0);
  new TH1D("Tau2pT_5", "Tau2Pt", 800, 0.0, 800.0);

  new TH1D("ZmT", "Z Transverse Mass", 1000, 0.0, 1000.0);
  new TH1D("ZmT_manual", "Z Transverse Mass", 1000, 0.0, 1000.0);
  new TH1D("h2MT_1", "h Transverse Mass", 500, 0.0, 1000.0);
  new TH1D("h2MToverPt_1", "h Transverse Mass/Lep pT", 100, 0.0, 10.0);
  new TH1D("h2MT_2", "h Transverse Mass", 1000, 0.0, 1000.0);
  new TH1D("h2MToverPt_2", "h Transverse Mass/Lep pT", 100, 0.0, 10.0);
  new TH1D("h2MT_3", "h Transverse Mass", 500, 0.0, 1000.0);
  new TH1D("h2MToverPt_3", "h Transverse Mass/Lep pT", 100, 0.0, 10.0);
  new TH1D("h1_l1l2DR", "DeltaR between mu and tau from h2", 100, 0.0, 5.0);
  new TH1D("h1_l1l2DPhi", "DeltaPhi between mu and tau from h2", 100, -5.0, 5.0);
  new TH1D("h1_l1l2DEta", "DeltaEta between mu and tau from h2", 100, -5.0, 5.0);

  new TH1D("h1l1M", "", 1000, 0.0, 5.0);
  new TH1D("h1l2M", "", 1000, 0.0, 5.0);
  new TH1D("ZhDR", "DeltaR between Z & h2", 100, 0.0, 5.0);
  new TH1D("ZhPhi", "DeltaPhi between Z & h2", 100, -5.0, 5.0);
  new TH1D("ZhDEta", "DeltaEta between Z & h2", 100, -5.0, 5.0);

  new TH1D("h1Mass", "hMass", 500, 0.0, 500.0);
  new TH1D("HT", "Jet Scalar Sum pT", 400, 0.0, 600.0);
  new TH1D("LT", "Lepton Scalar Sum pT", 500, 0.0, 1500.0);
  new TH1D("ST", "All Scalar Sum pT", 500, 0.0, 1500.0);
  new TH1D("Jet1Pt_stage2", "pT of 1st Jet",400, 0.0, 800.0);
  new TH1D("met_PreCut", "MET", 500, 0.0, 1500.0);
  new TH1D("LT_PreCut", "Lepton Scalar Sum pT", 500, 0.0, 1500.0);
  new TH1D("nbJets", "Number of b Jets", 10, -0.5, 9.5);
  //new TH1D("ZhDEta", "DeltaEta between Z & h2", 100, -5.0, 5.0);
  //new TH1D("ZhDPhi", "DeltaPhi between Z & h2", 100, -5.0, 5.0);
  //new TH1D("ZhDR", "DeltaR between Z & h2", 100, 0.0, 5.0);

  new TH1D("ZMuDEta", "DeltaEta between Z & Mu from h2", 100, -5.0, 5.0);
  new TH1D("ZMuDPhi", "DeltaPhi between Z & Mu from h2", 100, -5.0, 5.0);
  new TH1D("ZMuDR", "DeltaR between Z & Mu from h2", 100, 0.0, 5.0);

  new TH1D("nTau", "Number of Taus", 10, -0.5, 9.5);
  new TH1D("nGenTau", "Number of Generator Taus", 10, -0.5, 9.5);
  new TH1D("genTauPt", "", 300, 0.,300.);
  new TH1D("genTauEta", "", 100, -5., 5.);
  new TH1D("genMuPt", "", 300, 0.,300.);
  new TH1D("genMuEta", "", 100, -5., 5.);
  new TH1D("genMuTauDEta", "", 100, -5.0, 5.0);
  new TH1D("genMuTauDPhi", "", 100, -5.0, 5.0);
  new TH1D("genMuTauDR", "DeltaR between mu and tau GEN", 100, 0.0, 5.0);

  //for ZtoTauTau_01Jet
  new TH1D("JetTauId_Z", "", 2, -0.5, 1.5);
  new TH1D("TauPt_Z", "", 300, 0.,300.);
  new TH1D("TauEta_Z", "", 100, -5., 5.);
  new TH1D("TauPhi_Z", "", 100, -5., 5.);
  new TH1D("nRawTau_Z", "", 6, -0.5, 5.5);


  new TH1D("METovHT_1", "", 500, 0.0, 5.0);
  new TH1D("METovHT_2", "", 500, 0.0, 5.0);
  new TH1D("METovHT_noTa", "", 500, 0.0, 5.0);
  new TH1D("METovHT_Ta", "", 500, 0.0, 5.0);
  new TH1D("LT_noTa", "", 1000, 0.0, 1000.0);
  new TH1D("LT_Ta", "", 1000, 0.0, 1000.0);
  new TH1D("h2M_1_corr", "h Transverse Mass", 500, 0.0, 1000.0);
  new TH1D("h2M_3_corr", "h Transverse Mass", 500, 0.0, 1000.0);
  new TH1D("test", "", 500, -25.0, 25.0);
  new TH1D("muCharge", "", 100, -5.0, 5.0);
  
  new TH1D("ZMuDEta+", "DeltaEta between Z & Mu from h2", 100, -5.0, 5.0);
  new TH1D("ZMuDPhi+", "DeltaPhi between Z & Mu from h2", 100, -5.0, 5.0);
  new TH1D("ZMuDR+", "DeltaR between Z & Mu from h2", 100, 0.0, 5.0);
  new TH1D("ZTauDEta+", "DeltaEta between Z & Tau from h2", 100, -5.0, 5.0);
  new TH1D("ZTauDPhi+", "DeltaPhi between Z & Tau from h2", 100, -5.0, 5.0);
  new TH1D("ZTauDR+", "DeltaR between Z & Tau from h2", 100, 0.0, 5.0);

  new TH1D("ZMuDEta-", "DeltaEta between Z & Mu from h2", 100, -5.0, 5.0);
  new TH1D("ZMuDPhi-", "DeltaPhi between Z & Mu from h2", 100, -5.0, 5.0);
  new TH1D("ZMuDR-", "DeltaR between Z & Mu from h2", 100, 0.0, 5.0);
  new TH1D("ZTauDEta-", "DeltaEta between Z & Tau from h2", 100, -5.0, 5.0);
  new TH1D("ZTauDPhi-", "DeltaPhi between Z & Tau from h2", 100, -5.0, 5.0);
  new TH1D("ZTauDR-", "DeltaR between Z & Tau from h2", 100, 0.0, 5.0);

  new TH1D("dEtaP", "DeltaEta", 100, 0.0, 10.0);
  new TH1D("dPhiP", "DeltaPhi", 100, 0.0, 10.0);
  new TH1D("dRP", "DeltaR", 100, 0.0, 10.0);
  new TH1D("dEtaM", "DeltaEta", 100, 0.0, 10.0);
  new TH1D("dPhiM", "DeltaPhi", 100, 0.0, 10.0);
  new TH1D("dRM", "DeltaR", 100, 0.0, 10.0);

  
  //Signal

  
  histf->cd();
  histf->ls();
}

// -------------------------------
// Clear vectors before event loop
// -------------------------------

void ExoHiggsAna::clearLists() {
  MuonColl.clear();
  ElectronColl.clear();
  TauColl.clear();
  JetColl.clear();
}

// -------------------
// The main event loop
// -------------------

void ExoHiggsAna::eventLoop(ExRootTreeReader *treeReader)
{
  /////////////////////////Objects for analysis////////////////////////////
  TClonesArray *BrGen     = treeReader->UseBranch("Particle");
  TClonesArray *BrElec    = treeReader->UseBranch("Electron");
  TClonesArray *BrMuon    = treeReader->UseBranch("Muon");
  TClonesArray *BrJet     = treeReader->UseBranch("Jet");
  TClonesArray *BrMet     = treeReader->UseBranch("MissingET");
  TClonesArray *BrEvent   = treeReader->UseBranch("Event");
  TClonesArray *BrWeight  = treeReader->UseBranch("Weight");

  size_t nEntries = treeReader->GetEntries();
  cout << "** Chain contains " << nEntries << " events" << endl;

  size_t nEvents = (maxEvt_ < 0) ? nEntries : maxEvt_;
  cout << "** Starting Analysis with  " << nEvents << " events" << endl;

  ////////////////////////////////////********Event Loop*********//////////////////////////////////// 
  for (size_t iEntry = 0; iEntry < nEvents; ++iEntry) {
    bool verbose {false};
    if ((iEntry/10000) > 0 && iEntry%10000 == 0) verbose = true;
    if (verbose) std::cout<<"Events processed :"<<"\t"<<iEntry<<std::endl;

    clearLists(); // reset analysis related lists for each event
    treeReader->ReadEntry(iEntry);// Load selected branches with data from specified event

    //////////////////EventWeight and EventWeightSum//////////////////
    int evWt = 1;
    
    HepMCEvent *ev = (HepMCEvent*)BrEvent -> At(0);
    
    //    AnaUtil::fillHist1D("eventWt_raw", ev->Weight, 1.0);
    //   std::cout<<BrWeight->GetEntriesFast()<<"\t"<<ev->Weight<<std::endl;
    //    std::cout<<iEntry<<"\t"<<((Weight*)BrWeight->At(0))->Weight<<"\t"<<((Weight*)BrWeight->At(0))->Weight<<"\t"<<((Weight*)BrWeight->At(0))->Weight<<"\t"<<((Weight*)BrWeight->At(0))->Weight<<"\t"<<((Weight*)BrWeight->At(0))->Weight<<"\t"<<((Weight*)BrWeight->At(0))->Weight<<std::endl;
    //double weight = ((Weight*)BrWeight->At(0))->Weight;
    //AnaUtil::fillHist1D("eventWt", weight, 1.0);
    //evWt = (weight > 0.0) ? 1 : -1;
    AnaUtil::fillHist1D("eventWt_", evWt, 1.0);
    evtWtSum_ += evWt;
    /////////////////////////////////////////////////////////////////


    histf->cd();
    histf->cd("Analysis");

    if (readGenInfo()) {
      //if (iEntry < 10) dumpGenInfo(BrGen, std::cout);
      int iGenTau = 0;
      bool hasTauFromH2 {false};
      bool hasMuFromH2  {false};
      TLorentzVector gentaup4, genmup4;
      for (TObject *g: *BrGen){
	const GenParticle& genp = dynamic_cast<const GenParticle&> (*g);
	int pdgid = std::fabs(genp.PID);
	int status = genp.Status;
	if (pdgid == 15 && (status == 23 || status == 2)) {
	  int momid = genp.M1;
	  if (momid < 0) continue;
	  const GenParticle& momgenp = dynamic_cast<const GenParticle&> (*BrGen->At(momid));
	  if (std::fabs(momgenp.PID) != 25000000) continue;
	  iGenTau++;
	  gentaup4 = genp.P4();
	  AnaUtil::fillHist1D("genTauPt", gentaup4.Pt(), 1.0);
	  AnaUtil::fillHist1D("genTauEta", gentaup4.Eta(), 1.0);
	  hasTauFromH2 = true;
	}
	if (pdgid == 13 && status == 1) {
	  int mmid = -1;
	  int momid = getMotherId(genp, BrGen, mmid);
	  if (momid < 0) continue;
	  const GenParticle& momgenp = dynamic_cast<const GenParticle&> (*BrGen->At(momid));
	  if (std::fabs(momgenp.PID) != 25000000) continue;
	  genmup4 = genp.P4();
	  AnaUtil::fillHist1D("genMuPt", genmup4.Pt(), 1.0);
	  AnaUtil::fillHist1D("genMuEta", genmup4.Eta(), 1.0);
	  hasMuFromH2 = true;
	}
      }
      AnaUtil::fillHist1D ("nGenTau", iGenTau, 1.0);
      if (hasTauFromH2 && hasMuFromH2) {
	AnaUtil::fillHist1D ("genMuTauDEta", (genmup4.Eta() - gentaup4.Eta()), 1.0);
	AnaUtil::fillHist1D ("genMuTauDPhi", (genmup4.Phi() - gentaup4.Phi()), 1.0);
	AnaUtil::fillHist1D ("genMuTauDR", genmup4.DeltaR(gentaup4), 1.0);
      }
    }

    AnaUtil::fillHist1D("evtCutFlow", 0, evWt); //events processed
    int nelec = BrElec->GetEntriesFast();
    AnaUtil::fillHist1D("nRawEle", nelec, evWt);
    int nmuon = BrMuon->GetEntriesFast();
    AnaUtil::fillHist1D("nRawMuon", nmuon, evWt);
    int njet  = BrJet->GetEntriesFast();
    AnaUtil::fillHist1D("nRawJet", njet, evWt);
    //    int ntau = BrTau->GetEntries();
    //    AnaUtil::fillHist1D("nRawTau", ntau, evWt);
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

      if (mu.PT <= AnaUtil::cutValue(muonCutMap(), "pt")) continue;
      AnaUtil::fillHist1D ("MuonCutFlow", 1, evWt);
      if (std::abs(mu.Eta) >= AnaUtil::cutValue(muonCutMap(), "eta")) continue;
      AnaUtil::fillHist1D ("MuonCutFlow", 2, evWt);
      //if (mu.IsolationVarRhoCorr > AnaUtil::cutValue(muonCutMap(), "iso")) continue;
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

      if (ele.PT <= AnaUtil::cutValue(electronCutMap(), "pt")) continue;
      AnaUtil::fillHist1D ("ElectronCutFlow", 1, evWt);
      if (std::abs(ele.Eta) >= AnaUtil::cutValue(electronCutMap(), "eta")) continue;
      AnaUtil::fillHist1D ("ElectronCutFlow", 2, evWt);
      //if (ele.IsolationVarRhoCorr > AnaUtil::cutValue(electronCutMap(), "iso")) continue;
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


      if (jet.PT <= AnaUtil::cutValue(jetCutMap(), "pt")) continue;
      AnaUtil::fillHist1D ("JetCutFlow", 1, evWt);
      if (std::abs(jet.Eta) >= AnaUtil::cutValue(jetCutMap(), "eta")) continue;
      AnaUtil::fillHist1D ("JetCutFlow", 2, evWt);
      if (!(jetLeptonCleaning(jet, MuonColl, ElectronColl, AnaUtil::cutValue(jetCutMap(), "mindRlep")))) continue;
      AnaUtil::fillHist1D ("JetCutFlow", 3, evWt);
      JetColl.push_back(jet); 
    }
    std::sort(std::begin(JetColl), std::end(JetColl), PtComparator<Jet>()); //sorting Jets

    
    int iTau = 0;
    for (TObject *_tau: *BrJet){
      const Jet& tauj = dynamic_cast<const Jet&> (*_tau);
      AnaUtil::fillHist1D ("JetTauId", tauj.TauTag, evWt);
      if (tauj.TauTag == 0) continue;
      iTau++;
      AnaUtil::fillHist1D ("TauCutFlow", 0, evWt);
      AnaUtil::fillHist1D ("TauPt", tauj.PT, evWt);
      AnaUtil::fillHist1D ("TauEta", tauj.Eta, evWt);
      AnaUtil::fillHist1D ("TauPhi", std::abs(tauj.Phi), evWt);

      if (tauj.PT <= AnaUtil::cutValue(tauCutMap(), "pt")) continue;
      AnaUtil::fillHist1D ("TauCutFlow", 1, evWt);
      if (std::abs(tauj.Eta) >= AnaUtil::cutValue(tauCutMap(), "eta")) continue;
      AnaUtil::fillHist1D ("TauCutFlow", 2, evWt);
      //if (ele.IsolationVarRhoCorr > AnaUtil::cutValue(electronCutMap(), "iso")) continue;
      AnaUtil::fillHist1D ("TauCutFlow", 3, evWt);
      TauColl.push_back(tauj); 
    }
    std::sort(std::begin(TauColl), std::end(TauColl), PtComparator<Jet>()); //sorting Taus
    


    histf->cd();
    histf->cd("Analysis");

    AnaUtil::fillHist1D("nRawTau", iTau, evWt);
    
    int ijet = 0;
    for (auto& j: JetColl) {
      if (j.TauTag == 1) continue;
      ijet++;
      AnaUtil::fillHist1D ("nonTauJetPt", j.PT, evWt);
      if (ijet == 1) 	  AnaUtil::fillHist1D ("nonTauJet1Pt", j.PT, evWt);
      else if (ijet == 2) AnaUtil::fillHist1D ("nonTauJet2Pt", j.PT, evWt);
      else if (ijet == 3) AnaUtil::fillHist1D ("nonTauJet3Pt", j.PT, evWt);
      else if (ijet == 4) AnaUtil::fillHist1D ("nonTauJet4Pt", j.PT, evWt);
    }
    AnaUtil::fillHist1D("nNonTauJet", ijet, evWt);
    
    int nGoodMuon = MuonColl.size();
    int nGoodEle  = ElectronColl.size();
    int nGoodTau  = TauColl.size();
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

    AnaUtil::fillHist1D("nGoodLepton", nGoodMuon+nGoodEle, evWt);
    
    AnaUtil::fillHist1D("nGoodTau", nGoodTau, evWt);
    if (TauColl.size() > 0) AnaUtil::fillHist1D("Tau1Pt", TauColl[0].PT, evWt);
    if (TauColl.size() > 1) AnaUtil::fillHist1D("Tau2Pt", TauColl[1].PT, evWt);
    if (TauColl.size() > 2) AnaUtil::fillHist1D("Tau3Pt", TauColl[2].PT, evWt);
	
    
    AnaUtil::fillHist1D("nGoodJet", nGoodJet, evWt);
    if (JetColl.size() > 0) AnaUtil::fillHist1D("Jet1Pt", JetColl[0].PT, evWt);
    if (JetColl.size() > 1) AnaUtil::fillHist1D("Jet2Pt", JetColl[1].PT, evWt);
    if (JetColl.size() > 2) AnaUtil::fillHist1D("Jet3Pt", JetColl[2].PT, evWt);
    if (JetColl.size() > 3) AnaUtil::fillHist1D("Jet4Pt", JetColl[3].PT, evWt);
    if (JetColl.size() > 4) AnaUtil::fillHist1D("Jet5Pt", JetColl[4].PT, evWt);
    if (JetColl.size() > 5) AnaUtil::fillHist1D("Jet6Pt", JetColl[5].PT, evWt);
    
    /*
    //To check tau reconstruction efficiency from ZtoTauTau_01Jet sample
    std::vector<Jet> TauBox;
    int itau = 0;
    for (TObject *_tau: *BrJet){
      const Jet& tauj = dynamic_cast<const Jet&> (*_tau);
      AnaUtil::fillHist1D ("JetTauId_Z", tauj.TauTag, evWt);
      if (tauj.TauTag == 0) continue;
      itau++;
      AnaUtil::fillHist1D ("TauPt_Z", tauj.PT, evWt);
      AnaUtil::fillHist1D ("TauEta_Z", tauj.Eta, evWt);
      AnaUtil::fillHist1D ("TauPhi_Z", std::abs(tauj.Phi), evWt);
      TauBox.push_back(tauj); 
    }
    AnaUtil::fillHist1D("nRawTau_Z", itau, evWt);
    */

    if (nGoodMuon > 0) AnaUtil::fillHist1D ("mu1Eta", MuonColl[0].Eta, evWt);
    
    if (nGoodMuon < 2 && nGoodEle < 2) continue;
    AnaUtil::fillHist1D("evtCutFlow", 1, evWt);
    
    
    // leptons coming from the h1 should have higher pT
    double hLepPt = (nGoodMuon > 0) ? MuonColl[0].PT : 0.0;
    if (nGoodEle > 0 && ElectronColl[0].PT > hLepPt)
      hLepPt = ElectronColl[0].PT;
    AnaUtil::fillHist1D("hLepPt", hLepPt, evWt);

    //    if (hLepPt < AnaUtil::cutValue(evselCutMap(), "hLepPtMin")) continue;
    if (hLepPt < 70) continue;
    AnaUtil::fillHist1D("evtCutFlow", 2, evWt);


    // Find Z candidates
    std::vector<ZCandidate> ZCandList;
    histf->cd();
    histf->cd("LepPairSelection"); //To fill the histograms mentioned in leptonPairSelector function
    if (nGoodMuon >= 2) ZSelector(MuonColl, ZCandList);
    if (nGoodEle  >= 2) ZSelector(ElectronColl, ZCandList);
    histf->cd();
    histf->cd("Analysis");
    AnaUtil::fillHist1D("nZcand", ZCandList.size(), evWt);

    if (ZCandList.empty()) continue;
    AnaUtil::fillHist1D("evtCutFlow", 3, evWt);
    
    // Z candidates                                        
    if (ZCandList.size() > 1) 
      std::sort(std::begin(ZCandList), std::end(ZCandList), dmComparator);
    
    // The first Z candidate
    const ZCandidate& ZCand = ZCandList[0];
    AnaUtil::fillHist1D ("Z1CandMass", ZCand.mass, evWt);
    if (ZCand.mass < AnaUtil::cutValue(evselCutMap(), "ZMassLow") ||
        ZCand.mass > AnaUtil::cutValue(evselCutMap(), "ZMassHigh")) continue;
    AnaUtil::fillHist1D("evtCutFlow", 4, evWt);

    TLorentzVector _l1p4;
    _l1p4.SetPtEtaPhiE(ZCand.l1P4.Pt(), 0.0, ZCand.l1P4.Phi(), ZCand.l1P4.E());
    TLorentzVector _l2p4;
    _l2p4.SetPtEtaPhiE(ZCand.l2P4.Pt(), 0.0, ZCand.l2P4.Phi(), ZCand.l2P4.E());
    AnaUtil::fillHist1D ("ZmT", (ZCand.l1P4 + ZCand.l2P4).Mt(), evWt);
    AnaUtil::fillHist1D ("ZmT_manual", (_l1p4 + _l2p4).M(), evWt);
    
    // The second Z candidate                                                   
    if (ZCandList.size() > 1) AnaUtil::fillHist1D("Z2CandMass", ZCandList[1].mass, evWt);

    // Now require at least four tight isolated leptons
    AnaUtil::fillHist1D("nGoodLeptons", nGoodMuon + nGoodEle, evWt);
    if (nGoodMuon + nGoodEle < 3) continue;
    AnaUtil::fillHist1D("evtCutFlow", 5, evWt);

    if (nGoodMuon < 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 6, evWt);
    
    if (nGoodTau == 1) AnaUtil::fillHist1D ("Tau1pT_5", TauColl[0].PT, evWt);
    if (nGoodTau > 1)  AnaUtil::fillHist1D ("Tau2pT_5", TauColl[1].PT, evWt);

    AnaUtil::fillHist1D ("hMuPt", MuonColl[0].PT, evWt);
    //    if (MuonColl[0].PT < 70) continue;
    //    AnaUtil::fillHist1D("evtCutFlow", 6, evWt);

    AnaUtil::fillHist1D ("METovHT_1", (met_->MET)/(MuonColl[0].PT + ZCand.l1P4.Pt() + ZCand.l2P4.Pt()), evWt);

    //Tight requirement: nTau == 1
    if (nGoodTau != 1) continue;
    AnaUtil::fillHist1D("evtCutFlow", 7, evWt);

    AnaUtil::fillHist1D ("METovHT_2", (met_->MET)/(MuonColl[0].PT + ZCand.l1P4.Pt() + ZCand.l2P4.Pt() + TauColl[0].PT), evWt);
	
    AnaUtil::fillHist1D ("hMass_almost", (MuonColl[0].P4()+TauColl[0].P4()).M(), evWt);
    AnaUtil::fillHist1D ("muTauDR", MuonColl[0].P4().DeltaR(TauColl[0].P4()), evWt);
    
    TLorentzVector tap4 = TauColl[0].P4();
    TLorentzVector mup4 = MuonColl[0].P4();
    TLorentzVector metp4;
    metp4.SetPtEtaPhiE(met_->MET, 0.0, met_->Phi, met_->MET);
    double h2MT_1 = Calculate_TotalMT(mup4, tap4, metp4);

    AnaUtil::fillHist1D ("h2MT_1", h2MT_1, evWt);
    AnaUtil::fillHist1D ("h2MToverPt_1", h2MT_1/(tap4.Pt()+mup4.Pt()), evWt);

    double corr = 1.0/(1 + metp4.Px()/tap4.Px());
    //    if (corr > 0.0 && corr <= 10.00) std::cout<<"corrGood: "<<corr<<"  "<<metp4.Px()<<"  "<<tap4.Px()<<std::endl;
    //    else std::cout<<"corrBad: "<<corr<<"  "<<metp4.Px()<<"  "<<tap4.Px()<<std::endl;
    TLorentzVector Ntap4;
    Ntap4.SetPtEtaPhiE((1.0/corr)*TauColl[0].P4().Pt(), TauColl[0].P4().Eta(), TauColl[0].P4().Phi(), TauColl[0].P4().E());
    double h2MT_1_corr = Calculate_TotalMT(mup4, Ntap4, metp4);
    AnaUtil::fillHist1D ("h2M_1_corr", h2MT_1_corr, evWt);

    
    // Now find the ele-mu candidate          
    std::vector<ZCandidate> leptonPairCandList;
    histf->cd();
    histf->cd("LepPairSelection"); //To fill the histograms mentioned in leptonPairSelector function
    auto result = leptonPairSelector(MuonColl, TauColl, ZCand, leptonPairCandList, 170.0, evWt, AnaUtil::cutValue(evselCutMap(), "minDRLP"));
    histf->cd();
    histf->cd("Analysis");
    if (leptonPairCandList.empty()) continue;
    AnaUtil::fillHist1D("evtCutFlow", 8, evWt);

    if (leptonPairCandList.size() > 1)
      std::sort(std::begin(leptonPairCandList), std::end(leptonPairCandList), massComparator);
    const ZCandidate& h1Cand = leptonPairCandList[0];
    AnaUtil::fillHist1D("h1Mass", h1Cand.mass, evWt);

    AnaUtil::fillHist1D ("ZhDEta", (ZCand.p4.Eta()-h1Cand.p4.Eta()), evWt);
    AnaUtil::fillHist1D ("ZhDPhi", TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1Cand.p4.Phi()), evWt);
    AnaUtil::fillHist1D ("ZhDR", (ZCand.p4).DeltaR(h1Cand.p4), evWt);

    /*    TLorentzVector _lp4;
    _lp4.SetPtEtaPhiE((h1Cand.p4).Pt(), 0.0, (h1Cand.p4).Phi(), (h1Cand.p4).E());
    TLorentzVector _metp4;
    _metp4.SetPtEtaPhiE(met_->MET, 0.0, met_->Phi, met_->MET);
    TLorentzVector _hp4= (_lp4 + _metp4);
    AnaUtil::fillHist1D ("_hmT", _hp4.M(), evWt);
    AnaUtil::fillHist1D ("_hmToverPt", _hp4.M()/(h1Cand.p4).Pt(), evWt);
    */


    TLorentzVector h1l1p4 = h1Cand.l1P4;
    TLorentzVector h1l2p4 = h1Cand.l2P4;
    TLorentzVector _metp4;
    _metp4.SetPtEtaPhiE(met_->MET, 0.0, met_->Phi, met_->MET);

    double h2MT_2 = Calculate_TotalMT_h2(h1l1p4, h1l2p4, _metp4);
    
    AnaUtil::fillHist1D ("h2MT_2", h2MT_2, evWt);
    AnaUtil::fillHist1D ("h2MToverPt_2", h2MT_2/(h1l1p4.Pt()+h1l2p4.Pt()), evWt);

    //    if (h1l2p4.Pt() < 70) continue;
    double h2MT_3 = Calculate_TotalMT(h1l1p4, h1l2p4, _metp4);
    
    AnaUtil::fillHist1D ("h2MT_3", h2MT_3, evWt);
    AnaUtil::fillHist1D ("h2MToverPt_3", h2MT_3/(h1l1p4.Pt()+h1l2p4.Pt()), evWt);

    double corr_3 = 1.0/(1 + metp4.Px()/h1l2p4.Px());
    //    if (corr > 0.0 && corr <= 10.00) std::cout<<"corrGood: "<<corr<<"  "<<metp4.Px()<<"  "<<tap4.Px()<<std::endl;
    //    else std::cout<<"corrBad: "<<corr<<"  "<<metp4.Px()<<"  "<<tap4.Px()<<std::endl;


    AnaUtil::fillHist1D("muCharge", h1Cand.l1Charge, evWt);
    
    if (h1Cand.l1Charge > 0.0) {
      AnaUtil::fillHist1D ("ZMuDEta+", (ZCand.p4.Eta()-h1l1p4.Eta()), evWt);
      AnaUtil::fillHist1D ("ZMuDPhi+", TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l1p4.Phi()), evWt);
      AnaUtil::fillHist1D ("ZMuDR+", (ZCand.p4).DeltaR(h1l1p4), evWt);
      AnaUtil::fillHist1D ("ZTauDEta+", (ZCand.p4.Eta()-h1l2p4.Eta()), evWt);
      AnaUtil::fillHist1D ("ZTauDPhi+", TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l2p4.Phi()), evWt);
      AnaUtil::fillHist1D ("ZTauDR+", (ZCand.p4).DeltaR(h1l2p4), evWt);
      double dEtaP = std::fabs((ZCand.p4.Eta()-h1l2p4.Eta()) - (ZCand.p4.Eta()-h1l1p4.Eta()));
      double dPhiP = std::fabs(TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l2p4.Phi()) - TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l1p4.Phi()));
      double dRP   = std::fabs((ZCand.p4).DeltaR(h1l2p4) - (ZCand.p4).DeltaR(h1l1p4));
      AnaUtil::fillHist1D ("dEtaP", dEtaP, evWt);
      AnaUtil::fillHist1D ("dPhiP", dPhiP, evWt);
      AnaUtil::fillHist1D ("dRP", dRP, evWt);
    }
    else {
      AnaUtil::fillHist1D ("ZMuDEta-", (ZCand.p4.Eta()-h1l1p4.Eta()), evWt);
      AnaUtil::fillHist1D ("ZMuDPhi-", TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l1p4.Phi()), evWt);
      AnaUtil::fillHist1D ("ZMuDR-", (ZCand.p4).DeltaR(h1l1p4), evWt);
      AnaUtil::fillHist1D ("ZTauDEta-", (ZCand.p4.Eta()-h1l2p4.Eta()), evWt);
      AnaUtil::fillHist1D ("ZTauDPhi-", TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l2p4.Phi()), evWt);
      AnaUtil::fillHist1D ("ZTauDR-", (ZCand.p4).DeltaR(h1l2p4), evWt);
      double dEtaM = std::fabs((ZCand.p4.Eta()-h1l2p4.Eta()) - (ZCand.p4.Eta()-h1l1p4.Eta()));
      double dPhiM = std::fabs(TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l2p4.Phi()) - TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l1p4.Phi()));
      double dRM   = std::fabs((ZCand.p4).DeltaR(h1l2p4) - (ZCand.p4).DeltaR(h1l1p4));
      AnaUtil::fillHist1D ("dEtaM", dEtaM, evWt);
      AnaUtil::fillHist1D ("dPhiM", dPhiM, evWt);
      AnaUtil::fillHist1D ("dRM", dRM, evWt);

    }



    TLorentzVector tap4_3;
    tap4_3.SetPtEtaPhiE((1.0/corr_3)*h1l2p4.Pt(), h1l2p4.Eta(), h1l2p4.Phi(), h1l2p4.E());
    double h2MT_3_corr = Calculate_TotalMT(mup4, tap4_3, metp4);
    AnaUtil::fillHist1D ("h2M_3_corr", h2MT_3_corr, evWt);

    AnaUtil::fillHist1D ("test", (metp4.Px()/metp4.Py() - h1l2p4.Px()/h1l2p4.Py()), evWt);
    
    AnaUtil::fillHist1D ("h1_l1l2DR", h1l1p4.DeltaR(h1l2p4), evWt);
    AnaUtil::fillHist1D ("h1_l1l2DPhi", TVector2::Phi_mpi_pi(h1l1p4.Phi()-h1l2p4.Phi()), evWt);
    AnaUtil::fillHist1D ("h1_l1l2DEta", (h1l1p4.Eta()-h1l2p4.Eta()), evWt);

    AnaUtil::fillHist1D ("ZMuDEta", (ZCand.p4.Eta()-h1l1p4.Eta()), evWt);
    AnaUtil::fillHist1D ("ZMuDPhi", TVector2::Phi_mpi_pi(ZCand.p4.Phi()-h1l1p4.Phi()), evWt);
    AnaUtil::fillHist1D ("ZMuDR", (ZCand.p4).DeltaR(h1l1p4), evWt);

    AnaUtil::fillHist1D ("LT_noTa", (ZCand.l1P4.Pt() + ZCand.l2P4.Pt() + h1Cand.l1P4.Pt()), evWt);
    AnaUtil::fillHist1D ("LT_Ta", (ZCand.l1P4.Pt() + ZCand.l2P4.Pt() + h1Cand.l1P4.Pt() + h1Cand.l2P4.Pt()), evWt);
    AnaUtil::fillHist1D ("METovHT_noTa", (met_->MET)/(ZCand.l1P4.Pt() + ZCand.l2P4.Pt() + h1Cand.l1P4.Pt()), evWt);
    AnaUtil::fillHist1D ("METovHT_Ta", (met_->MET)/(ZCand.l1P4.Pt() + ZCand.l2P4.Pt() + h1Cand.l1P4.Pt() + h1Cand.l2P4.Pt()), evWt);

    if (h1l1p4.DeltaR(h1l2p4) > 2.0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 9, evWt);    
    

    bool hasbJet {false};
    int nbJet {0};

    //hT
    double hT = 0.0;
    if (JetColl.size() > 0){
      for (auto &j: JetColl){
	if (j.BTag   == 1) {
	  hasbJet = true;
	  nbJet++;
	}
	hT += j.PT;
      }
      AnaUtil::fillHist1D ("HT", hT, evWt);
    }
    
    //lT
    
    double lT = 0.0;
    for (auto &e: ElectronColl){
      lT += e.PT;
    }
    for (auto &m: MuonColl){
      lT += m.PT;
    }
    AnaUtil::fillHist1D ("LT", lT, evWt);
    AnaUtil::fillHist1D ("ST", (lT+hT), evWt);
    
    if (JetColl.size() > 0){
      AnaUtil::fillHist1D ("Jet1Pt_stage2", JetColl[0].PT, evWt);
      // if (JetColl[0].PT > 80) continue;
      //      AnaUtil::fillHist1D("evtCutFlow", 8, evWt);
    }

    AnaUtil::fillHist1D ("met_PreCut", met_-> MET, evWt);
    //    if ( (met_-> MET) > AnaUtil::cutValue(evselCutMap(), "maxMET")) continue;
    //    AnaUtil::fillHist1D ("LT_PreCut", (lT+hT), evWt);
    AnaUtil::fillHist1D ("LT_PreCut", lT, evWt);
    if (lT < 200.0) continue;
    AnaUtil::fillHist1D("evtCutFlow", 10, evWt);

    AnaUtil::fillHist1D ("nbJets", nbJet, evWt);
    if (hasbJet) continue;
    AnaUtil::fillHist1D("evtCutFlow", 11, evWt);

    //_________________________________________End of PreSelection_______________________________________________//

    histf->cd();

    //_____________________________________________MVA Skimming__________________________________________________//

    if (skimObj_) {
      TreeVariables varList;

      varList.met          = met_->MET;
      varList.hT           = hT;
      
      skimObj_->fill(varList);
    }
    //_________________________________________MVA Evaluation_______________________________________________//

    histf->cd();
    double mvaOut = -999.9;

    if (_readMVA) {
      InputVariables varlist;

      varlist.met          = met_->MET;
      varlist.hT           = hT;

      mvaOut = _mvaObj->evaluate(_MVAnetwork, varlist);
    }

    histf->cd();
    histf->cd("Analysis");

    AnaUtil::fillHist1D ("mvaOutput", mvaOut, evWt);
    
    
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
      "0) Events processed: ",
      "1) nEle || nMu >=2 : ",
      "2) hLepPt >= 70    : ",
      "3) no ZCand        : ",
      "4) 75<ZCandMass<105: ",
      "5) nMu+nEle >= 3   : ",
      "6) nMu >= 1        : ",	
      "7) nTau = 1        : ",
      "8) have hCand      : ",
      "9) h2_l1l2dR <= 2.0: ",
      "10) lepSumPt > 200 : ",
      "11)has no b Jet    : "
	  };
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection");  
  double lumiFac = 1.0;
  lumiFac = lumiWt(evtWtSum_, maxEvt_);
  cout << endl
       << "evtWeightSum: " << setw(10) << setprecision(0) << evtWtSum_ << endl
       << "      lumiWt: " << setw(10) << setprecision(5) << lumiFac
       << endl;
  AnaUtil::scaleHistogram("evtCutFlow", lumiFac);
  //  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection (Weighted)", "Events", 1);
  AnaUtil::showEfficiency("evtCutFlow", evLabels, "Event Selection (lumiScaled)");
}


void ExoHiggsAna::closeHistFiles(){
  histf->cd();
  histf->Write();
  histf->Close();
}
void ExoHiggsAna::closeFiles(){
  closeHistFiles();
  if (skimObj_ != nullptr) skimObj_->close();
  delete treeReader;
  delete chain;
}
bool ExoHiggsAna::readJob(const string& jobFile, int& nFiles)
{
  // Open the file containing the datacards
  ifstream fin(jobFile.c_str(), ios::in);    
  if (!fin) {
    cerr << "==> Input File: <<" << jobFile << ">> could not be opened!" << endl;
    return false;
  }

  static constexpr int BUF_SIZE = 256;
  char buf[BUF_SIZE];
  while (fin.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
    string line(buf);
    if (line.empty() || line == "START") continue;   
    
    // enable '#' and '//' style comments
    if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
    if (line == "END") break;
    
    // Split the line into words
    vector<string> tokens;
    AnaUtil::tokenize(line, tokens);
    int vsize = tokens.size();
    assert(vsize > 1);
    
    const string& key   = tokens.at(0);
    const string& value = tokens.at(1);
    if (key == "dataType") {
      string vtmp(value);
      std::transform(vtmp.begin(), vtmp.end(), vtmp.begin(), ::toupper);
      vector<string> dt;
      AnaUtil::tokenize(vtmp, dt, "#");
      if (dt.size()) {
        isMC_ = (dt.at(0) == "MC") ? true : false;
        if (isMC_ && dt.size() > 1) {
          isSignal_ = (dt.at(1) == "SIGNAL") ? true : false;
        }
      }
    }
    else if (key == "readGenInfo") 
      readGenInfo_ = std::stoi(value.c_str()) > 0 ? true : false;
    else if (key == "inputFile")
      AnaUtil::buildList(tokens, fileList_);
    else if (key == "maxEvent") 
      maxEvt_ = std::stoi(value.c_str());
    else if (key == "histFile")
      histFile_ = value;
    else if (key == "readMVA")
      _readMVA = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "MVAnetwork")
      _MVAnetwork = value;
    else if (key == "MVAxmlFile")
      _MVAxmlFile = value;
    else if (key == "createMVATree")
      _createMVATree = (atoi(value.c_str()) > 0) ? true : false;
    else if (key == "mvaInputFile")
      _mvaInputFile = value;
    else if (key == "eventId" && tokens.size() == 4)
      AnaUtil::buildMap(tokens, eventIdMap_);
    else {
      if (0) cout << "==> " << line << endl;
    AnaUtil::storeCuts(tokens, hmap_);
    }						
  }
  // Close the file
  fin.close();
  //  if (!isSignal_) readGenInfo_ = false;
  chain = new TChain("Delphes");
  // Build the chain of root files
  for (const auto& fname: fileList_) {
    cout << ">>> INFO. Adding input file " << fname << " to TChain " << endl;
    ++nFiles;
    chain->Add(fname.c_str());
    int nevt = static_cast<int>(chain->GetEntries());
    if (maxEvt_ > 0 && nevt >= maxEvt_) break;
  }
  if (!nFiles) {
    cerr << ">>> WARN. Input Root file list is empty! exiting ..." << endl;
    return false;
  } 
  //  cout<<nFiles<<endl;
  treeReader = new ExRootTreeReader(chain);
  if(!treeReader) return false;

  printJob();  
  return true;
}

void ExoHiggsAna::printJob(ostream& os) const{
  os << endl;
  AnaUtil::showCuts(hmap_, os);
}

//******************************************Some More Functions (Make AnaBase and move them)***********************************************//

double ExoHiggsAna::lumiWt(double evtWeightSum, int maxEvents) const{
  //  double nevt = (evtWeightSum > -1) ? evtWeightSum : AnaUtil::cutValue(lumiWtMap(), "nevents");
  double nevt = (evtWeightSum > -1) ? evtWeightSum : maxEvents;
  std::cout << "-- intLumi: " << AnaUtil::cutValue(lumiWtMap(), "intLumi")
	    << " xsec: " << AnaUtil::cutValue(lumiWtMap(), "xsec")
	    << " nevt: " << nevt << std::endl;
  return (AnaUtil::cutValue(lumiWtMap(), "intLumi") * AnaUtil::cutValue(lumiWtMap(), "xsec") / nevt);
}

int ExoHiggsAna::getMotherId(const GenParticle& gp, TClonesArray *gen, int& mmid) const {
  int pdgid = gp.PID;
  int indx = gp.M1;
  if (indx < 0) return -1;
  GenParticle& mgp = dynamic_cast<GenParticle&> (*(gen->At(indx)));
  mmid = mgp.PID;
  if (std::abs(mmid) == 2212) return -1;
  while (mmid == pdgid) {
    indx = mgp.M1;
    if (indx < 0) continue;
    mgp  = dynamic_cast<GenParticle&>(*(gen->At(indx)));
    mmid = mgp.PID;
  }
  return indx;
}

bool ExoHiggsAna::jetLeptonCleaning(const Jet& jet, std::vector <Muon> muList_, std::vector <Electron> eleList_, double dR) {
  TLorentzVector jetP4 = jet.P4();
  for (auto& m: muList_){
    AnaUtil::fillHist1D ("jetMuDR", jetP4.DeltaR(m.P4()), 1.0);
    if (jetP4.DeltaR(m.P4()) <= dR) return false;
  }
  for (auto& e: eleList_){
    AnaUtil::fillHist1D ("jetEleDR", jetP4.DeltaR(e.P4()), 1.0);
    if (jetP4.DeltaR(e.P4()) <= dR) return false;
  }
  return true;
}

void ExoHiggsAna::dumpGenInfo(TClonesArray *gen, ostream& os) {
  int ngen = gen->GetEntries();
  if (ngen == 0) return;

  os << setprecision(2);
  os << " -- # GenParticle: " << ngen << endl;
  os << "indx    status    pdgId     eta      phi      pt     energy  moIndx"
     << "      moID                   daughterID"
     << endl;
  int indx = 0;
  for (TObject* obj: *gen) {
    //const GenParticle& gp = *(dynamic_cast<GenParticle*> (obj));                                                                                                      
    const GenParticle& gp = dynamic_cast<const GenParticle&> (*obj);
    //int abs_pid = std::abs(gp.PID);                                                                                                                                   
    std::ostringstream mID;
    vector<int> m;
    m.push_back(gp.M1);
    m.push_back(gp.M2);
    for (int mi: m) {
      if (mi < 0 || mi >= ngen) continue;
      const GenParticle& mgp = dynamic_cast<const GenParticle&>(*(gen->At(mi)));
      mID << " " << mgp.PID;
    }
    string ms = mID.str();
    if (!ms.length()) ms = " -";

    std::ostringstream dID;
    vector<int> d;
    d.push_back(gp.D1);
    d.push_back(gp.D2);
    for (int di: d) {
      if (di < 0 || di >= ngen) continue;
      const GenParticle& dgp = dynamic_cast<const GenParticle&>(*(gen->At(di)));
      //  const GenParticle& dgp = (GenParticle*)gen->At(di);                                                                                                           
      double energy = dgp.E;
      int pdgid = dgp.PID;
      if (std::abs(pdgid) == 21 && energy <= 10) continue;
      dID << " " << dgp.PID;
    }
    string ds = dID.str();
    if (!ds.length()) ds = " -";
    os << setw(4)  << indx++
       << setw(8)  << gp.Status
       << setw(10) << gp.PID
       << setw(10) << gp.Eta
       << setw(9)  << gp.Phi
       << setw(9)  << gp.PT
       << setw(9)  << gp.E
       << setw(8)  << gp.M1
       << setw(10) << ms
       << setw(28) << ds
       << endl;
  }
}
//******************************************Some More Functions (Make AnaBase and move these)***********************************************//
