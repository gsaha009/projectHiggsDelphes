#ifndef __ExoHiggsAna__hh
#define __ExoHiggsAna__hh

#define NEL(x) (sizeof((x))/sizeof((x)[0]))

#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <bitset>

#include "TROOT.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TString.h"
#include "THStack.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"

#include "AnaUtil.h"

#include "classes/DelphesClasses.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

template <class T>
class PtComparator {
 public:
  bool operator()(const T &a, const T &b) const {
    return a.PT > b.PT;
  }
};

template <class T>
class MassComparator {
 public:
  bool operator()(const T &a, const T &b) const {
    TLorentzVector l1,l2;
    l1.SetPtEtaPhiE(a.pt, a.eta, a.phi, a.energy);
    l2.SetPtEtaPhiE(b.pt, b.eta, b.phi, b.energy);
    return l1.M() > l2.M();
  }
};



class ExoHiggsAna {
    
public:
  ExoHiggsAna();
  virtual ~ExoHiggsAna();
  
  void eventLoop(ExRootTreeReader *treeReader);  // the main analysis 
  bool beginJob(const std::string& inputFile, std::string& outputFile);
  void endJob();
  //virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  //  void printJob(std::ostream& os=std::cout) const;
  void bookHistograms();
  void closeHistFiles();
  void closeFiles();
  void clearLists();
  //  bool hasZcandidate (const std::vector<LeptonCandidate>& lepColl, double puWt);
  
  //  std::vector <Jet> bJetList;
  //  std::vector <Jet> jetList;
  //  std::vector <GenParticle> Gens_L;
  std::vector <Muon> MuonColl;
  std::vector <Electron> ElectronColl;
  std::vector <Jet> JetColl;
  //  std::vector <std::pair<Jet, Jet>> jetPairSelected;	
  //  std::set <Jet> signalLikeJets;

  TFile *histf;
  ExRootTreeReader *treeReader;
  TChain *chain;

};
#endif
