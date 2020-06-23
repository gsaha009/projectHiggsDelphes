#ifndef __ExoHiggsAna_ZZ4L__hh
#define __ExoHiggsAna_ZZ4L__hh

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
#include "ZCandidate.h"

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


template <typename T>
void ZSelector(const std::vector<T>& LeptonList, std::vector<ZSpace::ZCandidate>& candList) {
  for (unsigned int i = 0; i < LeptonList.size(); ++i) {
    const auto& ip = LeptonList[i];

    TLorentzVector lep1P4 (ip.P4());
    for (unsigned int j = i+1; j < LeptonList.size(); ++j) {
      const auto& jp = LeptonList[j];
      // require opposite charges                                                           
      if ((ip.Charge * jp.Charge) > 0) continue;
      TLorentzVector lep2P4 (jp.P4());
      if (lep1P4.DeltaR(lep2P4) < 0.02) continue; // new 
                       
      ZSpace::ZCandidate ztmp;
      if (typeid(jp) == typeid(Muon)) {
	ztmp.flavour = static_cast<int>(ZSpace::ZType::mumu);
      }
      else if (typeid(jp) == typeid(Electron)) {
	ztmp.flavour = static_cast<int>(ZSpace::ZType::ee);
      }
      else
	ztmp.flavour = static_cast<int>(ZSpace::ZType::unkwn);

      ztmp.l1Index = i;
      ztmp.l1P4 = lep1P4;
      ztmp.l1Charge = ip.Charge;

      ztmp.l2Index = j;
      ztmp.l2P4 = lep2P4;
      ztmp.l2Charge = jp.Charge;

      TLorentzVector p4 = lep1P4 + lep2P4;
      double Zmass = p4.M();
      ztmp.p4 = p4;
      ztmp.mass = Zmass;
      ztmp.massDiff = std::fabs(Zmass - ZSpace::MZnominal);
      ztmp.dEtall = ztmp.l1P4.Eta() - ztmp.l2P4.Eta();
      ztmp.dPhill = TVector2::Phi_mpi_pi(ztmp.l1P4.Phi() - ztmp.l2P4.Phi());
      ztmp.dRll   = ztmp.l1P4.DeltaR(ztmp.l2P4);

      candList.push_back(ztmp);
    }
  }
}


// find an opposite signed lepton-lepton pair with different flavors              
template <typename T1, typename T2>
  std::bitset<6> leptonPairSelector(const std::vector<T1>& Lepton1List, const std::vector<T2>& Lepton2List, const ZSpace::ZCandidate& eventZ, std::vector<ZSpace::ZCandidate>& leptonPairList, double modelMass, int EVwt, double minDR)
{
  int flags [] {0,0,0,0,0,0};
  for (size_t i = 0; i < Lepton1List.size(); ++i) {
    const auto& ielem = Lepton1List[i];
    TLorentzVector lep1P4 = ielem.P4();
    // Distinct object                                                      
    if (AnaUtil::sameObject(lep1P4, eventZ.l1P4) || AnaUtil::sameObject(lep1P4, eventZ.l2P4)) {
      ++flags[0];
      continue;
    }

    // Angular separation of the lepton with the leptons from Z               
    AnaUtil::fillHist1D("dR_l1Zl1", lep1P4.DeltaR(eventZ.l1P4), EVwt);
    AnaUtil::fillHist1D("dR_l1Zl2", lep1P4.DeltaR(eventZ.l2P4), EVwt);
    AnaUtil::fillHist1D("dR_l1Z", lep1P4.DeltaR(eventZ.p4), EVwt);
    if (lep1P4.DeltaR(eventZ.l1P4) < 0.1 || lep1P4.DeltaR(eventZ.l2P4) < 0.1) {
      ++flags[1];
      continue;
    }
    for (size_t j = 0; j < Lepton2List.size(); ++j) {
      const auto& jelem = Lepton2List[j];
      TLorentzVector lep2P4 = jelem.P4();
      // Distinct object                                                      
      if (AnaUtil::sameObject(lep2P4, eventZ.l1P4) || AnaUtil::sameObject(lep2P4, eventZ.l2P4)) {
	++flags[2];
	continue;
      }

      // Angular separation of the lepton with the leptons from Z
      AnaUtil::fillHist1D("dR_l2Zl1", lep2P4.DeltaR(eventZ.l1P4), EVwt);
      AnaUtil::fillHist1D("dR_l2Zl2", lep2P4.DeltaR(eventZ.l2P4), EVwt);
      AnaUtil::fillHist1D("dR_l2Z", lep2P4.DeltaR(eventZ.p4), EVwt);
      if (lep2P4.DeltaR(eventZ.l1P4) < 0.1 || lep2P4.DeltaR(eventZ.l2P4) < 0.1) {
	++flags[3];
	continue;
      }
     
      // Require oppositely charged leptons                  
      if (ielem.Charge * jelem.Charge > 0) {
	++flags[4];
	continue;
      }
      
      // Angular separation between the two leptons being considered 
      AnaUtil::fillHist1D("dR_l1l2", lep1P4.DeltaR(lep2P4), EVwt);
      if (lep1P4.DeltaR(lep2P4) < minDR) {
	++flags[5];
	continue;
      }
      ZSpace::ZCandidate llp;
      if (typeid(ielem) == typeid(Muon))
	llp.flavour = static_cast<int>(ZSpace::llType::mue);
      else if (typeid(ielem) == typeid(Electron))
	llp.flavour = static_cast<int>(ZSpace::llType::emu);
      else
	llp.flavour = static_cast<int>(ZSpace::ZType::unkwn);
     
      //lep1 candidate
      llp.l1Index = i;
      llp.l1P4 = lep1P4;
      llp.l1Charge = ielem.Charge;
      //lep2 candidate
      llp.l2Index = j;
      llp.l2P4 = lep2P4;
      llp.l2Charge = jelem.Charge;

      TLorentzVector p4 = llp.l1P4 + llp.l2P4;
      llp.p4   = p4;
      llp.mass = p4.M();
      llp.massDiff = (modelMass > 0) ? std::fabs(llp.mass - modelMass) : -999;
      llp.dEtall   = llp.l1P4.Eta() - llp.l2P4.Eta();
      llp.dPhill   = TVector2::Phi_mpi_pi(llp.l1P4.Phi() - llp.l2P4.Phi());
      llp.dRll     = llp.l1P4.DeltaR(llp.l2P4);

      leptonPairList.push_back(llp);
    }
  }
  std::bitset<6> result;
  int i = 0;
  for (auto v: flags) {
    ((v > 0) ? result.set(i, 1) : AnaUtil::fillHist1D("leptonPairCutFlow", i, EVwt));
    i++;
  }
  return result;
}



//Main Analysis Class Description
class ExoHiggsAna_ZZ4L {
    
public:
  ExoHiggsAna_ZZ4L();
  virtual ~ExoHiggsAna_ZZ4L();
  
  void eventLoop(ExRootTreeReader *treeReader);  // the main analysis 
  bool beginJob(const std::string& inputFile, std::string& outputFile);
  void endJob();
  //virtual bool readJob(const std::string& jobFile, int& nFiles) override;
  //  void printJob(std::ostream& os=std::cout) const;
  void bookHistograms();
  void closeHistFiles();
  void closeFiles();
  void clearLists();
  int getMotherId(const GenParticle& gp, TClonesArray *gen, int& mmid) const;

  static bool dmComparator(const ZSpace::ZCandidate& a, const ZSpace::ZCandidate& b) {
    return (a.massDiff < b.massDiff);
  }
  static bool massComparator(const ZSpace::ZCandidate& a, const ZSpace::ZCandidate& b) {
    return (a.mass > b.mass);
  }

  std::vector <Muon> MuonColl;
  std::vector <Electron> ElectronColl;
  std::vector <Jet> JetColl;

  TFile *histf;
  ExRootTreeReader *treeReader;
  TChain *chain;

};
#endif
