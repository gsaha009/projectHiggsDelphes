#include <iostream>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "MVASkim.h"

using std::string;
using std::cout;
using std::endl;

MVASkim::MVASkim(const string& filename) {
  _mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  _tree = new TTree("RTree", "RTree");
  
  _tree->Branch("hLepPt", &_varList.hLepPt, "hLepPt/F");
  _tree->Branch("met", &_varList.met, "met/F");
  _tree->Branch("mTovSt", &_varList.mTovSt, "mTovSt/F");
  _tree->Branch("muTaDR", &_varList.muTaDR, "muTaDR/F");
  _tree->Branch("muZDR", &_varList.muZDR, "muZDR/F");
  _tree->Branch("colTest", &_varList.colTest, "colTest/F");
  _tree->Branch("colMass", &_varList.colMass, "colMass/F");
  _tree->Branch("tauPt", &_varList.tauPt, "tauPt/F");
  _tree->Branch("ST", &_varList.ST, "ST/F");

  _mvaFile->ls();
}
MVASkim::~MVASkim() {
  if (_tree) delete _tree;  
  if (_mvaFile) delete _mvaFile;
}
void MVASkim::fill(const TreeVariables& varList) {
  memcpy(&_varList, &varList, sizeof(varList));
  _mvaFile->cd();
  _tree->Fill();
}
void MVASkim::close() {
  //_mvaFile = TFile::Open(filename.c_str(), "RECREATE", "Skimmed Tree");
  _mvaFile->cd();
  _tree->Print();
  _tree->Write();
  _mvaFile->Save();
  _mvaFile->Write();
  //_mvaFile->Close();
}
