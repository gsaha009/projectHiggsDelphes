#ifndef __MVASkim__h
#define __MVASkim__h

#include <fstream>
#include <string>

class TTree;
class TFile;

typedef struct  
{
  float hLepPt;
  float met;
  float mTovSt;
  float muTaDR;
  float muZDR;
  float colTest;
  float colMass;
  float tauPt;
  float ST;

} TreeVariables;

class MVASkim {
    
public:

  MVASkim(const std::string& filename);
  virtual ~MVASkim();

  void fill(const TreeVariables& varList);
  void close();

  TFile* _mvaFile;
  TTree* _tree;

  TreeVariables _varList;
};
#endif
