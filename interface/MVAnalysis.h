#ifndef __MVAnalysis__h
#define __MVAnalysis__h

#include <fstream>
#include <string>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

typedef struct  
{
  float hLepPt;
  float met;
  float mTovSt;
  float muTaDR;
  float muZDR;
  float tauPt;
  
} InputVariables;

class MVAnalysis {
    
public:

  MVAnalysis(const std::string& mva_algo, const std::string& xmlfile);
  virtual ~MVAnalysis() {}

  double evaluate(const std::string& tag, const InputVariables& varList);

  InputVariables varList_;
  //  std::unique_ptr<TMVA::Reader> reader_;
  TMVA::Reader* reader_;
};
#endif
