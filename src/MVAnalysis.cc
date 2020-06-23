#include <iostream>
#include <memory>
#include "MVAnalysis.h"

using std::string;
using std::cout;
using std::endl;

MVAnalysis::MVAnalysis(const string& mva_algo, const string& xmlfile) 
{
  //  reader_ = std::make_unique<TMVA::Reader>(new TMVA::Reader("!Silent"));
  reader_ = new TMVA::Reader("!Silent");

  reader_->AddVariable("hLepPt", &varList_.hLepPt);
  reader_->AddVariable("met", &varList_.met);
  reader_->AddVariable("mTovSt", &varList_.mTovSt);
  reader_->AddVariable("muTaDR", &varList_.muTaDR);
  reader_->AddVariable("muZDR", &varList_.muZDR);
  reader_->AddVariable("tauPt", &varList_.tauPt);
  
 reader_->BookMVA(mva_algo.c_str(), xmlfile);
}
double MVAnalysis::evaluate(const string& mva_algo, const InputVariables& varList) {
  memcpy(&varList_, &varList, sizeof(varList)); // use move syntax here
  return reader_->EvaluateMVA(mva_algo.c_str());
}
