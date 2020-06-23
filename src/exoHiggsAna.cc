#include <iostream>
#include <string>
#include <memory> 
#include <utility>
#include <functional>
#include <utility>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <climits>
#include <cassert>
#include <cstdlib>
#include <sstream>
#include "TROOT.h"
#include "TStopwatch.h"
#include "TSystem.h"
#include "TApplication.h"

#include "TString.h"
#include "string.h"

#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "/home/sinpcms/GOURAB/Delphes-3.4.2/classes/DelphesClasses.h"

#include "/home/sinpcms/GOURAB/Delphes-3.4.2/external/ExRootAnalysis/ExRootTreeReader.h"
#include "/home/sinpcms/GOURAB/Delphes-3.4.2/external/ExRootAnalysis/ExRootTreeWriter.h"
#include "/home/sinpcms/GOURAB/Delphes-3.4.2/external/ExRootAnalysis/ExRootTreeBranch.h"
#include "/home/sinpcms/GOURAB/Delphes-3.4.2/external/ExRootAnalysis/ExRootResult.h"
#include "/home/sinpcms/GOURAB/Delphes-3.4.2/external/ExRootAnalysis/ExRootUtilities.h"


using namespace std;

#include "ExoHiggsAna.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;

int main(int argc, char* argv[]) {
  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " jobFile " << endl;
    exit(1);
  }     
  //  string inFile (argv[1]);
  //  string outFile(argv[2]);
  string jobFile(argv[1]);

  gROOT->SetBatch(kTRUE);
  
  // Create analysis object 
  cout << "=== Start of Analysis === " << endl;

  ExoHiggsAna anaH;

  TStopwatch timer;
  cout << "==> Start event loop now with " << endl;
  timer.Start();

  int nFiles = 0;
  bool succeed = anaH.readJob(jobFile, nFiles);
  if (!succeed) exit(2);

  //Initialize analysis
  if (!anaH.beginJob()) exit(4);

  // Analysis over
  anaH.endJob();
  anaH.closeFiles();
  cout << "=== End of Analysis === " << endl;
  timer.Stop();
  cout << "Realtime/CpuTime = " << timer.RealTime() << "/" << timer.CpuTime() << endl;
  return 0;
}
