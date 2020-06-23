#ifndef __ANAUTIL__HH
#define __ANAUTIL__HH

#include <string>
#include <vector>
#include <map>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TH1K.h"

using std::ostream;

namespace AnaUtil {
  template <typename T>
  T deltaPhiT(T phi1, T phi2) {
    T result = phi1 - phi2;
    while (result > TMath::Pi()) result -= 2*TMath::Pi();
    while (result <= -TMath::Pi()) result += 2*TMath::Pi();
    return result;
  }  
  double deltaPhi(double phia, double phib);
  double deltaPhi(const TLorentzVector& a, const TLorentzVector& b);
  double deltaR  (const TLorentzVector& a, const TLorentzVector& b);
  bool sameObject(const TLorentzVector& lv1, const TLorentzVector& lv2);
  void buildList(const std::vector<std::string>& tokens, std::vector<std::string>& list);
  double cutValue(const std::map<std::string, double>& m, const std::string& cname);
  const std::map<std::string, double>& cutMap(const std::map<std::string, std::map<std::string, double>>& hmap, const std::string& mkey);
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters=" ");
  void buildMap(const std::vector<std::string>& tokens, std::map<std::string, int>& hmap);
  void storeCuts(const std::vector<std::string>& tokens, std::map<std::string, std::map<std::string, double>>& hmap);
  void showCuts(const std::map<std::string, std::map<std::string, double> >& hmap, ostream& os=std::cout);
  

  // ------------------------------------------------------------------------
  // Convenience routine for filling 1D histograms. We rely on root to keep 
  // track of all the histograms that are booked all over so that we do not 
  // have to use any global variables to save the histogram pointers. Instead, 
  // we use the name of the histograms and gROOT to retrieve them from the 
  // Root object pool whenever necessary. This is the closest one can go to 
  // hbook and ID based histogramming
  // -------------------------------------------------------------------------
  template <class T> 
  TLorentzVector getP4(const T& obj) {
    TLorentzVector lv;
    lv.SetPtEtaPhiE(obj.pt, obj.eta, obj.phi, obj.energy);
    return lv;
  }

  void printP4(const TLorentzVector& lv, const std::string& tag, std::ostream& os=std::cout);

  TH1* getHist1D(const char* hname);
  TH1* getHist1D(const std::string& hname);
  template <class T>
  bool fillHist1D(const char* hname, T value, double w=1.0) {
    TH1* h = getHist1D(hname);
    if (h == nullptr) return false;
    h->Fill(value, w);
    return true;
  }

  template <class T>
    bool fillHist1D(const std::string& hname, T value, double w=1.0) {
    return fillHist1D(hname.c_str(), value, w);
  }

  // ---------------------------------------------
  // Convenience routine for filling 2D histograms
  // ---------------------------------------------
  TH2* getHist2D(const char* hname);
  TH2* getHist2D(const std::string& hname);
  template <class T1, class T2>
  bool fillHist2D(const char* hname, T1 xvalue, T2 yvalue, double w=1.0) {
    TH2* h = getHist2D(hname);
    if (h == nullptr) return false;
    h->Fill(xvalue, yvalue, w);
    return true;
  }

  template <class T1, class T2>
  bool fillHist2D(const std::string& hname, T1 xvalue, T2 yvalue, double w=1.0) {
    return fillHist2D(hname.c_str(), xvalue, yvalue, w);
  }
  // ---------------------------------------------
  // Convenience routine for filling 3D histograms
  // ---------------------------------------------
  TH3* getHist3D(const char* hname);
  TH3* getHist3D(const std::string& hname);

  template <class T1, class T2, class T3>
  bool fillHist3D(const char* hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) {
    TH3* h = getHist3D(hname);
    if (h == nullptr) return false;
    h->Fill(xvalue, yvalue, zvalue, w);
    return true;
  }
  template <class T1, class T2, class T3>
  bool fillHist3D(const std::string& hname, T1 xvalue, T2 yvalue, T3 zvalue, double w=1.0) {
    return fillHist3D(hname.c_str(), xvalue, yvalue, zvalue, w);
  }

  // --------------------------------------------------
  // Convenience routine for filling profile histograms
  // --------------------------------------------------
  TProfile* getProfile(const char* hname);
  TProfile* getProfile(const std::string& hname);

  bool fillProfile(const char *hname, float xvalue, float yvalue, double w=1.0);
  bool fillProfile(const std::string& hname, float xvalue, float yvalue, double w=1.0);

  void scaleHistogram(const std::string& hname, double fac);
  void showEfficiency(const std::string& hname, 
		      const std::vector<std::string>& slist, 
		      const std::string& header, 
		      const std::string& tag="Events", std::ostream& os=std::cout);
  void showCount(const std::string& hname, 
		 const std::vector<std::string>& slist, 
		 const std::string& tag, int prec=0, std::ostream& os=std::cout);

  
}
#endif
