#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include <iostream>
#include <string>
#include <iomanip>
#include <any>
#include <chrono>
#include <map>

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TString.h"

class Classifier {
public:
  Classifier(TString name="h_classifier",TString title="Selection Classifier");
  void add(TString entry,float value=1);
  void write(TFile& output);
private:
  TString _name;
  TString _title;
  std::vector<TString>    _entries;
  std::map<TString,float> _classifier;
};

#endif // CLASSIFIER_H
