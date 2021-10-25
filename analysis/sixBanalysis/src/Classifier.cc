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

#include "Classifier.h"

Classifier::Classifier(TString name,TString title)
{
  _name = name;
  _title = title;
}

void Classifier::add(TString entry,float value)
{
  if (_classifier.count(entry) == 0) {
    _entries.push_back(entry);
    _classifier[entry] = 0;
  }
  _classifier[entry] += value;
}

void Classifier::write(TFile& output)
{
  output.cd();
  unsigned int nentries = _entries.size();
  TH1F classifier(_name,_title,nentries,0,nentries);
  for (unsigned int i = 0; i < nentries; i++)
    {
      TString entry = _entries[i];
      float value = _classifier[entry];

      classifier.SetBinContent(i+1,value);
      classifier.GetXaxis()->SetBinLabel(i+1,entry);
    }

  classifier.Write();
}
