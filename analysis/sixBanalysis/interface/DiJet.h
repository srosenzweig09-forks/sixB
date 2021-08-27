#ifndef DIJET_H
#define DIJET_H

#include "Jet.h"

class DiJet
{
public:
  DiJet(Jet& j1,Jet& j2);
  
  void set_signalId(int id) {  signalId = id; }
  
  int get_signalId() const { return signalId; }
  float Pt() const         { return p4.Pt(); }
  float Eta() const        { return p4.Eta(); }
  float Phi() const        { return p4.Phi(); }
  float M() const          { return p4.M(); }
  float E() const          { return p4.E(); }

private:
  int signalId = -1;
  p4_t p4;
  
};

#endif
