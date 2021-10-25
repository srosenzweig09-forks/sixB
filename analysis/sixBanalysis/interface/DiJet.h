#ifndef DIJET_H
#define DIJET_H

#include "Jet.h"

class DiJet
{
public:
  DiJet(Jet& j1,Jet& j2);
  
  void set_signalId(int id) {  signalId = id; }
  void set_id1(int id) {  id1 = id; }
  void set_id2(int id) {  id2 = id; }
  void set_2j_score(float score) { n_2j_score = score; }
  
  int get_signalId() const   { return signalId; }
  int get_id1() const   { return id1; }
  int get_id2() const   { return id2; }
  float get_2j_score() const { return n_2j_score; }
  float Pt() const           { return p4.Pt(); }
  float Eta() const          { return p4.Eta(); }
  float Phi() const          { return p4.Phi(); }
  float M() const            { return p4.M(); }
  float E() const            { return p4.E(); }
  float dR() const           { return dr_; }

private:
  int signalId = -1;
  int id1 = -1;
  int id2 = -1;
  float n_2j_score = -1;
  float dr_ = -1;
  p4_t p4;
  
};

#endif
