#ifndef ELECTRON_H
#define ELECTRON_H

#include "Candidate.h"

class Electron : public Candidate
{
public:
  Electron () : Candidate(){typeId_=11;}
  Electron (int idx, NanoAODTree* nat) : Candidate(idx, nat){typeId_=11; buildP4();}
  ~Electron(){};
  std::unique_ptr<Candidate> clone() const override{
    Electron *clonedElectron = new Electron(this->getIdx(), this->getNanoAODTree());
    clonedElectron->setP4(this->P4());
    return std::unique_ptr<Electron> (clonedElectron);
  }

  float get_E() const { return this->P4().E(); }
  float get_m() const { return this->P4().M(); }
  float get_pt() const { return this->P4().Pt(); }
  float get_eta() const { return this->P4().Eta(); }
  float get_phi() const { return this->P4().Phi(); }
  float get_dxy() const { return get_property((*this), Electron_dxy); }
  float get_dz() const { return get_property((*this), Electron_dz); }
  float get_charge() const { return get_property((*this), Electron_charge); }
  float get_pfRelIso03_all() const { return get_property((*this), Electron_pfRelIso03_all); }
  bool get_mvaFall17V2Iso_WPL() const { return get_property((*this), Electron_mvaFall17V2Iso_WPL); }
  bool get_mvaFall17V2Iso_WP90() const { return get_property((*this), Electron_mvaFall17V2Iso_WP90); }
  bool get_mvaFall17V2Iso_WP80() const { return get_property((*this), Electron_mvaFall17V2Iso_WP80); }
  
private:
  void buildP4() override; 
};

#endif
