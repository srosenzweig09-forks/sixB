/*
 * Author: Suzanne Rosenzweig, s.rosenzweig@cern.ch
 * Purpose: To filter the selected events
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <numeric>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

#include "TH1F.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TChain.h"
#include "TString.h"

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

#include "Math/VectorUtil.h"
#include "Math/Vector3D.h"
#include "Math/Functions.h"
#include "Math/Vector4D.h"
typedef ROOT::Math::PtEtaPhiMVector p4_t;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

void info(string output) {
  string info = "[\e[38;5;11mINFO\e[39m] .. ";
  cout << info << output << endl;
}

void info(string output, string filename) {
  string info = "[\e[38;5;11mINFO\e[39m] .. ";
  string finfo = "[\e[38;5;63mFILE\e[39m]    -> ";
  cout << info << output;
  if (filename.size() < 30) {
    cout << " : " << filename << endl;
  }
  else {
    cout << endl << finfo << filename << endl;
  }
}

int main(int argc, char** argv){
  cout << endl << ".. starting program" << endl << endl;

  po::options_description desc("skim options");
  desc.add_options()
  // ("help", "produce help message")
  // ("input" , po::value<string>()->required(), "input file list")
  ("is-data", po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a data sample (default is false)")
  ("is-signal", po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a HH signal sample (default is false)")
  ("score-file",  po::value<string>()->default_value(""), "score file name without extension")
  ("MX", po::value<int>()->default_value(0), "MX value if is-signal flag is passed")
  ("MY", po::value<int>()->default_value(0), "MY value if is-signal flag is passed")
  ("inner-radius", po::value<int>()->default_value(30), "signal region maximum distance from mH")
  ("outer-radius", po::value<int>()->default_value(40), "control region maximum distance from mH")
  ;

  po::variables_map opts;
  try {
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), opts);
    po::notify(opts);
  }
  catch (po::error& e) {
    cerr << "** [ERROR] " << e.what() << endl;
    return 1;
  }

  const bool is_data = opts["is-data"].as<bool>();
  const bool is_signal = (is_data ? false : opts["is-signal"].as<bool>());

  const int MX = opts["MX"].as<int>();
  const int MY = opts["MY"].as<int>();

  const string sfile = opts["score-file"].as<string>();

  // cout << "[INFO] ... opening file list : " << opts["input"].as<string>() << endl;
  // if ( access( opts["input"].as<string>(), F_OK ) == -1 ){
  //   cerr << "** [ERROR] The input file list does not exist, aborting" << endl;
  //   return 1;        
  // }

  string base = "root://cmseos.fnal.gov//store/user/srosenzw/sixb/ntuples/Summer2018UL/";

  string txtfile, infile, outfile;
  fstream newfile;

  if (is_signal) {
    string mx_my_str;

    if (sfile == "") {
      if (MX==0 || MY==0) {
        cerr << "** [ERROR] flagged as signal but no MX and/or MY value(s) provided!" << endl;
        return 1;
      }

      ostringstream masses;
      masses << "MX_" << MX << "_MY_" << MY;
      mx_my_str = masses.str();
    }
    else {mx_my_str = sfile;}
    

    ostringstream fname;
    fname << "scores/" << mx_my_str << ".txt";
    txtfile = fname.str();

    ostringstream filename;
    filename << base << "NMSSM_presel/NMSSM_XYH_YToHH_6b_" << mx_my_str << "/ntuple.root";
    infile = filename.str();

    ostringstream outname;
    outname << "gnn/" << mx_my_str << ".root";
    outfile = outname.str();
  }

  if (is_data) {
    txtfile = "scores/data.txt";
    outfile = "gnn/data.root";
    infile = "root://cmseos.fnal.gov//store/user/srosenzw/sixb/ntuples/Summer2018UL/presel/JetHT_Data_UL/JetHT_Run2018_full/ntuple.root";
  }

  info("scores read from", txtfile);
  info("preselections read from", infile);
  info("gnn selections saved to", outfile);

  fstream myfile(txtfile, std::ios_base::in);

  TFile outputFile(outfile.c_str(), "recreate");
  // TFile *fout = new TFile(outfile,"RECREATE");
  TTree *t1   = new TTree("sixBtree","sixBtree");
  TString tree = "sixBtree";
  TChain *ch  = new TChain(tree);
  ch->AddFile(Form(infile.c_str()));

  TTreeReader reader(ch);
  TTreeReaderValue<int> n_jet(reader,"n_jet");
  TTreeReaderArray<float> jet_pt(reader,"jet_pt");
  TTreeReaderArray<float> jet_eta(reader,"jet_eta");
  TTreeReaderArray<float> jet_phi(reader,"jet_phi");
  TTreeReaderArray<float> jet_m(reader,"jet_m");
  TTreeReaderArray<float> jet_btag(reader,"jet_btag");
  TTreeReaderArray<int> jet_signalId(reader,"jet_signalId");
  // TTreeReaderArray<float> jet_ptRegressed(reader,"jet_ptRegressed");
  // TTreeReaderArray<float> jet_mRegressed(reader,"jet_mRegressed");
  
  info("declaring vars");
  int n_jets, itriH;
  bool vsr_mask, vcr_mask, asr_mask, acr_mask;
  float score, pt6bsum, dR6bmin, jets_avg_btag, DeltaM_AR, DeltaM_VR, HX_pt, HY1_pt, HY2_pt, HX_m, HY1_m, HY2_m,  HX_dR, HY1_dR, HY2_dR, HY1_HY2_dEta, HY2_HX_dEta, HX_HY1_dEta, HX_HY1_dPhi, HY1_HY2_dPhi, HY2_HX_dPhi, HY1_HY2_dR, HY2_HX_dR, HX_HY1_dR, HY1_costheta, HY2_costheta, HX_costheta, X_m;
  vector<int> jets_signalId;
  vector<float> scores;
  vector<float> jets_pt, jets_eta, jets_phi, jets_btag, jets_m, jets_mRegressed, jets_ptRegressed;
  vector<float> dijets_pt, dijets_eta, dijets_phi, dijets_m, dijets_dR, temp_dR;
  vector<pair<float,int>> triH_d_hhh_ar, triH_d_hhh_vr;
  vector<p4_t> jets_p4, tri_dijet_sys, H_p4, dijets;
  p4_t jet_p4;
  
  info("reading branches");
  t1->Branch("n_jet", &n_jets);
  t1->Branch("jet_pt", &jets_pt);
  t1->Branch("jet_eta", &jets_eta);
  t1->Branch("jet_phi", &jets_phi);
  t1->Branch("jet_m", &jets_m);
  // t1->Branch("jet_mRegressed", &jets_mRegressed);
  t1->Branch("jet_btag", &jets_btag);
  t1->Branch("jet_avg_btag", &jets_avg_btag);
  t1->Branch("jet_signalId", &jets_signalId);
  t1->Branch("jet_gnn_score", &scores);
  t1->Branch("dijet_pt", &dijets_pt);
  t1->Branch("dijet_eta", &dijets_eta);
  t1->Branch("dijet_phi", &dijets_phi);
  t1->Branch("dijet_m", &dijets_m);
  t1->Branch("dijet_dR", &dijets_dR);
  t1->Branch("DeltaM_AR", &DeltaM_AR);
  t1->Branch("DeltaM_VR", &DeltaM_VR);
  t1->Branch("pt6bsum", &pt6bsum);
  t1->Branch("dR6bmin", &dR6bmin);
  t1->Branch("HX_pt", &HX_pt);
  t1->Branch("HY1_pt", &HY1_pt);
  t1->Branch("HY2_pt", &HY2_pt);
  t1->Branch("HX_m", &HX_m);
  t1->Branch("HY1_m", &HY1_m);
  t1->Branch("HY2_m", &HY2_m);
  t1->Branch("HX_dR", &HX_dR);
  t1->Branch("HY1_dR", &HY1_dR);
  t1->Branch("HY2_dR", &HY2_dR);
  t1->Branch("HY1_HY2_dEta", &HY1_HY2_dEta);
  t1->Branch("HY2_HX_dEta", &HY2_HX_dEta);
  t1->Branch("HX_HY1_dEta", &HX_HY1_dEta);
  t1->Branch("HX_HY1_dPhi", &HX_HY1_dPhi);
  t1->Branch("HY1_HY2_dPhi", &HY1_HY2_dPhi);
  t1->Branch("HY2_HX_dPhi", &HY1_HY2_dPhi);
  t1->Branch("HY1_HY2_dR", &HY1_HY2_dR);
  t1->Branch("HY2_HX_dR", &HY2_HX_dR);
  t1->Branch("HX_HY1_dR", &HX_HY1_dR);
  t1->Branch("HY1_costheta", &HY1_costheta);
  t1->Branch("HY2_costheta", &HY2_costheta);
  t1->Branch("HX_costheta", &HX_costheta);
  t1->Branch("X_m", &X_m);
  t1->Branch("vsr_mask", &vsr_mask);
  t1->Branch("vcr_mask", &vcr_mask);
  t1->Branch("asr_mask", &asr_mask);
  t1->Branch("acr_mask", &acr_mask);
  // t1->Branch("jet_ptRegressed",& jets_ptRegressed);
  // t1->Branch("jet_mRegressed",& jets_mRegressed);
  
  
  const std::vector<std::vector<int>> dijet_pairings = {
      {0, 1},{0, 2},{0, 3},{0, 4},{0, 5},
      {1, 2},{1, 3},{1, 4},{1, 5},
      {2, 3},{2, 4},{2, 5},
      {3, 4},{3, 5},
      {4, 5}
    };

  const std::vector<std::vector<int>> triH_pairings = {
      {0,  9, 14},
      {0, 10, 13},
      {0, 11, 12},
      {1,  6, 14},
      {1,  7, 13},
      {1,  8, 12},
      {2,  5, 14},
      {2,  7, 11},
      {2,  8, 10},
      {3,  5, 13},
      {3,  6, 11},
      {3,  8,  9},
      {4,  5, 12},
      {4,  6, 10},
      {4,  7,  9}
    };

  
  int a_center = 125;
  int v_center = 125 + (30 + 40)/sqrt(2);

  int eventCount = 0;
  int numEvents = 0;
  info("starting loop");
  while(reader.Next()){
    eventCount++;
    if (eventCount % 100000 == 0) {
      std::cout << eventCount << " events read!" << std::endl;
    }
    // if (eventCount == 10000) break;

    scores.clear();

    // jets_ptRegressed.clear();
    // jets_mRegressed.clear();
    jets_pt.clear();
    jets_ptRegressed.clear();
    jets_eta.clear();
    jets_phi.clear();
    jets_m.clear();
    jets_mRegressed.clear();
    jets_btag.clear();
    jets_signalId.clear();
    jets_p4.clear();

    triH_d_hhh_ar.clear();
    triH_d_hhh_vr.clear();
    
    dijets.clear();
    dijets_pt.clear();
    dijets_eta.clear();
    dijets_phi.clear();
    dijets_m.clear();
    dijets_dR.clear();
    temp_dR.clear();

    H_p4.clear();

    vsr_mask = false;
    vcr_mask = false;
    asr_mask = false;
    acr_mask = false;
    
    n_jets = *n_jet;

    // Read scores from txt file
    for (int i=0; i<n_jets; i++) {
      myfile >> score;
      scores.emplace_back(score);
    }

    vector<size_t> idx(scores.size());
    iota(idx.begin(), idx.end(), 0);

    stable_sort(idx.begin(), idx.end(),
       [&scores](size_t i1, size_t i2) {return scores[i1] > scores[i2];});

    jets_avg_btag = 0;
    pt6bsum = 0;
    for (int i=0; i<6; i++) {
      // jets_ptRegressed.emplace_back(jet_ptRegressed[idx[i]]);
      // jets_mRegressed.emplace_back(jet_mRegressed[idx[i]]);
      jets_pt.emplace_back(jet_pt[idx[i]]);
      jets_eta.emplace_back(jet_eta[idx[i]]);
      jets_phi.emplace_back(jet_phi[idx[i]]);
      jets_m.emplace_back(jet_m[idx[i]]);
      jets_btag.emplace_back(jet_btag[idx[i]]);
      jets_avg_btag += jets_btag[i];

      jets_signalId.emplace_back(jet_signalId[idx[i]]);

      jet_p4.SetCoordinates (jets_pt[i], jets_eta[i], jets_phi[i], jets_m[i]);
      jets_p4.emplace_back(jet_p4);
      
      pt6bsum += jets_pt[i];
    }
    jets_avg_btag = jets_avg_btag / 6;

    dR6bmin = 999;
    for (const std::vector<int> ijets : dijet_pairings)
    { 
      int ij1 = ijets[0]; int ij2 = ijets[1];

      p4_t dijet_p4 = jets_p4[ij1] + jets_p4[ij2];
      dijets.emplace_back(dijet_p4);

      double dijet_dR = ROOT::Math::VectorUtil::DeltaR(jets_p4[ij1], jets_p4[ij2]);
      temp_dR.emplace_back(dijet_dR);

      if (dijet_dR < dR6bmin) {
        dR6bmin = dijet_dR;
      }
    }
    
    for (unsigned int i = 0; i < triH_pairings.size(); i++) {
      tri_dijet_sys.clear();
      for (int id : triH_pairings[i]) tri_dijet_sys.push_back( dijets[id] );

      std::sort(tri_dijet_sys.begin(),tri_dijet_sys.end(),[](p4_t dj1,p4_t dj2){ return dj1.Pt()>dj2.Pt(); });

      float r_AR = sqrt(pow(tri_dijet_sys[0].M()-a_center,2) + pow(tri_dijet_sys[1].M()-a_center,2) + pow(tri_dijet_sys[2].M()-a_center,2));
      float r_VR = sqrt(pow(tri_dijet_sys[0].M()-v_center,2) + pow(tri_dijet_sys[1].M()-v_center,2) + pow(tri_dijet_sys[2].M()-v_center,2));
      
      triH_d_hhh_ar.push_back( std::make_pair(r_AR,i) );
      triH_d_hhh_vr.push_back( std::make_pair(r_VR,i) );
    }
    
    std::sort(triH_d_hhh_ar.begin(),triH_d_hhh_ar.end(),[](std::pair<float,int> h1,std::pair<float,int> h2){ return h1.first<h2.first; });
    DeltaM_AR = triH_d_hhh_ar[0].first;
    
    std::sort(triH_d_hhh_vr.begin(),triH_d_hhh_vr.end(),[](std::pair<float,int> h1,std::pair<float,int> h2){ return h1.first<h2.first; });
    DeltaM_VR = triH_d_hhh_vr[0].first;

    if (DeltaM_AR > 40 && DeltaM_VR > 40) {continue;}
    
    if (DeltaM_AR < 40) {
      acr_mask = true;
      if (DeltaM_AR < 30) {
        asr_mask = true;
        acr_mask = false;
      }
      itriH = triH_d_hhh_ar[0].second;
    }
    else if (DeltaM_VR < 40) {
      vcr_mask = true;
      if (DeltaM_VR < 30) {
        vsr_mask = true;
        vcr_mask = false;
      }
      itriH = triH_d_hhh_vr[0].second;
    }

    numEvents++;

    std::vector<int> idijets = triH_pairings[itriH];

    std::sort(idijets.begin(),idijets.end(),[dijets](int id1,int id2){ return dijets[id1].Pt() > dijets[id2].Pt(); });

    for (int i=0; i < idijets.size(); i++)
    {
      H_p4.emplace_back(dijets[idijets[i]]);
      dijets_pt.emplace_back(dijets[idijets[i]].Pt());
      dijets_eta.emplace_back(dijets[idijets[i]].Eta());
      dijets_phi.emplace_back(dijets[idijets[i]].Phi());
      dijets_m.emplace_back(dijets[idijets[i]].M());
      dijets_dR.emplace_back(temp_dR[idijets[i]]);
    }


    HX_m = dijets_m[0];
    HY1_m = dijets_m[1];
    HY2_m = dijets_m[2];

    HX_pt = dijets_pt[0];
    HY1_pt = dijets_pt[1];
    HY2_pt = dijets_pt[2];

    HX_dR = dijets_dR[0];
    HY1_dR = dijets_dR[1];
    HY2_dR = dijets_dR[2];

    HX_HY1_dEta = abs(H_p4[0].Eta() - H_p4[1].Eta());
    HY1_HY2_dEta = abs(H_p4[1].Eta() - H_p4[2].Eta());
    HY2_HX_dEta = abs(H_p4[2].Eta() - H_p4[0].Eta());
    
    HX_HY1_dPhi = H_p4[0].Phi() - H_p4[1].Phi();
    HY1_HY2_dPhi = H_p4[1].Phi() - H_p4[2].Phi();
    HY2_HX_dPhi = H_p4[2].Phi() - H_p4[0].Phi();

    if (HX_HY1_dEta < -M_PI) {HX_HY1_dEta += 2*M_PI;}
    if (HX_HY1_dEta < +M_PI) {HX_HY1_dEta -= 2*M_PI;}
    if (HY1_HY2_dPhi < -M_PI) {HY1_HY2_dPhi += 2*M_PI;}
    if (HY1_HY2_dPhi < +M_PI) {HY1_HY2_dPhi -= 2*M_PI;}
    if (HY2_HX_dPhi < -M_PI) {HY2_HX_dPhi += 2*M_PI;}
    if (HY2_HX_dPhi < +M_PI) {HY2_HX_dPhi -= 2*M_PI;}

    HX_HY1_dR = ROOT::Math::VectorUtil::DeltaR(H_p4[0], H_p4[1]);
    HY1_HY2_dR = ROOT::Math::VectorUtil::DeltaR(H_p4[1], H_p4[2]);
    HY2_HX_dR = ROOT::Math::VectorUtil::DeltaR(H_p4[2], H_p4[0]);
    
    HY1_costheta = ROOT::Math::VectorUtil::CosTheta(H_p4[0], H_p4[1]);
    HY2_costheta = ROOT::Math::VectorUtil::CosTheta(H_p4[1], H_p4[2]);
    HX_costheta = ROOT::Math::VectorUtil::CosTheta(H_p4[2], H_p4[0]);

    p4_t X = H_p4[0] + H_p4[1] + H_p4[2];
    X_m = X.M();

    t1->Fill();

  } // end event loop
  
  t1->Write();

  std::cout << endl << eventCount << " Total Events" << endl;
  std::cout << numEvents << " Events Passed" << endl;

  return 0;
} // end function

