// skim_ntuple.exe --input input/PrivateMC_2018/NMSSM_XYH_YToHH_6b_MX_600_MY_400.txt --cfg config/skim_ntuple_2018.cfg  --output prova.root --is-signal
// skim_ntuple.exe --input input/Run2_UL/2018/TTJets.txt --cfg config/skim_ntuple_2018_ttbar.cfg  --output prova_ttbar.root
// skim_ntuple.exe --input input/Run2_UL/2018/SingleMuon_Run2.txt --cfg config/skim_ntuple_2018_ttbar.cfg  --output prova_singlemu_ttbarskim.root --is-data

#include <iostream>
#include <string>
#include <iomanip>
#include <any>
#include <chrono>

#include <boost/program_options.hpp>
namespace po = boost::program_options;

#include "DirectionalCut.h"
#include "CfgParser.h"
#include "NanoAODTree.h"
#include "NormWeightTree.h"
#include "SkimUtils.h"
namespace su = SkimUtils;

#include "OutputTree.h"
#include "jsonLumiFilter.h"

#include "TrgEff_functions.h"
#include "Skim_functions.h"
#include "SixB_functions.h"
#include "EightB_functions.h"
#include "TTBar_functions.h"
#include "JetTools.h"
#include "BtagSF.h"
#include "EventShapeCalculator.h"
#include "Cutflow.h"
#include "EvalNN.h"

#include "Timer.h"
#include "DebugUtils.h"

#include "TFile.h"
#include "TROOT.h"
#include "TH1F.h"

using namespace std;

std::vector<std::string> split_by_delimiter(std::string input, std::string delimiter)
{
  std::vector<std::string> tokens;
  if(input == "")
    return tokens;

  size_t pos = 0;
  while ((pos = input.find(delimiter)) != std::string::npos)
    {
      tokens.push_back(input.substr(0, pos));
      input.erase(0, pos + delimiter.length());
    }
  tokens.push_back(input); // last part splitted

  return tokens;
}

Variation string_to_jer_variation (std::string s)
{
  if (s == "nominal")
    return Variation::NOMINAL;
  if (s == "up")
    return Variation::UP;
  if (s == "down")
    return Variation::DOWN;
  throw std::runtime_error(string("Cannot parse the variation ") + s);
}



// -----------------------------------
// shortcuts to access cfgParser options with defaults

template <typename T>
T readCfgOptWithDefault(CfgParser& config, std::string optName, T default_ret){
  throw std::runtime_error("please provide an override template implementation for this type");
}

template <>
bool readCfgOptWithDefault<bool>(CfgParser& config, std::string optName, bool default_ret){
  if (config.hasOpt(optName))
    return config.readBoolOpt(optName);
  else
    return default_ret;
}

template <>
int readCfgOptWithDefault<int>(CfgParser& config, std::string optName, int default_ret){
  if (config.hasOpt(optName))
    return config.readIntOpt(optName);
  else
    return default_ret;
}

template <>
float readCfgOptWithDefault<float>(CfgParser& config, std::string optName, float default_ret){
  if (config.hasOpt(optName))
    return config.readFloatOpt(optName);
  else
    return default_ret;
}

template <>
double readCfgOptWithDefault<double>(CfgParser& config, std::string optName, double default_ret){
  if (config.hasOpt(optName))
    return config.readDoubleOpt(optName);
  else
    return default_ret;
}

template <>
std::string readCfgOptWithDefault<std::string>(CfgParser& config, std::string optName, std::string default_ret){
  if (config.hasOpt(optName))
    return config.readStringOpt(optName);
  else
    return default_ret;
}

// -----------------------------------


int main(int argc, char** argv)
{
  cout << "[INFO] ... starting program" << endl;
  const auto start_prog_t = chrono::high_resolution_clock::now();

  ////////////////////////////////////////////////////////////////////////
  // Declare command line options
  ////////////////////////////////////////////////////////////////////////
    
  po::options_description desc("Skim options");
  desc.add_options()
    ("help", "produce help message")
    // required
    ("cfg"   , po::value<string>()->required(), "skim config")
    ("input" , po::value<string>()->required(), "input file list")
    ("output", po::value<string>()->required(), "output file LFN")
    // optional
    // ("xs"        , po::value<float>(), "cross section [pb]")
    ("maxEvts"   , po::value<int>()->default_value(-1), "max number of events to process")
    ("pickEvt"   , po::value<string>()->default_value(""), "run on this run:lumi:event number only (for debug). Use wildcard * to match all")
    ("puWeight"  , po::value<string>()->default_value(""), "PU weight file name")
    ("seed"      , po::value<int>()->default_value(12345), "seed to be used in systematic uncertainties such as JEC, JER, etc")
    // ("kl-rew-list"  , po::value<std::vector<float>>()->multitoken()->default_value(std::vector<float>(0), "-"), "list of klambda values for reweight")
    // ("kl-rew"    , po::value<float>(),  "klambda value for reweighting")
    // ("kl-map"    , po::value<string>()->default_value(""), "klambda input map for reweighting")
    // ("kl-histo"  , po::value<string>()->default_value("hhGenLevelDistr"), "klambda histogram name for reweighting")
    ("jes-shift-syst",  po::value<string>()->default_value("nominal"), "Name of the JES (scale) source uncertainty to be shifted. Usage as <name>:<up/down>. Pass -nominal- to not shift the jets")
    ("jer-shift-syst",  po::value<string>()->default_value("nominal"), "Name of the JER (resolution) source uncertainty to be shifted. Usage as <up/down>. Pass -nominal- to not shift the jets")
    ("bjer-shift-syst", po::value<string>()->default_value("nominal"), "Name of the b regressed JER (resolution) source uncertainty to be shifted. Usage as <up/down>. Pass -nominal- to not shift the jets")
    // pairing variables
    // ("bbbbChoice"    , po::value<string>()->default_value("BothClosestToDiagonal"), "bbbb pairing choice")
    // ("mh1mh2"        , po::value<float>()->default_value(1.05), "Ratio Xo/Yo or 1/slope of the diagonal") 
    // ("option"        , po::value<int>()->default_value(0), "Option: 0=Nominal, 1=Alternative 1, 2=Alternative 2") 
    // flags
    ("is-data",       po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a data sample (default is false)")
    ("is-signal",     po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "mark as a HH signal sample (default is false)")
    //
    ("save-p4",       po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "save the tlorentzvectors in the output")
    //
    ("no-genw-tree",  po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "disable the storage of the genweight tree for normalizations")
    //
    ("debug",         po::value<bool>()->zero_tokens()->implicit_value(true)->default_value(false), "debug this event (verbose printing)")
    ;

  po::variables_map opts;
  try {
    po::store(parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), opts);
    if (opts.count("help")) {
      cout << desc << "\n";
      return 1;
    }
    po::notify(opts);
  }    
  catch (po::error& e) {
    cerr << "** [ERROR] " << e.what() << endl;
    return 1;
  }

  ////////////////////////////////////////////////////////////////////////
  // Read config and other cmd line options for skims
  ////////////////////////////////////////////////////////////////////////

  const bool is_data = opts["is-data"].as<bool>();
  cout << "[INFO] ... is a data sample? " << std::boolalpha << is_data << std::noboolalpha << endl;

  const bool is_signal = (is_data ? false : opts["is-signal"].as<bool>());
  cout << "[INFO] ... is a signal sample? " << std::boolalpha << is_signal << std::noboolalpha << endl;

  const bool save_genw_tree = (is_data ? false : !opts["no-genw-tree"].as<bool>());
  cout << "[INFO] ... will save the gen weight NormTree ? " << std::boolalpha << save_genw_tree << std::noboolalpha << endl;

  CfgParser config;
  if (!config.init(opts["cfg"].as<string>())){
    cerr << "** [ERROR] no config file was provided" << endl;
    return 1;
  }
  cout << "[INFO] ... using config file " << opts["cfg"].as<string>() << endl;

  enum SkimTypes
    {
      kTrgEff,
      ksixb,
      kttbar,
      keightb,
      // kshapecr,
      khiggscr,
      kpass,
      kpresel,
      knull
    };
  
  string year = config.readStringOpt("parameters::year");
  string skim_type_name = config.readStringOpt("configurations::skimType");
  cout << "[INFO] ... skim type " << skim_type_name << endl;
  const SkimTypes skim_type = (skim_type_name == "sixb"    ? ksixb    :
                               skim_type_name == "ttbar"   ? kttbar   :
                               skim_type_name == "eightb"  ? keightb  :
			                        //  skim_type_name == "shapecr" ? kshapecr :
                               skim_type_name == "pass"    ? kpass    :
                               skim_type_name == "presel"  ? kpresel  :
                               skim_type_name == "trgeff"  ? kTrgEff  :
                               knull
                              );
  if (skim_type == knull)
    throw std::runtime_error("skim type not recognized");

  ////////////////////////////////////////////////////////////////////////
  // Prepare event loop
  ////////////////////////////////////////////////////////////////////////

  cout << "[INFO] ... opening file list : " << opts["input"].as<string>().c_str() << endl;
  if ( access( opts["input"].as<string>().c_str(), F_OK ) == -1 ){
    cerr << "** [ERROR] The input file list does not exist, aborting" << endl;
    return 1;        
  }

  // Joining all the NANOAOD input file in a TChain in order to be used like an unique three
  TChain ch("Events");
  int nfiles = su::appendFromFileList(&ch, opts["input"].as<string>());
    
  if (nfiles == 0){
    cerr << "** [ERROR] The input file list contains no files, aborting" << endl;
    return 1;
  }
  cout << "[INFO] ... file list contains " << nfiles << " files" << endl;

  cout << "[INFO] ... creating tree reader" << endl;

  // The TChain is passed to the NanoAODTree_SetBranchImpl to parse all the branches
  NanoAODTree nat (&ch);

  ////////////////////////////////////////////////////////////////////////
  // Trigger information
  ////////////////////////////////////////////////////////////////////////
  
  cout << "[INFO] ... loading " << config.readStringListOpt("triggers::makeORof").size() << " triggers" << endl;
  
  const bool apply_trigger     = config.readBoolOpt("triggers::applyTrigger");
  const bool save_trg_decision = config.readBoolOpt("triggers::saveDecision");
  cout << "[INFO] ... is the OR decision of these triggers applied? " << std::boolalpha << apply_trigger << std::noboolalpha << endl;
  cout << "[INFO] ... will save the trigger decision? " << std::boolalpha << save_trg_decision << std::noboolalpha << endl;

  // std::vector<std::string> triggerAndNameVector;
  // if(apply_trigger) triggerAndNameVector = config.readStringListOpt("triggers::makeORof");
  // <triggerName , < objectBit, minNumber> >
  // std::map<std::string, std::map< std::pair<int,int>, int > > triggerObjectAndMinNumberMap;

  std::vector<std::string> triggerAndNameVector = config.readStringListOpt("triggers::makeORof");
  std::vector<std::string> triggerVector;

  cout << "[INFO] ... listing the triggers applied" << endl;
  for (auto & trigger : triggerAndNameVector)
    {
      if(trigger == "")
        continue;
        
      std::vector<std::string> triggerTokens = split_by_delimiter(trigger, ":");
      if (triggerTokens.size() != 2)
        throw std::runtime_error("** skim_ntuple : could not parse trigger entry " + trigger + " , aborting");

      triggerVector.push_back(triggerTokens[1]);
      cout << "   - " << triggerTokens[0] << "  ==> " << triggerTokens[1] << endl;

      // if(!config.hasOpt( Form("triggers::%s_ObjectRequirements",triggerTokens[0].data()) ))
      // {
      //     cout<<Form("triggers::%s_ObjectRequirements",triggerTokens[0].data())<<std::endl;
      //     cout<<"Trigger "<< triggerTokens[1] <<" does not have ObjectRequirements are not defined";
      //     continue;
      // }

      // triggerObjectAndMinNumberMap[triggerTokens[1]] = std::map< std::pair<int,int>, int>();   

      // std::vector<std::string> triggerObjectMatchingVector = config.readStringListOpt(Form("triggers::%s_ObjectRequirements",triggerTokens[0].data()));

      // for (auto & triggerObject : triggerObjectMatchingVector)
      // {

      //     std::vector<std::string> triggerObjectTokens;
      //     while ((pos = triggerObject.find(delimiter)) != std::string::npos)
      //     {
      //         triggerObjectTokens.push_back(triggerObject.substr(0, pos));
      //         triggerObject.erase(0, pos + delimiter.length());
      //     }
      //     triggerObjectTokens.push_back(triggerObject); // last part splitted
      //     if (triggerObjectTokens.size() != 3)
      //     {
      //         throw std::runtime_error("** skim_ntuple : could not parse trigger entry " + triggerObject + " , aborting");
      //     }

      //     triggerObjectAndMinNumberMap[triggerTokens[1]][std::pair<int,int>(atoi(triggerObjectTokens[0].data()),atoi(triggerObjectTokens[1].data()))] = atoi(triggerObjectTokens[2].data());
      // }
    }

  nat.triggerReader().setTriggers(triggerVector);

  ////////////////////////////////////////////////////////////////////////
  // Prepare the output
  ////////////////////////////////////////////////////////////////////////
 
  string outputFileName = opts["output"].as<string>();
  cout << "[INFO] ... saving output to file : " << outputFileName << endl;
  TFile outputFile(outputFileName.c_str(), "recreate");
  OutputTree ot(opts["save-p4"].as<bool>(),
                map<string, bool>{
          // {"jet_coll",   config.readBoolOpt("configurations::saveJetColl")},
          // {"shape_brs",  config.readBoolOpt("configurations::saveShapes")},
          {"muon_coll",   readCfgOptWithDefault<bool>(config, "configurations::saveMuonColl", false)},
          {"ele_coll",    readCfgOptWithDefault<bool>(config, "configurations::saveEleColl", false)},  
          {"jet_coll",    readCfgOptWithDefault<bool>(config, "configurations::saveJetColl", false)},
          {"shape_brs",   readCfgOptWithDefault<bool>(config, "configurations::saveShapes", false)},
	        {"dijets_coll", readCfgOptWithDefault<bool>(config, "configurations::saveDiJets", false)},
	        {"sixb_brs", (skim_type == ksixb)},
	        {"eightb_brs", (skim_type == keightb)},
	        {"ttbar_brs", (skim_type == kttbar)},
		      {"trgeff_brs", (skim_type == kTrgEff)},
	        {"sig_gen_brs", (is_signal)},
	        {"gen_brs", (!is_data)},
		  });
  
  if (save_trg_decision) {
    for (auto& tname : triggerVector)
      ot.declareUserIntBranch(tname,   0);
  }

  std::string pu_weight_file;
  if (!is_data){
    if (opts["puWeight"].as<string>().size() != 0){ // a valid option is passed from cmd line
      cout << "[INFO] Using custom PU weight file passed from cmd line options" << endl;
      pu_weight_file = opts["puWeight"].as<string>();
    }
    else { // revert to default in skim cfg
      cout << "[INFO] Using PU weight file from skim cfg" << endl;
      // pu_weight_file = config.readStringOpt("parameters::PUweightFile");
      pu_weight_file = readCfgOptWithDefault<string>(config, "parameters::PUweightFile", "");
    }
  }

  NormWeightTree nwt;
  map<string, string> pu_data{
      {"filename", pu_weight_file},
      {"name_PU_w", "PUweights"},
      {"name_PU_w_up", "PUweights_up"},
      {"name_PU_w_do", "PUweights_down"},
  };

  // // just a test
  // nwt.add_weight("test1", {"test1_up", "test1_down"});
  // nwt.add_weight("test2", {"test2_A", "test2_B", "test2_C"});
  // nwt.add_weight("test3", {});

  ////////////////////////////////////////////////////////////////////////
  // All pre-running configurations (corrections, methods from cfg, etc)
  ////////////////////////////////////////////////////////////////////////

  jsonLumiFilter jlf;
  if (is_data)
    jlf.loadJSON(config.readStringOpt("data::lumimask")); // just read the info for data, so if I just skim MC I'm not forced to parse a JSON

  // -----------

  Timer loop_timer;

  // -----------

  Skim_functions* skf;

  switch (skim_type)
  {
  case kTrgEff:
    skf = new TrgEff_functions();
    break;
  case ksixb:
    skf = new SixB_functions();
    break;
  case keightb:
    skf = new EightB_functions();
    break;
  case kttbar:
    skf = new TTBar_functions();
    break;
  case kpresel:
    skf = new SixB_functions();
    break;
  default:
    skf = new Skim_functions();
    break;
  }
  skf->Print();
  skf->set_timer(&loop_timer);
  // skf->set_debug(true);
  // -----------
    
  const std::vector<double> btag_WPs = config.readDoubleListOpt("configurations::bTagWPDef");
  const int nMinBtag = config.readIntOpt("configurations::nMinBtag");
  const int bTagWP   = config.readIntOpt("configurations::bTagWP");

  skf->set_btag_WPs(config.readDoubleListOpt("configurations::bTagWPDef"));

  cout << "[INFO] ... events must contain >= " << nMinBtag << " jets passing WP (0:L, 1:M, 2:T) : " << bTagWP << endl;
  cout << "[INFO] ... the WPs are: (L/M/T) : " << btag_WPs.at(0) << "/" << btag_WPs.at(1) << "/" << btag_WPs.at(2) << endl;

  BtagSF btsf;
  if (!is_data)
    {
      string btsffile = config.readStringOpt("parameters::DeepJetScaleFactorFile");
      btsf.init_reader("DeepJet", btsffile);
      btsf.set_WPs(btag_WPs.at(0), btag_WPs.at(1), btag_WPs.at(2));
    }
  
  bool blind = false;
  if (config.hasOpt("configurations::blinded"))
    blind = config.readBoolOpt("configurations::blinded");

  // -----------

  JetTools jt;

  string jes_shift = opts["jes-shift-syst"].as<string>();
  bool do_jes_shift = (jes_shift != "nominal");
  cout << "[INFO] ... shifting jet energy scale? " << std::boolalpha << do_jes_shift << std::noboolalpha << endl;
  bool dir_jes_shift_is_up;
  if (do_jes_shift && !is_data){
    string JECFileName = config.readStringOpt("parameters::JECFileName");
    auto tokens = split_by_delimiter(opts["jes-shift-syst"].as<string>(), ":");
    if (tokens.size() != 2)
      throw std::runtime_error(string("Cannot parse the jes shift name : ") + opts["jes-shift-syst"].as<string>());
    string jes_syst_name = tokens.at(0);
    dir_jes_shift_is_up   = (tokens.at(1) == "up"   ? true  :
                             tokens.at(1) == "down" ? false :
                             throw std::runtime_error(string("Could not parse jes direction ") + tokens.at(1)));
    cout << "       ... jec file name           : " << JECFileName << endl;
    cout << "       ... jet energy scale syst   : " << jes_syst_name << endl;
    cout << "       ... jet energy scale is up? : " << std::boolalpha << dir_jes_shift_is_up << std::noboolalpha << endl;
    jt.init_jec_shift(JECFileName, jes_syst_name);
  }

  // FIXME: block below to be run only if !is_data?
  string JERScaleFactorFile = config.readStringOpt("parameters::JERScaleFactorFile");
  string JERResolutionFile  = config.readStringOpt("parameters::JERResolutionFile");
  const int rndm_seed = opts["seed"].as<int>();
  cout << "[INFO] ... initialising JER corrector with the following parameters" << endl;
  cout << "       ... SF file         : " << JERScaleFactorFile << endl;
  cout << "       ... resolution file : " << JERResolutionFile << endl;
  cout << "       ... rndm seed       : " << rndm_seed << endl;
  jt.init_smear(JERScaleFactorFile, JERResolutionFile, rndm_seed);

  cout << "[INFO] ... jet resolution syst is    : " << opts["jer-shift-syst"].as<string>() << endl;
  cout << "[INFO] ... b regr resolution syst is : " << opts["bjer-shift-syst"].as<string>() << endl;
  const Variation jer_var  = string_to_jer_variation(opts["jer-shift-syst"].as<string>());
  const Variation bjer_var = string_to_jer_variation(opts["bjer-shift-syst"].as<string>());

  // ------------------------------------------------------------------
  // FIXME: move to another file to make code slimmer here?
  // ----------- configure the sixB functions
  cout << "[INFO] ... configurations read from the config file" << endl;

  skf->initialize_params_from_cfg(config);
  skf->initialize_functions(outputFile);

  ////////////////////////////////////////////////////////////////////////
  // Execute event loop
  ////////////////////////////////////////////////////////////////////////

  // max nEv debug
  const int maxEvts = opts["maxEvts"].as<int>();
  if (maxEvts >= 0)
    cout << "[INFO] ... running on : " << maxEvts << " events" << endl;

  // single Ev debug
  std::vector<std::string> RunLumiEvtStr =  split_by_delimiter(opts["pickEvt"].as<string>(), ":");
  const int pickRunNr       = (RunLumiEvtStr.size() != 3 ? -1 : RunLumiEvtStr.at(0) == "*" ? -1 : std::stoi (RunLumiEvtStr.at(0)) );
  const int pickLumiNr      = (RunLumiEvtStr.size() != 3 ? -1 : RunLumiEvtStr.at(1) == "*" ? -1 : std::stoi (RunLumiEvtStr.at(1)) );
  const long long pickEvtNr = (RunLumiEvtStr.size() != 3 ? -1 : RunLumiEvtStr.at(2) == "*" ? -1 : std::stoll(RunLumiEvtStr.at(2)) );
  const bool doPickEvt = (pickEvtNr >= 0 || pickRunNr >= 0 || pickLumiNr >= 0);
  if (RunLumiEvtStr.size() > 0 )
    cout << "[INFO] ... running on a single run:lumi:event for debug: " << opts["pickEvt"].as<string>() << " --> run:lumi:evt " << pickRunNr << ":" << pickLumiNr << ":" << pickEvtNr << endl;
  // if (RunLumiEvtStr.size() > 0 && RunLumiEvtStr.size() != 3){
  //   cout << "[ERROR] ... pickEvt string must be formatted as evt:run:lumi " << endl;
  // }

  const bool debug = opts["debug"].as<bool>();
  if (debug) skf->set_debug(debug);

  const auto start_loop_t = chrono::high_resolution_clock::now();

  Cutflow cutflow;
  Cutflow cutflow_Unweighted("h_cutflow_unweighted", "Unweighted selection cutflow");
  
  for (int iEv = 0; true; ++iEv)
  {
    if (maxEvts >= 0 && iEv >= maxEvts)
      break;

    loop_timer.start_lap();

    if (!nat.Next()) break;
    if (iEv % 10000 == 0 || debug) {
      cout << "... processing event " << iEv << endl;
      // auto bsize  = ot.getTree()->GetBranch("Run")->GetBasketSize();
      // cout << "... tree basket size (branch Run) : " << bsize  << endl;
    }

    // use the tree content to initialise weight tree in the first event
    if (iEv == 0 && !is_data && save_genw_tree){
      nwt.init_weights(nat, pu_data); // get the syst structure from nanoAOD
      su::init_gen_weights(ot, nwt);  // and forward it to the output tree
    }

    // pick a specific event for debug
    if ( doPickEvt && (
          (pickEvtNr >= 0 && (long long) *(nat.event) != pickEvtNr) ||
          (pickRunNr >= 0 && (int) *(nat.run) != pickRunNr) ||
          (pickLumiNr >= 0 && (long long) *(nat.luminosityBlock) != pickLumiNr) ) )
      continue;

    // apply certification json file
    if (is_data && !jlf.isValid(*nat.run, *nat.luminosityBlock)){
      continue; // not a valid lumi
    }

    EventInfo ei;
    ot.clear();
    loop_timer.click("Input read");
    
    // save the normalization weights.
    // IMPORTANT NOTE: this must occurr *before* any MC event is filtered because of analysis selections to correctly compute the normalization
    if (!is_data && save_genw_tree){
      nwt.read_weights(nat);
      // example to fill user weights
      // auto& w1 = nwt.get_weight("test1");
      // auto& w2 = nwt.get_weight("test2");
      // auto& w3 = nwt.get_weight("test3");
      // w1.w = iEv;
      // w2.w = 10*iEv;
      // w3.w = 100*iEv;
      // w1.syst_val = {iEv + 1., iEv - 1.};
      // w2.syst_val = {10. * iEv - 10, 10. * iEv - 20, 10. * iEv - 30};
      // w3.syst_val = {};
      nwt.fill();
      loop_timer.click("Norm weight read + fill");
    }
    // ------- events can start be filtered from here (after saving all gen weights)
    cutflow.add("total", nwt);
    cutflow_Unweighted.add("total");
        
    // trigger requirements
    if (apply_trigger && !(nat.getTrgOr()) )
      continue;
    cutflow.add("trigger", nwt);
    cutflow_Unweighted.add("trigger");
  
    if (save_trg_decision) {
      auto listOfPassedTriggers = nat.getTrgPassed();
      for (auto& t : listOfPassedTriggers)
        ot.userInt(t) = 1; // defaults are left to 0
    }
    loop_timer.click("Trigger");

    // global event info
    skf->copy_event_info(nat, ei, !is_data);
    loop_timer.click("Global info");

    //==================================
    // Apply METFilters
    //==================================
    bool bMETFilters = *nat.Flag_goodVertices && *nat.Flag_globalSuperTightHalo2016Filter && *nat.Flag_HBHENoiseFilter && *nat.Flag_HBHENoiseIsoFilter && *nat.Flag_EcalDeadCellTriggerPrimitiveFilter && *nat.Flag_BadPFMuonFilter && *nat.Flag_eeBadScFilter && (*nat.Flag_ecalBadCalibFilter || (year=="2016"));
    // if (!bMETFilters) continue;
    loop_timer.click("MET Filters");
    // cutflow.add("met filters", nwt);
    // cutflow_Unweighted.add("met filters");
    
    //==================================
    // Apply muon selection or veto
    //==================================
    std::vector<Muon> selected_muons = skf->select_muons(config, nat, ei);
    ei.n_muon    = selected_muons.size();
    ei.muon_list = selected_muons;
        
    bool applyMuonVeto = config.readBoolOpt("configurations::applyMuonVeto");
    bool applyMuonSelection = config.readBoolOpt("configurations::applyMuonSelection");
    if (applyMuonVeto)
      {
	if (selected_muons.size() != 0) continue;
	cutflow.add("#mu veto", nwt);
	cutflow_Unweighted.add("#mu veto");
	loop_timer.click("#mu veto");
      }
    else if (applyMuonSelection)
      {
	const DirectionalCut<int> cfg_nMuons(config, "configurations::nMuonsCut");
	if (!cfg_nMuons.passedCut(selected_muons.size())) continue;
	cutflow.add("#mu selection", nwt);
	cutflow_Unweighted.add("#mu selection");
	loop_timer.click("#mu selection");
      }
    
    //=====================================
    // Apply electron selection or veto
    //=====================================
    std::vector<Electron> selected_electrons = skf->select_electrons(config, nat, ei);
    ei.n_ele    = selected_electrons.size();
    ei.ele_list = selected_electrons;
    
    bool applyEleVeto = config.readBoolOpt("configurations::applyEleVeto");
    bool applyEleSelection = config.readBoolOpt("configurations::applyEleSelection");
    if (applyEleVeto)
      {
	if (selected_electrons.size() !=0) continue;
	cutflow.add("e veto", nwt);
	cutflow_Unweighted.add("e veto");
	loop_timer.click("e veto");
      }
    else if (applyEleSelection)
      {
	const DirectionalCut<int> cfg_nElectrons(config, "configurations::nEleCut");
	if (!cfg_nElectrons.passedCut(selected_electrons.size())) continue;
	cutflow.add("e selection", nwt);
	cutflow_Unweighted.add("e selection");
	loop_timer.click("e selection");
      }
    
    //======================================
    // Save signal specific GEN info
    //======================================
    if (is_signal)
    {
      skf->select_gen_particles(nat, ei);        // find gen level X, Y, H, b
      skf->match_genbs_to_genjets(nat, ei);      // match the b quarks found above to the genjets
      skf->match_genbs_genjets_to_reco(nat, ei); // match the genjets found above to the reco jets
      loop_timer.click("Signal gen level");
    }
    
    //======================================
    // Jet Selection
    //======================================
    std::vector<Jet> all_jets = skf->get_all_jets(nat); // dump all nanoAOD jets into a vector<Jet>
    ei.nfound_all = skf->n_gjmatched_in_jetcoll(nat, ei, all_jets);
    ei.nfound_all_h = skf->n_ghmatched_in_jetcoll(nat, ei, all_jets);
    loop_timer.click("All jets copy");
    
    if (!is_data){
      if (do_jes_shift)
        all_jets = jt.jec_shift_jets(nat, all_jets, dir_jes_shift_is_up); // apply JEC scale shift to jets
      all_jets = jt.smear_jets(nat, all_jets, jer_var, bjer_var);         // apply JER smearing to jets
      loop_timer.click("JEC + JER");
    }
    
    // Apply preselections to jets (min pT / max eta / PU ID / PF ID)
    std::vector<Jet> presel_jets = skf->preselect_jets(nat, ei, all_jets);


    
    ei.nfound_presel = skf->n_gjmatched_in_jetcoll(nat, ei, presel_jets);
    //std::cout << "Number of selected jets found matched with GEN-level objects = "<<ei.nfound_presel<<std::endl;
    
    ei.nfound_presel_h = skf->n_ghmatched_in_jetcoll(nat, ei, presel_jets);
    //std::cout << "Number of selected jets found matched to Higgs objects = "<<ei.nfound_presel_h<<std::endl;
    
    ei.n_jet = presel_jets.size();
    //std::cout << "Number of selected jets: "<<ei.n_jet<<std::endl;
    
    loop_timer.click("Jet selection");
    if (debug) dumpObjColl(presel_jets, "==== PRESELECTED JETS ===");
    
    if (is_signal) 
    {
      skf->match_signal_recojets(nat,ei,presel_jets);
      std::vector<GenJet> all_genjets = skf->get_all_genjets(nat);
      skf->match_genjets_to_reco(nat,ei,all_genjets,presel_jets);
      skf->match_signal_genjets(nat,ei,all_genjets);
      ei.genjet_list = all_genjets;
    }
    
    ei.jet_list = presel_jets;
    
    //=================================
    // Apply analysis-specific skims
    //=================================
    if (skim_type == kTrgEff) {
      const DirectionalCut<int> cfg_nJets(config, "presel::njetsCut");
      if (!cfg_nJets.passedCut(presel_jets.size())) continue;
      cutflow.add("npresel_jets >= 6", nwt);	
      cutflow_Unweighted.add("npresel_jets >= 6");
      // Sorted jets passed cuts
      //std::vector<Jet> sorted_jets = skf->select_jets(nat, ei, presel_jets);
      
      // Work on-going here...

      }
    else if (skim_type == keightb)
    {
      if (presel_jets.size() < 8)
        continue;
      cutflow.add("npresel_jets>=8", nwt);
    
      // std::vector<DiJet> dijets = skf->make_dijets(nat, ei, presel_jets);
      // ei.dijet_list = dijets;
      
      std::vector<Jet> selected_jets = skf->select_jets(nat, ei, presel_jets);
      loop_timer.click("Eight B Selection");
      
      if (selected_jets.size() < 8)
        continue;
      cutflow.add("nselect_jets>=8", nwt);
      skf->pair_jets(nat, ei, selected_jets);
      skf->compute_seljets_btagmulti(nat, ei);
      loop_timer.click("Eight b jet pairing");
      
      if (is_signal)
      {
	if (presel_jets.size() < 8)
	  continue;
	cutflow.add("npresel_jets>=8", nwt);
	cutflow_Unweighted.add("npresel_jets>=8");
	// std::vector<DiJet> dijets = skf->make_dijets(nat, ei, presel_jets);
	// ei.dijet_list = dijets;
	
	std::vector<Jet> selected_jets = skf->select_jets(nat, ei, presel_jets);
	loop_timer.click("Eight B Selection");
	
	if (selected_jets.size() < 8)
	  continue;
	cutflow.add("nselect_jets>=8", nwt);
	cutflow_Unweighted.add("nselect_jets>=8");
	skf->pair_jets(nat, ei, selected_jets);
	skf->compute_seljets_btagmulti(nat, ei);
	loop_timer.click("Eight b jet pairing");
	
	if (is_signal)
	  {
	    skf->compute_seljets_genmatch_flags(nat, ei);
	    loop_timer.click("Eight b pairing flags");
	  }
	
	skf->compute_event_shapes(nat, ei, selected_jets);
	loop_timer.click("Event shapes calculation");
      }
    
      skf->compute_event_shapes(nat, ei, selected_jets);
      loop_timer.click("Event shapes calculation");
    }
    else if (skim_type == ksixb) {

      // Preselected jets are all jets in the event sorted in pT
      const DirectionalCut<int> cfg_nJets(config, "presel::njetsCut");
      // const bool applyPresel = config.readBoolOpt("presel::apply");
      // if (applyPresel) {
      if (!cfg_nJets.passedCut(presel_jets.size())) continue;
      cutflow.add("npresel_jets >= 6", nwt);
      cutflow_Unweighted.add("npresel_jets >= 6");
      // }
      
      //=============================================
      // Jets for pairing selection (either 6 or 0)
      //=============================================
      std::vector<Jet> selected_jets = skf->select_jets(nat, ei, presel_jets);

      if (!cfg_nJets.passedCut(selected_jets.size())) continue;
      cutflow.add("selected jets >= 6", nwt);
      cutflow_Unweighted.add("selected jets >= 6");

      loop_timer.click("Six b jet selection");
      
      ei.nfound_select = skf->n_gjmatched_in_jetcoll(nat, ei, selected_jets);
      ei.nfound_select_h = skf->n_ghmatched_in_jetcoll(nat, ei, selected_jets);

      if (debug) {
        dumpObjColl(selected_jets, "==== SELECTED 6b JETS ===");
      }
        
      
      //================================================
      // Proceed with the pairing of the 6 selected jets
      //=================================================
      skf->pair_jets(nat, ei, selected_jets);
      // cutflow.add("Jets for pairing selection");
      // cutflow_Unweighted.add("Jets for pairing selection");
      loop_timer.click("Six b jet pairing");
      ei.jet_list = selected_jets;

      
      if (is_signal) {
        skf->compute_seljets_genmatch_flags(nat, ei);
        loop_timer.click("Six b pairing flags");
      }
    
      skf->compute_event_shapes(nat, ei, selected_jets);
      loop_timer.click("Event shapes calculation");
    }
      
    else if (skim_type == kttbar) {
      if (presel_jets.size() < 2)
        continue;
      cutflow.add("npresel_jets>=2", nwt);
      cutflow_Unweighted.add("npresel_jets>=2");
      
      std::vector<Jet> ttjets = skf->select_jets(nat, ei, presel_jets); // ttjets sorted by DeepJet
      double deepjet1 = get_property(ttjets.at(0), Jet_btagDeepFlavB);
      double deepjet2 = get_property(ttjets.at(1), Jet_btagDeepFlavB);
      int nbtag = 0;
      if (deepjet1 > btag_WPs.at(bTagWP)) nbtag += 1;
      if (deepjet2 > btag_WPs.at(bTagWP)) nbtag += 1;
      if (nbtag < nMinBtag)
        continue;
      cutflow.add("ttbar_jet_cut", nwt);
      cutflow_Unweighted.add("ttbar_jet_cut");
      if (!is_data)
        ei.btagSF_WP_M = btsf.get_SF_allJetsPassWP({ttjets.at(0), ttjets.at(1)}, BtagSF::btagWP::medium);
      loop_timer.click("ttbar b jet selection");
    }
    else if (skim_type == kttbar) {
      if (presel_jets.size() < 2)
        continue;
      cutflow.add("npresel_jets>=2", nwt);
      
      std::vector<Jet> ttjets = skf->select_jets(nat, ei, presel_jets); // ttjets sorted by DeepJet
      double deepjet1 = get_property(ttjets.at(0), Jet_btagDeepFlavB);
      double deepjet2 = get_property(ttjets.at(1), Jet_btagDeepFlavB);
      int nbtag = 0;
      if (deepjet1 > btag_WPs.at(bTagWP)) nbtag += 1;
      if (deepjet2 > btag_WPs.at(bTagWP)) nbtag += 1;
      if (nbtag < nMinBtag)
        continue;
      cutflow.add("ttbar_jet_cut", nwt);
      if (!is_data)
        ei.btagSF_WP_M = btsf.get_SF_allJetsPassWP({ttjets.at(0), ttjets.at(1)}, BtagSF::btagWP::medium);
      loop_timer.click("ttbar b jet selection");
    }
    
    if (blind && is_data && skf->is_blinded(nat, ei, is_data))
      continue;
    
    if (!is_data && save_genw_tree){
      su::copy_gen_weights(ot, nwt);
      loop_timer.click("Read and copy gen weights");
    }
    
    su::fill_output_tree(ot, nat, ei);
    loop_timer.click("Output tree fill");
  } // Closes Event Loop 
  
  const auto end_loop_t = chrono::high_resolution_clock::now();

  outputFile.cd();
  cutflow.write(outputFile);
  cutflow_Unweighted.write(outputFile);
  ot.write();
  if (!is_data)
    nwt.write();
  const auto end_prog_t = chrono::high_resolution_clock::now();

  // timing statistics
  cout << endl;
  cout << "[INFO] : summary of skim loop execution time" << endl;
  loop_timer.print_summary();
  cout << endl;
  cout << "[INFO] : total elapsed time : " << chrono::duration_cast<chrono::milliseconds>(end_prog_t - start_prog_t).count()/1000.   << " s" << endl;
  cout << "       : startup time       : " << chrono::duration_cast<chrono::milliseconds>(start_loop_t - start_prog_t).count()/1000. << " s" << endl;
  cout << "       : loop time          : " << chrono::duration_cast<chrono::milliseconds>(end_loop_t - start_loop_t).count()/1000.   << " s" << endl;
  cout << "       : post-loop time     : " << chrono::duration_cast<chrono::milliseconds>(end_prog_t - end_loop_t).count()/1000.     << " s" << endl;

  cout << endl;
  cout << "[INFO] ... skim finished" << endl;
}
