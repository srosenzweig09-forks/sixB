#include "EvalNN.h"

#include <iostream>
#include <string>
#include <vector>

using namespace std;

void test_6j_multi_classifier()
{
  string f_classifier = "/uscms/home/srosenzw/nobackup/workarea/higgs/sixb_analysis/CMSSW_10_2_18/src/sixb/6jet_classifier/models/20210828_multiclass_nsignal/";
  cout << "[INFO] Testing Classifier: " << f_classifier << endl;

  EvalNN network(f_classifier);
  // network.set_output("training/Nadam/dense_2/bias/v");
  network.set_debug(true);
    
  vector<float> inputs = {120.153,  114.914,  75.994,  75.099,  61.896,  61.564, -0.525,
                          -0.278, -0.992,  0.311, -1.057, -1.844, -1.140,  1.660,  0.878,
                          2.717,  3.042, -1.256,  0.778,  0.700,  0.401,  0.996,  0.003,
                          0.996,  108.143,  124.577,  80.956,  79.714,  67.812,  42.099};
  vector<float> outputs = network.evaluate(inputs);
}

void test_6j_classifier()
{
  string f_classifier = "models/6jet_classifier/";

  cout << "[INFO] Testing Classifier: " << f_classifier << endl;

  // tensorflow::GraphDef* graphDef = tensorflow::loadGraphDef(f_2j_classifier);
  EvalNN network(f_classifier);
  network.set_debug(true);
    
  vector<float> inputs = {120.153,  114.914,  75.994,  75.099,  61.896,  61.564, -0.525,
                          -0.278, -0.992,  0.311, -1.057, -1.844, -1.140,  1.660,  0.878,
                          2.717,  3.042, -1.256,  0.778,  0.700,  0.401,  0.996,  0.003,
                          0.996,  108.143,  124.577,  80.956,  79.714,  67.812,  42.099};
    
  vector<float> outputs = network.evaluate(inputs);
}

void test_2j_classifier()
{
  
  string f_classifier = "models/2jet_classifier/";

  cout << "[INFO] Testing Classifier: " << f_classifier << endl;

  // tensorflow::GraphDef* graphDef = tensorflow::loadGraphDef(f_classifier);
  EvalNN network(f_classifier);
  network.set_debug(true);
    
  vector<float> inputs = {167.819839, 94.005722,
			  1.685547, 1.035645,
			  -1.925537, -0.478943,
			  0.922363, 0.997559,
			  1.585878};
    
  vector<float> outputs = network.evaluate(inputs);
}

int main() {
  test_6j_multi_classifier();
  test_6j_classifier();
  // test_2j_classifier();
}
