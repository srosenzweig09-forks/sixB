import modules.Sample as sam
import collections
import ROOT
import copy
import numpy as np
import hep_ml

ROOT.gInterpreter.ProcessLine('#include <cmath>')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                      SAMPLES
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
year = "2018"

from modules.Samples_NMSSM_XYH_YToHH_6b import dsetGroups
samples = dsetGroups["NMSSM_XYH_YToHH_6b"][year] + dsetGroups["QCD"][year] + dsetGroups["TT"][year]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#                  ROOT FUNCTIONS
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
average_btag = """
double average_btag(std::vector<double> bdiscs) {
  double sum = 0.0;
  for (unsigned int i=0; i<bdiscs.size(); i++)
    {
      sum += bdiscs.at(i);
    }
  return sum/bdiscs.size();
}
"""

HTb = """
double HTb(std::vector<double> pts) {
  double sum = 0.0;
  for (unsigned int i=0; i<pts.size(); i++)
   { 
    sum += pts.at(i);
  }
  return sum;
}
"""

n_medium_btag = """
int n_medium_btag(double btag1, double btag2, double btag3, double btag4, double btag5, double btag6) {
    double btagwp = 0.2783;
    std::vector<double> btags = {btag1,btag2,btag3,btag4,btag5,btag6};
    int n_medium = 0;
    for (double btag : btags) {
        if (btag > btagwp) {
            n_medium++;
        }
    }
    return n_medium;
}
"""

Dm_cand = """
double Dm_cand(double HX_m, double HY1_m, double HY2_m, double center)
{
  double dm_HX  = std::abs(HX_m - center);
  double dm_HY1 = std::abs(HY1_m - center);
  double dm_HY2 = std::abs(HY2_m - center);

  dm_HX = dm_HX * dm_HX;
  dm_HY1 = dm_HY1 * dm_HY1;
  dm_HY2 = dm_HY2 * dm_HY2;

  double sum = dm_HX + dm_HY1 + dm_HY2;
  double dm_cand = sqrt(sum);
  return dm_cand;
}
"""

Dm_cand_ver = """
double Dm_cand_ver(double HX_m, double HY1_m, double HY2_m, double center)
{
  double dm_HX  = std::abs(HX_m - center);
  double dm_HY1 = std::abs(HY1_m - center);
  double dm_HY2 = std::abs(HY2_m - center);

  dm_HX = dm_HX * dm_HX;
  dm_HY1 = dm_HY1 * dm_HY1;
  dm_HY2 = dm_HY2 * dm_HY2;

  double sum = dm_HX + dm_HY1 + dm_HY2;
  double dm_cand = sqrt(sum);
  return dm_cand;
}
"""

## all variables to declare to the gInterpreter should be listed in a "declarations" variable
declarations = [n_medium_btag, HTb, average_btag, Dm_cand, Dm_cand_ver]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#              New columns to declare
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
new_columns                  = collections.OrderedDict()
new_columns["leptonveto"]    = "n_muon == 0 && n_ele == 0"
new_columns["HTb"]           = 'HTb({HY1_b1_pt, HY1_b2_pt, HY2_b1_pt, HY2_b2_pt, HX_b1_pt, HX_b2_pt})'
new_columns['n_medium_btag'] = 'n_medium_btag(HY1_b1_btag, HY1_b2_btag, HY2_b1_btag, HY2_b2_btag, HX_b1_btag, HX_b2_btag)'
new_columns['average_btag']  = 'average_btag({HY1_b1_btag, HY1_b2_btag, HY2_b1_btag, HY2_b2_btag, HX_b1_btag, HX_b2_btag})'
new_columns['Dm_cand']       = 'Dm_cand(HX_m, HY1_m, HY2_m, 125.0)'
new_columns['Dm_cand_ver']   = 'Dm_cand_ver(HX_m, HY1_m, HY2_m, 174.0)'
new_columns['A_SR_mask']     = 'Dm_cand <= 30.0'
new_columns['A_CR_mask']     = '(Dm_cand > 30.0) && (Dm_cand <= 40.0)' # Analysis CR
new_columns['ls_mask']       = 'average_btag < 0.6'
new_columns['hs_mask']       = 'average_btag >= 0.6'
new_columns['A_CRls_mask']   = 'A_CR_mask && ls_mask'
new_columns['A_CRhs_mask']   = 'A_CR_mask && hs_mask'
new_columns['A_SRls_mask']   = 'A_SR_mask && ls_mask'
new_columns['A_SRhs_mask']   = 'A_SR_mask && hs_mask'
#new_columns['V_SR_mask']     = 'Dm_cand_ver <= 30.0'
#new_columns['V_CR_mask']     = '(Dm_cand_ver > 30.0) && (Dm_cand_ver <= 40.0)' # Validation CR
#new_columns['V_CRls_mask']   = 'V_CR_mask && ls_mask'
#new_columns['V_CRhs_mask']   = 'V_CR_mask && hs_mask'
#new_columns['V_SRls_mask']   = 'V_SR_mask && ls_mask'
#new_columns['V_SRhs_mask']   = 'V_SR_mask && hs_mask'
  
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#              Perform plots for all selections
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
selections_defs = {
    'baseline'            : 'n_jet > 0',
    'lepton_veto'         : 'leptonveto == 1',
    'A_SR_highbtag'       : 'A_SRhs_mask == 1',
    'A_CR_highbtag'       : 'A_CRhs_mask == 1',
    'A_SR_lowbtag'        : 'A_SRls_mask == 1',
    'A_CR_lowbtag'        : 'A_CRls_mask == 1',
#    'V_SR_highbtag'       : 'V_SRhs_mask == 1',
#    'V_CR_highbtag'       : 'V_CRhs_mask == 1',
#    'V_SR_lowbtag'        : 'V_SRls_mask == 1',
#    'V_CR_lowbtag'        : 'V_CRls_mask == 1',
}

###################### WEIGHTS DEFINITION #######################
# weights that appear in norm_weights undergo a rescaling by sum(w) when filling the plots
# these weights must be defined both in the norm_tree and in the event_tree
# norm_weights = ['genWeight', 'PUWeight']
norm_weights = ['genWeight']

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Histograms to be created
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minMX = 375.0
maxMX = 2000.0
nbins = int(65)

x = 375
binList = [x]
for i in range(0, 65):
    x+= 25.0
    binList.append(x)
mBins = tuple(binList)

histos_descs = [ {'var': 'HTb', 'weightlist': ['genWeight'], 'bins': (100, 0.0, 2000),},
                 {'var': 'X_m', 'weightlist': ['genWeight'], 'nbins': nbins, 'min': minMX, 'max': maxMX},
#                 {'var': 'n_jet', 'weightlist': ['genWeight'], 'bins' : (10, 6, 16),},
#                 {'var': 'HY1_b1_pt', 'weightlist': ['genWeight'], 'bins' : (80, 0.0, 800),},
#                 {'var': 'HY2_b1_pt', 'weightlist': ['genWeight'], 'bins' : (80, 0.0, 800),},
#                 {'var': 'HX_b1_pt', 'weightlist': ['genWeight'], 'bins' : (80, 0.0, 800),},
#                 {'var': 'HY1_b2_pt', 'weightlist': ['genWeight'], 'bins' : (80, 0.0, 800),},
#                 {'var': 'HY2_b2_pt', 'weightlist': ['genWeight'], 'bins' : (80, 0.0, 800),},
#                 {'var': 'HX_b2_pt', 'weightlist': ['genWeight'], 'bins' : (80, 0.0, 800),},
#                 {'var': 'HX_m', 'weightlist': ['genWeight'], 'bins': (300, 0.0, 600),},
#                 {'var': 'HY1_m', 'weightlist': ['genWeight'], 'bins': (300, 0.0, 600),},
#                 {'var': 'HY2_m', 'weightlist': ['genWeight'], 'bins': (300, 0.0, 600),},
#                 {'var': 'HX_b1_btag', 'weightlist': ['genWeight'], 'bins': (50, 0.0, 1.0),},
#                 {'var': 'HX_b2_btag', 'weightlist': ['genWeight'], 'bins': (50, 0.0, 1.0),},
#                 {'var': 'HY1_b1_btag', 'weightlist': ['genWeight'], 'bins': (50, 0.0, 1.0),},
#                 {'var': 'HY1_b2_btag', 'weightlist': ['genWeight'], 'bins': (50, 0.0, 1.0),},
#                 {'var': 'HY2_b1_btag', 'weightlist': ['genWeight'], 'bins': (50, 0.0, 1.0),},
#                 {'var': 'HY2_b2_btag','weightlist': ['genWeight'], 'bins': (50, 0.0, 1.0),},
#                 {'var': "average_btag", 'weightlist': ['genWeight'], 'bins': (100, 0.0, 1.0),},
]

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Create a histogram for each selection declared
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
histos = []
for hd in histos_descs:
    for s in selections_defs.keys():
        d2 = copy.deepcopy(hd)
        d2['sel'] = s
        histos.append(d2)
