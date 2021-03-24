#ifndef NORMWEIGHTTREE_H
#define NORMWEIGHTTREE_H

/**
** class  : NormWeightTree
** author : L. Cadamuro (UF)
** date   : 17/03/2021
** brief  : a tree of all the gen level weights for normalisation and syst tracking
**/

#include "BaseOutTree.h"
#include "NanoAODTree.h"
#include "PUWeightsReader.h"
#include <string>
#include <vector>
#include <map>

class NormWeightTree : public BaseOutTree {
    public:

        struct weight_t{
            double w;
            std::string name;
            std::vector<double>      syst_val;
            std::vector<std::string> syst_name;
        };

        NormWeightTree(std::string name = "NormWeightTree", std::string title = "NormWeightTree");
        ~NormWeightTree(){};

        // create the branches and initialise all the weight reading functions
        // pu_data is a map that must contain the following information to initialise PU reader (as strings)
        // - filename     : path/to/pu/histos
        // - name_PU_w    : nominal histogram name
        // - name_PU_w_up : up syst histogram name
        // - name_PU_w_do : down syst histogram name
        void init_weights(NanoAODTree& nat, std::map<std::string, std::string> pu_data);

        void read_weights(NanoAODTree& nat);


         // other custom weights added - they will be saved internally and will take part to the computation of the normalisation weight
         // FIXME: TODO
        // void add_weight(std::string wname, std::vector <std::string> wsystsname = {});

        // specific accessors from the nanoAOD tree for the weights
        void read_gen_weight   (NanoAODTree &nat);
        void read_pu_weight    (NanoAODTree &nat);
        void read_pdf_weight   (NanoAODTree &nat);
        void read_scale_weight (NanoAODTree &nat);
        void read_ps_weight    (NanoAODTree &nat);

        // specific initialisations
        void init_gen_weight   ();
        void init_pu_weight    (std::string filename, std::string name_PU_w, std::string name_PU_w_up, std::string name_PU_w_do); // need to read from the nat the number of variation
        void init_pdf_weight   (NanoAODTree &nat); // need to read from the nat the number of variation
        void init_scale_weight (NanoAODTree &nat); // need to read from the nat the number of variation
        void init_ps_weight    (NanoAODTree &nat); // need to read from the nat the number of variation

    private:

        weight_t gen_w_;
        weight_t pu_w_;
        weight_t pdf_w_;
        weight_t scale_w_;
        weight_t ps_w_;

        // special classes for specific weight access
        PUWeightsReader pu_wread_;

        void make_branches(weight_t& weight);

        // extra weights
        // bool extra_weights_branched_; // needed to control when branches are created for these
        std::vector<weight_t> extra_weights_;

        int verbosity_;
};
#endif