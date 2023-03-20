#ifndef __majority_vote_circuit_h
#define __majority_vote_circuit_h

#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/gadgetlib1/gadgets/merkle_tree/merkle_tree_check_read_gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp>

#include "methods/MV.h"

using namespace libsnark;

template <typename FieldT>
class MajorityVoteCircuit: public gadget<FieldT> {

public: 
    MV &mv;
    std::vector<unsigned> &predicted_data;

public: 

    MajorityVoteCircuit(protoboard<FieldT> &pb, MV &_mv, std::vector<unsigned> &_predicted_data,const std::string &annotation = "") : gadget<FieldT>(pb, annotation), mv(_mv), predicted_data(_predicted_data) {
        
    }

    ~MajorityVoteCircuit() {}

    void generate_r1cs_constraints()
    {
        
    }

    void generate_r1cs_witness()
    {

    }
};

#endif