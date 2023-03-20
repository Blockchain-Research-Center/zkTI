#ifndef __zkp_groth16_h
#define __zkp_groth16_h

#include <iostream>

#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/examples/run_r1cs_gg_ppzksnark.hpp>

using namespace libsnark;

// mostly a wrapper of the Groth16 implementation from libsnark
template<typename ppT>
bool run_r1cs_gg_ppzksnark(const protoboard<libff::Fr<ppT>> &pb, const std::string& proof_out_filename="")
{
    std::cout << "Generate keys..." << std::endl;
    r1cs_gg_ppzksnark_keypair<ppT> keypair = r1cs_gg_ppzksnark_generator<ppT>(pb.get_constraint_system());
    r1cs_gg_ppzksnark_processed_verification_key<ppT> pvk = r1cs_gg_ppzksnark_verifier_process_vk<ppT>(keypair.vk);
    std::cout << "Done." << std::endl;

    std::cout << "Proving..." << std::endl;
    r1cs_gg_ppzksnark_proof<ppT> proof = r1cs_gg_ppzksnark_prover<ppT>(keypair.pk, pb.primary_input(), pb.auxiliary_input());
    std::cout << "Proof generated." << std::endl;

    std::cout << "Verifying..." << std::endl;
    const bool ans = r1cs_gg_ppzksnark_verifier_strong_IC<ppT>(keypair.vk, pb.primary_input(), proof);
    std::cout << "The verification result is: " << (ans ? "PASS" : "FAIL") << std::endl;

    if (proof_out_filename != "") {
        std::cout << "Saving proof to file: " << proof_out_filename << std::endl;
        std::ofstream out(proof_out_filename);
        out << proof;
        out.close();
    }

    return ans;
}

#endif 