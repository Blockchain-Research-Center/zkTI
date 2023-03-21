#ifndef __zkTI_majority_vote_circuit_h
#define __zkTI_majority_vote_circuit_h

#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/gadgetlib1/gadgets/merkle_tree/merkle_tree_check_read_gadget.hpp>
#include <libsnark/gadgetlib1/gadgets/hashes/sha256/sha256_gadget.hpp>

#include "methods/MV.h"
#include "gadgets.h"

using namespace libsnark;

#define add_r1cs(x, y, z) this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(x, y, z))

template <typename FieldT>
class MajorityVoteCircuit : public gadget<FieldT>
{

public:
    MV &mv;
    std::vector<unsigned> &predicted_data;
    int task_number, worker_number, label_number;
    int label_class_number;

private:
    pb_variable<FieldT> zero_var;

    // (primary) public input
    std::vector<pb_variable<FieldT>> predicted_variables;
    std::vector<pb_variable<FieldT>> label_class_variables;
    std::vector<std::vector<pb_variable<FieldT>>> label_class_total_counts_variable;

    // (auxiliary) secret witness
    std::vector<std::vector<pb_variable<FieldT>>> answer_variables;

    // (auxiliary) auxiliary witness and gadgets
    // label equality gadgets
    std::vector<std::vector<std::vector<EqualityCheckGadget<FieldT>>>> equality_gadgets; // task * worker * label
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> equality_results_variables;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> label_class_counts_variables;
    // // majority gadgets
    // std::vector<MajorityGadget<FieldT>> majority_gadgets;
    // std::vector<pb_variable<FieldT>> major_variables;
    // // result equality gadgets
    // std::vector<std::vector<EqualityCheckGadget<FieldT>>> result_checker_gadgets;


    // helper function
    // initialize all the variables pb needed
    void init_pb_vars()
    {
        auto &prefix = this->annotation_prefix;

        // zero variable
        zero_var.allocate(this->pb, prefix + std::string("zero_var"));

        // predicted variables;
        // TODO

        // label class variables
        for (int i = 0; i < label_class_number; i++)
        {
            pb_variable<FieldT> var;
            var.allocate(this->pb, std::string("label_class_") + std::to_string(i));
            label_class_variables.push_back(var);
        }

        // label class total counts variables
        for (int i = 0; i < task_number; i++)
        {
            std::vector<pb_variable<FieldT>> row_x;
            for (int j = 0; j < label_class_number; j++)
            {
                pb_variable<FieldT> var;
                var.allocate(this->pb, std::string("label_class_total_counts_") + std::to_string(i) + "_" + std::to_string(j));
                row_x.push_back(var);
            }
            label_class_total_counts_variable.push_back(row_x);
        }

        // answer variables
        for (int i = 0; i < task_number; i++)
        {
            std::vector<pb_variable<FieldT>> row_x;
            for (int j = 0; j < worker_number; j++)
            {
                pb_variable<FieldT> var;
                var.allocate(this->pb, std::string("answer_variables_") + std::to_string(i) + "_" + std::to_string(j));
                row_x.push_back(var);
            }
            answer_variables.push_back(row_x);
        }

        // equality results variables
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<std::vector<pb_variable<FieldT>>> row_x;
            for (int j = 0; j < worker_number; ++j)
            {
                std::vector<pb_variable<FieldT>> row_y;
                for (int k = 0; k < label_class_number; ++k)
                {
                    pb_variable<FieldT> var;
                    var.allocate(this->pb, std::string("equality_result_") + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
                    row_y.push_back(var);
                }
                row_x.push_back(row_y);
            }
            equality_results_variables.push_back(row_x);
        }

        // label class counts variables
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<std::vector<pb_variable<FieldT>>> row_x;
            for (int j = 0; j < worker_number + 1; ++j)
            {
                ;
                std::vector<pb_variable<FieldT>> row_y;
                for (int k = 0; k < label_class_number; ++k)
                {
                    pb_variable<FieldT> var;
                    var.allocate(this->pb, std::string("label_class_count_") + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
                    row_y.push_back(var);
                }
                row_x.push_back(row_y);
            }
            label_class_counts_variables.push_back(row_x);
        }
    }

    void assign_public_inputs()
    {
        // answer variables
        // for (int i = 0; i < task_number; i++)
        // {
        //     for (int j = 0; j < worker_number; j++) {
        //         this->pb.val(answer_variables[i][j]) = this->mv.answer_data[i][j];
        //     }
        // }

        // predicted variables
        // TODO

        // label class variables
        for (int i = 0; i < label_class_number; i++)
        {
            this->pb.val(label_class_variables[i]) = i;
        }

        // label class total count variables
        for (int i = 0; i < task_number; i++)
        {
            for (int j = 0; j < label_class_number; j++)
            {
                this->pb.val(label_class_total_counts_variable[i][j]) = mv.label_counts[i][j];
            }
        }
    }

public:
    MajorityVoteCircuit(protoboard<FieldT> &pb, MV &_mv, std::vector<unsigned> &_predicted_data, const std::string &annotation = "") : gadget<FieldT>(pb, annotation), mv(_mv), predicted_data(_predicted_data)
    {
        task_number = mv.answer_data.size();
        worker_number = mv.answer_data[0].size();
        label_number = task_number * worker_number;
        label_class_number = mv.label_class_number;

        // initialize variables
        init_pb_vars();

        // assign public inputs
        assign_public_inputs();

        // initialize equality gadgets
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<std::vector<EqualityCheckGadget<FieldT>>> row_x;
            for (int j = 0; j < worker_number; ++j)
            {
                ;
                std::vector<EqualityCheckGadget<FieldT>> row_y;
                for (int k = 0; k < label_class_number; ++k)
                {
                    auto gadget_name = annotation + std::string("equality_gadget_") + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                    // max label class number is 16
                    row_y.push_back(EqualityCheckGadget<FieldT>(this->pb, answer_variables[i][j], label_class_variables[k], equality_results_variables[i][j][k], 4, gadget_name));
                }
                row_x.push_back(row_y);
            }
            equality_gadgets.push_back(row_x);
        }
    }

    ~MajorityVoteCircuit() {}

    void generate_r1cs_constraints()
    {
        // label equality
        for (int i = 0; i < task_number; i++)
        {
            for (int j = 0; j < worker_number; j++)
            {
                for (int k = 0; k < label_class_number; k++)
                {
                    equality_gadgets[i][j][k].generate_r1cs_constraints();
                    // count_i+1 - count_i = result
                    add_r1cs(equality_results_variables[i][j][k], 1, label_class_counts_variables[i][j + 1][k] - label_class_counts_variables[i][j][k]);
                }
            }
        }
    }

    void generate_r1cs_witness()
    {
        // answer variables
        for (int i = 0; i < task_number; i++)
        {
            for (int j = 0; j < worker_number; j++)
            {
                this->pb.val(answer_variables[i][j]) = mv.answer_data[i][j];
            }
        }

        // equality results variables
        for (int i = 0; i < task_number; i++)
        {
            for (int j = 0; j < worker_number; j++)
            {
                for (int k = 0; k < label_class_number; k++)
                {
                    this->pb.val(equality_results_variables[i][j][k]) = (mv.answer_data[i][j] == k);
                }
            }
        }

        // label class counts variables
        for (int i = 0; i < task_number; i++)
        {
            int counts[label_class_number] = {0};
            for (int j = 0; j < worker_number + 1; j++)
            {
                // first item is 0
                if (j == 0)
                {
                    for (int k = 0; k < label_class_number; k++)
                    {
                        this->pb.val(label_class_counts_variables[i][0][k]) = 0;
                    }
                }
                else
                {
                    for (int k = 0; k < label_class_number; k++)
                    {
                        counts[k] = (mv.answer_data[i][j - 1] == k) ? counts[k] + 1 : counts[k];
                        this->pb.val(label_class_counts_variables[i][j][k]) = counts[k];
                    }
                }
            }
        }

        for (int i = 0; i < task_number; i++)
        {
            for (int j = 0; j < worker_number; j++)
            {
                for (int k = 0; k < label_class_number; k++)
                {
                    equality_gadgets[i][j][k].generate_r1cs_witness(mv.answer_data[i][j], k);
                }
            }
        }
    }
};

#endif