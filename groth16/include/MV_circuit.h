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
    std::vector<pb_variable<FieldT>> predicted_class_variables;
    std::vector<pb_variable<FieldT>> label_class_variables;
    std::vector<pb_variable<FieldT>> correct;

    // (auxiliary) secret witness
    std::vector<std::vector<pb_variable<FieldT>>> answer_variables;

    // (auxiliary) auxiliary witness and gadgets
    // label equality gadgets
    std::vector<std::vector<std::vector<EqualityCheckGadget<FieldT>>>> equality_gadgets; // task * worker * label
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> equality_results_variables;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> label_class_counts_variables;
    std::vector<std::vector<EqualityCheckGadget<FieldT>>> predicted_value_equality_gadgets;
    std::vector<std::vector<pb_variable<FieldT>>> predicted_value_equality_result_variables;
    std::vector<std::vector<pb_variable<FieldT>>> predicted_value_counts_variables;
    // comparison gadgets
    std::vector<std::vector<ComparisonGadget<FieldT>>> comparison_gadgets;
    std::vector<std::vector<pb_variable<FieldT>>> comparison_result_variables;;
    std::vector<std::vector<pb_variable<FieldT>>> comparison_diff_variables;
    // result equality gadgets
    std::vector<EqualityCheckGadget<FieldT>> result_checker_gadgets;

    // helper function
    // initialize all the variables pb needed
    void init_pb_vars()
    {
        auto &prefix = this->annotation_prefix;

        // zero variable
        zero_var.allocate(this->pb, prefix + std::string("zero_var"));

        // predicted variables;
        init_one_dimension_vec(predicted_class_variables, task_number, std::string("predicted"));
        // label class variables
        init_one_dimension_vec(label_class_variables, label_class_number, std::string("label_class"));
        // answer variables
        init_two_dimension_vec(answer_variables, task_number, worker_number, std::string("answer_variables"));
        // equality results variables
        init_three_dimension_vec(equality_results_variables, task_number, worker_number, label_class_number, std::string("equality_result"));
        // label class counts variables
        init_three_dimension_vec(label_class_counts_variables, task_number, worker_number + 1, label_class_number, std::string("label_class_count"));
        // comparison result variables
        init_two_dimension_vec(comparison_result_variables, task_number, label_class_number, std::string("comprison_result"));
        // comparison diff variables
        init_two_dimension_vec(comparison_diff_variables, task_number, label_class_number, std::string("comprison_diff"));
        init_two_dimension_vec(predicted_value_counts_variables, task_number, worker_number + 1, std::string("predicted_value_counts"));
        init_two_dimension_vec(predicted_value_equality_result_variables, task_number, worker_number, std::string("predicted_value_equality_result"));
        // correct
        init_one_dimension_vec(correct, task_number, std::string("correct"));

    }

    void assign_public_inputs()
    {
        // predicted variables
        for(int i = 0; i < task_number; i++) {
            this->pb.val(predicted_class_variables[i]) = predicted_data[i];
        }   

        // label class variables
        for (int i = 0; i < label_class_number; i++)
        {
            this->pb.val(label_class_variables[i]) = i;
        }

        // correct
        for (int i = 0; i < task_number; i++)
        {
            this->pb.val(correct[i]) = 1;
        }
    }

    void init_one_dimension_vec(std::vector<pb_variable<FieldT>> &pb_vec, int x_n, const std::string& prefix_name) {
        for (int i = 0; i < x_n; i++)
        {
            pb_variable<FieldT> var;
            var.allocate(this->pb, prefix_name + "_" + std::to_string(i));
            pb_vec.push_back(var);
        }
    }

    void init_two_dimension_vec(std::vector<std::vector<pb_variable<FieldT>>> &pb_vec, int x_n, int y_n, const std::string& prefix_name) {
        for (int i = 0; i < x_n; i++) {
            std::vector<pb_variable<FieldT>> row_x;
            for (int j = 0; j < y_n; j++)
            {
                pb_variable<FieldT> var;
                var.allocate(this->pb, prefix_name + "_" + std::to_string(i) + "_" + std::to_string(j));
                row_x.push_back(var);
            }
            pb_vec.push_back(row_x);
        }
    }

    void init_three_dimension_vec(std::vector<std::vector<std::vector<pb_variable<FieldT>>>> &pb_vec, int x_n, int y_n, int z_n, const std::string& prefix_name) {
        for (int i = 0; i < x_n; ++i)
        {
            std::vector<std::vector<pb_variable<FieldT>>> row_x;
            for (int j = 0; j < y_n; ++j)
            {
                std::vector<pb_variable<FieldT>> row_y;
                for (int k = 0; k < z_n; ++k)
                {
                    pb_variable<FieldT> var;
                    var.allocate(this->pb, prefix_name + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
                    row_y.push_back(var);
                }
                row_x.push_back(row_y);
            }
            pb_vec.push_back(row_x);
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

        // initialize comparison gadgets
        for(int i = 0; i < task_number; i++) {
            std::vector<ComparisonGadget<FieldT>> row_x;
            for(int j = 0; j < label_class_number; j++) {
                auto gadget_name = annotation + std::string("comparison_gadget_") + std::to_string(i) + "_" + std::to_string(j);
                row_x.push_back(ComparisonGadget<FieldT>(this->pb, label_class_counts_variables[i][worker_number][j], label_class_counts_variables[i][worker_number][predicted_data[i]], comparison_result_variables[i][j], comparison_diff_variables[i][j], gadget_name)); 
            }
            comparison_gadgets.push_back(row_x);
        }

        for (int i = 0; i < task_number; ++i)
        {
            std::vector<EqualityCheckGadget<FieldT>> row_x;
            for (int j = 0; j < worker_number; ++j)
            {
                auto gadget_name = annotation + std::string("predicted_value_equality_gadget_") + std::to_string(i) + "_" + std::to_string(j);
                row_x.push_back(EqualityCheckGadget<FieldT>(this->pb, answer_variables[i][j], predicted_class_variables[i], predicted_value_equality_result_variables[i][j], 4, gadget_name));
            }
            predicted_value_equality_gadgets.push_back(row_x);
        }

        // initialize result checker
        for (int i = 0; i < task_number; i++) {
            auto gadget_name = annotation + std::string("result_checker_") + std::to_string(i);
            result_checker_gadgets.push_back(EqualityCheckGadget<FieldT>(this->pb, predicted_value_counts_variables[i][worker_number], label_class_counts_variables[i][worker_number][predicted_data[i]], correct[i], 16, gadget_name)); 
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

        // predicted label equality
        for (int i = 0; i < task_number; i++) {
            for(int j = 0; j < worker_number; j++) {
                predicted_value_equality_gadgets[i][j].generate_r1cs_constraints();
                // count_i+1 - count_i = result
                add_r1cs(predicted_value_equality_result_variables[i][j], 1, predicted_value_counts_variables[i][j + 1] - predicted_value_counts_variables[i][j]);
            }
        }

        // compare and make sure the max count
        for(int i = 0; i < task_number; i++) {
            for(int j = 0; j < label_class_number; j++) {
                // compare with each label count, verify the maximum
                comparison_gadgets[i][j].generate_r1cs_constraints();
                // check result
                add_r1cs(comparison_result_variables[i][j], 1, 1);
            }
        }

        // correct checker
        for(int i = 0; i < task_number; i++) {
            result_checker_gadgets[i].generate_r1cs_constraints();
            add_r1cs(correct[i], 1, 1);
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
                this->pb.val(predicted_value_equality_result_variables[i][j]) = (mv.answer_data[i][j] == predicted_data[i]);
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
            int pcount = 0;
            for (int j = 0; j < worker_number + 1; j++)
            {
                // first item is 0
                if (j == 0)
                {
                    this->pb.val(predicted_value_counts_variables[i][0]) = 0;
                    for (int k = 0; k < label_class_number; k++)
                    {
                        this->pb.val(label_class_counts_variables[i][0][k]) = 0;
                    }
                }
                else
                {
                    pcount = (mv.answer_data[i][j-1] == predicted_data[i]) ? pcount + 1 : pcount;
                    this->pb.val(predicted_value_counts_variables[i][j]) = pcount;
                    for (int k = 0; k < label_class_number; k++)
                    {
                        counts[k] = (mv.answer_data[i][j - 1] == k) ? counts[k] + 1 : counts[k];
                        this->pb.val(label_class_counts_variables[i][j][k]) = counts[k];
                    }
                }
            }

            result_checker_gadgets[i].generate_r1cs_witness(counts[predicted_data[i]], pcount);
        }

        // comparison diff and result
        for (int i = 0; i < task_number; i++) {
            for(int j = 0; j < label_class_number; j++) {
                int x = mv.label_counts[i][j];
                int y = mv.label_counts[i][predicted_data[i]];
                this->pb.val(comparison_result_variables[i][j]) = (x <= y ? 1 : 0);
                this->pb.val(comparison_diff_variables[i][j]) = (x <= y ? y-x : x-y);
            }
        }

        // gadget
        for (int i = 0; i < task_number; i++)
        {
            for (int j = 0; j < worker_number; j++)
            {
                predicted_value_equality_gadgets[i][j].generate_r1cs_witness(mv.answer_data[i][j], predicted_data[i]);
                comparison_gadgets[i][j].generate_r1cs_witness();
                for (int k = 0; k < label_class_number; k++)
                {
                    equality_gadgets[i][j][k].generate_r1cs_witness(mv.answer_data[i][j], k);
                }
            }
        }
    }
};

#endif