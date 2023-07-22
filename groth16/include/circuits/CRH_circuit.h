#ifndef __zkTI_crh_circuit_h
#define __zkTI_crh_circuit_h

#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include "../methods/Float.h"
#include "../methods/CRH.h"
#include "./gadgets.h"
#include "./gadgets2.h"
#include "../common.h"

using namespace libsnark;

#ifndef __BASE
#define __BASE 2
#endif

template <typename FieldT>
class CRHCircuit : public gadget<FieldT>
{

public:
    CRH &crh;
    std::vector<unsigned> &predicted_data;
    int task_number, worker_number, label_number;
    int label_class_number;
    float init_quality;

private:
    // (primary input) public input
    pb_variable<FieldT> base;
    pb_variable<FieldT> random_point, random_coefficients;
    std::vector<pb_variable<FieldT>> task_n_inv; // 1 / task_number
    std::vector<pb_variable<FieldT>> multiplies; // base**base**i
    std::vector<pb_variable<FieldT>> label_class_variables;

    // (auxiliary input) secret witness
    std::vector<std::vector<pb_variable<FieldT>>> answer_variables;    std::vector<std::vector<pb_variable<FieldT>>> worker_quality;
    std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>> task_label_prediction;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> updated_task_label_prediction;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> updated_worker_quality;
    std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>> add_result_variables;
    std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>> worker_chosen_label_prediction;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> equality_result_variables;

    // (auxiliary input) auxiliary witness and gadgets
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> sum;
    // label equality gadgets
    std::vector<std::vector<std::vector<EqualityCheckGadget<FieldT>*>>> equality_check_gadgets; // task * worker * label
    // add gates 1
    std::vector<std::vector<std::vector<FloatAddGadget<FieldT>*>>> add_gates_1;
    // add gates 2
    std::vector<std::vector<FloatAddGadget<FieldT>*>> add_gates_2;
    // divide gates
    std::vector<std::vector<FloatDivideGadget<FieldT>*>> divide_gates;

    // helper function
    // initialize all the variables pb needed
    void init_pb_vars()
    {
        auto &prefix = this->annotation_prefix;

        random_point.allocate(this->pb, std::string("random_point"));
        random_coefficients.allocate(this->pb, std::string("random_coefficients"));
        base.allocate(this->pb, std::string("base"));
        init_one_dimension_vec(this->pb, task_n_inv, 2, std::string("task_n_inv"));
        init_one_dimension_vec(this->pb, multiplies, 8, std::string("base_multiplies"));
        // label class variables
        init_one_dimension_vec(this->pb, label_class_variables, label_class_number, std::string("label_class"));

        init_two_dimension_vec(this->pb, worker_quality, worker_number, 2 * 2, "worker_quality");
        init_three_dimension_vec(this->pb, updated_worker_quality, worker_number, task_number + 1, 2, "updated_worker_quality");
        init_four_dimension_vec(this->pb, task_label_prediction, task_number, label_class_number, worker_number + 1, 2, "task_label_prediction");
        // task prediction result
        init_three_dimension_vec(this->pb, updated_task_label_prediction, task_number, label_class_number, 2, std::string("updated_task_label_prediction"));
        // answer variables
        init_two_dimension_vec(this->pb, answer_variables, task_number, worker_number, std::string("answer_variables"));
        // sum
        init_three_dimension_vec(this->pb, sum, task_number, label_class_number + 1, 2, std::string("sum"));
        // equality check result
        init_three_dimension_vec(this->pb, equality_result_variables, task_number, worker_number, label_class_number, std::string("equality_result"));
        // add result
        init_four_dimension_vec(this->pb, add_result_variables, task_number, label_class_number, worker_number + 1, 2, "add_result_variables");
        // worker choosed label prediction
        init_four_dimension_vec(this->pb, worker_chosen_label_prediction, worker_number, task_number, label_class_number + 1, 2, "worker_chosen_label_prediction");
        // multiplication result 2
        init_three_dimension_vec(this->pb, multiply_result_variables_2, worker_number, task_number, 2, std::string("multiply_result_variables_2"));

        // set public number
        this->pb.set_input_sizes(8);
    }

    // assign values for public inputs
    void assign_public_inputs()
    {
        this->pb.val(base) = __BASE;
        this->pb.val(random_point) = FieldT();
        this->pb.val(random_coefficients) = FieldT();
        this->pb.val(task_n_inv[0]) = Float(1.0f/task_number).value;
        this->pb.val(task_n_inv[1]) = Float(1.0f/task_number).scale; 

        for(int i = 0; i < 8; i++) {
            mpz_t base, tmp;
            mpz_init(tmp);
            mpz_init_set_ui(base, __BASE);
            mpz_pow_ui(tmp, base, (u64)pow(__BASE, i));
            this->pb.val(multiplies[i]) = libff::bigint<__limbs>(tmp);
        }

        // label class variables
        for (int i = 0; i < label_class_number; i++)
        {
            this->pb.val(label_class_variables[i]) = i;
        }

        // worker_quality_value
        for(int j = 0; j < worker_number; j++) {
            // 0.8
            this->pb.val(worker_quality[j][0]) = Float(init_quality).value;
            this->pb.val(worker_quality[j][1]) = Float(init_quality).scale;   

            // 1 - 0.8
            this->pb.val(worker_quality[j][2]) = Float(init_quality).one_minus().value;
            this->pb.val(worker_quality[j][3]) = Float(init_quality).one_minus().scale;         
        }

        // task_label_prediction
        for(int i = 0; i < task_number; i++) {
            for(int k = 0; k < label_class_number; k++) {
                // 1
                this->pb.val(task_label_prediction[i][k][0][0]) = Float::one().value;
                this->pb.val(task_label_prediction[i][k][0][1]) = Float::one().scale;
            }
        }

        // label class variables
        for (int i = 0; i < label_class_number; i++)
        {
            this->pb.val(label_class_variables[i]) = i;
        }
    }

public:
    CRHCircuit(protoboard<FieldT> &pb, CRH &_crh, std::vector<unsigned> &_predicted_data, const std::string &annotation = "") : gadget<FieldT>(pb, annotation), crh(_crh), predicted_data(_predicted_data)
    {
        task_number = crh.answer_data.size();
        // task_number = 10;
        worker_number = crh.answer_data[0].size();
        // worker_number = 2;
        label_number = task_number * worker_number;
        label_class_number = crh.label_class_number;
        init_quality = crh.init_quality.real_value;

        // initialize variables
        init_pb_vars();

        // assign public inputs
        assign_public_inputs();;

        // initialize equality gadgets
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<std::vector<EqualityCheckGadget<FieldT>*>> row_x;
            for (int j = 0; j < worker_number; ++j)
            {
                std::vector<EqualityCheckGadget<FieldT>*> row_y;
                for (int k = 0; k < label_class_number; ++k)
                {
                    auto gadget_name = annotation + std::string("equality_check_gadget_") + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);
                    // max label class number is 4 (2 bits)
                    row_y.push_back(new EqualityCheckGadget<FieldT>(this->pb, answer_variables[i][j], label_class_variables[k], equality_result_variables[i][j][k], 2, gadget_name));
                }
                row_x.push_back(row_y);
            }
            equality_check_gadgets.push_back(row_x);
        }

        // initialize multiply gadgets 1
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<std::vector<std::vector<FloatMultiplyGadget<FieldT>*>>> row_x;
            for (int j = 0; j < worker_number; ++j)
            {
                std::vector<std::vector<FloatMultiplyGadget<FieldT>*>> row_y;
                for (int k = 0; k < label_class_number; k++) {
                    std::vector<FloatMultiplyGadget<FieldT>*> row_k;
                    auto gadget_name = annotation + std::string("multiply_gates_1_") + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);   

                    row_k.push_back(new FloatMultiplyGadget<FieldT>(this->pb, task_label_prediction[i][k][j][0], task_label_prediction[i][k][j][1], worker_quality[j][0], worker_quality[j][1], multiply_result_variables[i][k][j][0][0], multiply_result_variables[i][k][j][0][1], gadget_name + "_0"));

                    row_k.push_back(new FloatMultiplyGadget<FieldT>(this->pb, task_label_prediction[i][k][j][0], task_label_prediction[i][k][j][1], worker_quality[j][2], worker_quality[j][3], multiply_result_variables[i][k][j][1][0], multiply_result_variables[i][k][j][1][1], gadget_name + "_1"));      

                    row_y.push_back(row_k);
                }
                row_x.push_back(row_y);
            }
            multiply_gates_1.push_back(row_x);
        }

        // initialize add gadgets 2
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<FloatAddGadget<FieldT>*> row_x;
            for (int k = 0; k < label_class_number; ++k)
            {
                auto gadget_name = annotation + std::string("add_gates_1_") + std::to_string(i) + "_" + std::to_string(k);   
                // we leave the implementation of permutation check inside gadget
                row_x.push_back(new FloatAddGadget<FieldT>(this->pb, sum[i][k][0], sum[i][k][1], task_label_prediction[i][k][worker_number][0], task_label_prediction[i][k][worker_number][1], sum[i][k+1][0], sum[i][k+1][1], multiplies, random_point, random_coefficients, gadget_name));      
            }
            add_gates_2.push_back(row_x);
        }

        // initialize divide gadgets
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<FloatDivideGadget<FieldT>*> row_x;
            for (int k = 0; k < label_class_number; ++k)
            {
                auto gadget_name = annotation + std::string("divide_gates_") + std::to_string(i) + "_" + std::to_string(k);   
                row_x.push_back(new FloatDivideGadget<FieldT>(this->pb, task_label_prediction[i][k][worker_number][0], task_label_prediction[i][k][worker_number][1], sum[i][label_class_number][0], sum[i][label_class_number][1], updated_task_label_prediction[i][k][0], updated_task_label_prediction[i][k][1], gadget_name));      
            }
            divide_gates.push_back(row_x);
        }

        // initialize multiply gadgets 2
        for (int j = 0; j < worker_number; ++j)
        {
            std::vector<FloatMultiplyGadget<FieldT>*> row_x;
            for (int i = 0; i < task_number; ++i)
            {
                auto gadget_name = annotation + std::string("multiply_gates_2_") + std::to_string(j) + "_" + std::to_string(i);   
                row_x.push_back(new FloatMultiplyGadget<FieldT>(this->pb, worker_chosen_label_prediction[j][i][label_class_number][0], worker_chosen_label_prediction[j][i][label_class_number][1], task_n_inv[0], task_n_inv[1], multiply_result_variables_2[j][i][0], multiply_result_variables_2[j][i][1], gadget_name));      
            }
            multiply_gates_2.push_back(row_x);
        }

        // initialize add gadgets 2
        for (int j = 0; j < worker_number; ++j)
        {
            std::vector<FloatAddGadget<FieldT>*> row_x;
            for (int i = 0; i < task_number; ++i)
            {
                auto gadget_name = annotation + std::string("add_gates_2_") + std::to_string(j) + "_" + std::to_string(i);   
                // we leave the implementation of permutation check inside gadget
                row_x.push_back(new FloatAddGadget<FieldT>(this->pb, updated_worker_quality[j][i][0], updated_worker_quality[j][i][1], multiply_result_variables_2[j][i][0], multiply_result_variables_2[j][i][1], updated_worker_quality[j][i+1][0], updated_worker_quality[j][i+1][1], multiplies, random_point, random_coefficients, gadget_name));      
            }
            add_gates_2.push_back(row_x);
        }
    }

    ~CRHCircuit() {}

    void generate_r1cs_constraints()
    {
        for(int i = 0; i < task_number; i++) {
            for(int j = 0; j < worker_number; j++) {
                for(int k = 0; k < label_class_number; k++) {
                    multiply_gates_1[i][j][k][0]->generate_r1cs_constraints();
                    multiply_gates_1[i][j][k][1]->generate_r1cs_constraints();
                    equality_check_gadgets[i][j][k]->generate_r1cs_constraints();
                    // b*(s_i+1-a) = (1-b)*(s_i+1-b)
                    add_r1cs(equality_result_variables[i][j][k], 2 * task_label_prediction[i][k][j+1][0] - multiply_result_variables[i][k][j][0][0] - multiply_result_variables[i][k][j][1][0], task_label_prediction[i][k][j+1][0] - multiply_result_variables[i][k][j][1][0]);
                    add_r1cs(equality_result_variables[i][j][k], 2 * task_label_prediction[i][k][j+1][1] - multiply_result_variables[i][k][j][0][1] - multiply_result_variables[i][k][j][1][1], task_label_prediction[i][k][j+1][1] - multiply_result_variables[i][k][j][1][1]);
                }
            }

            for(int k = 0; k < label_class_number; k++) {
                add_gates_2[i][k]->generate_r1cs_constraints();
                divide_gates[i][k]->generate_r1cs_constraints();
            }
        }

        for(int j = 0; j < worker_number; j++) {
            for(int i = 0; i < task_number; i++) {
                // for(int k = 0; k < label_class_number; k++) {
                //     // label = answer_data[i][j]
                //     // task_label_prediction[i][label]
                //     add_r1cs(updated_task_label_prediction[i][k][0], equality_result_variables[i][j][k], worker_chosen_label_prediction[j][i][k+1][0] - worker_chosen_label_prediction[j][i][k][0]);
                //     add_r1cs(updated_task_label_prediction[i][k][1], equality_result_variables[i][j][k], worker_chosen_label_prediction[j][i][k+1][1] - worker_chosen_label_prediction[j][i][k][1]);
                // }
                // multiply_gates_2[j][i]->generate_r1cs_constraints();
                // add_gates_2[j][i]->generate_r1cs_constraints();
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
                this->pb.val(answer_variables[i][j]) = crh.answer_data[i][j];
            }
        }

        // CRH
        std::vector<Float> worker_quality_value = std::vector<Float>(worker_number);
        for(int j = 0; j < worker_number; j++) {
            worker_quality_value[j] = Float(init_quality);
        }
        std::vector<std::vector<Float>> task_label_prediction_value = std::vector<std::vector<Float>>(task_number);
        for(int i = 0; i < task_number; i++) {
            task_label_prediction_value[i].resize(label_class_number);
            for(int k = 0; k < label_class_number; k++) {
                task_label_prediction_value[i][k] = Float::one();
            }

            for(int j = 0; j < worker_number; j++) {
                unsigned label = crh.answer_data[i][j];
                for(int k = 0; k < label_class_number; k++) {
                    Float a = worker_quality_value[j];
                    Float b = a.one_minus();
                    Float multi_value_0 = task_label_prediction_value[i][k] * a;
                    Float multi_value_1 = task_label_prediction_value[i][k] * b;
                    
                    multiply_gates_1[i][j][k][0]->generate_r1cs_witness(task_label_prediction_value[i][k], a, multi_value_0);
                    multiply_gates_1[i][j][k][1]->generate_r1cs_witness(task_label_prediction_value[i][k], b, multi_value_1);   

                    this->pb.val(multiply_result_variables[i][k][j][0][0]) = multi_value_0.value;
                    this->pb.val(multiply_result_variables[i][k][j][0][1]) = multi_value_0.scale;
                    this->pb.val(multiply_result_variables[i][k][j][1][0]) = multi_value_1.value;
                    this->pb.val(multiply_result_variables[i][k][j][1][1]) = multi_value_1.scale;

                    if(label == k) {
                        task_label_prediction_value[i][k] = multi_value_0;
                    } else {
                        task_label_prediction_value[i][k] = multi_value_1;
                    }

                    this->pb.val(equality_result_variables[i][j][k]) = (crh.answer_data[i][j] == k);
                    equality_check_gadgets[i][j][k]->generate_r1cs_witness(crh.answer_data[i][j], k);

                    this->pb.val(task_label_prediction[i][k][j+1][0]) = task_label_prediction_value[i][k].value;
                    this->pb.val(task_label_prediction[i][k][j+1][1]) = task_label_prediction_value[i][k].scale;
                }
            }

            Float sums_value = Float::zero();
            this->pb.val(sum[i][0][0]) = sums_value.value;
            this->pb.val(sum[i][0][1]) = sums_value.scale;
            for(int k = 0; k < label_class_number; k++) {
                Float added_value = sums_value + task_label_prediction_value[i][k];
                this->pb.val(sum[i][k+1][0]) = added_value.value;
                this->pb.val(sum[i][k+1][1]) = added_value.scale;
                add_gates_2[i][k]->generate_r1cs_witness(sums_value, task_label_prediction_value[i][k], added_value);                      
                sums_value = added_value;
            }

            for(int k = 0; k < label_class_number; k++) {
                Float divided_value = task_label_prediction_value[i][k] / sums_value;
                this->pb.val(updated_task_label_prediction[i][k][0]) = divided_value.value;
                this->pb.val(updated_task_label_prediction[i][k][1]) = divided_value.scale;
                divide_gates[i][k]->generate_r1cs_witness(task_label_prediction_value[i][k], sums_value, divided_value);
                task_label_prediction_value[i][k] = divided_value;
            }
        }

        for(int j = 0; j < worker_number; j++) {
            worker_quality_value[j] = Float::zero();
            this->pb.val(updated_worker_quality[j][0][0]) = worker_quality_value[j].value;
            this->pb.val(updated_worker_quality[j][0][1]) = worker_quality_value[j].scale;
            for(int i = 0; i < task_number; i++) {
                unsigned label = crh.answer_data[i][j];
                this->pb.val(worker_chosen_label_prediction[j][i][0][0]) = 0;
                this->pb.val(worker_chosen_label_prediction[j][i][0][1]) = 0;
                for(int k = 0; k < label_class_number; k++) {
                    this->pb.val(worker_chosen_label_prediction[j][i][k+1][0]) = this->pb.val(worker_chosen_label_prediction[j][i][k][0]) + (label == k ? task_label_prediction_value[i][k].value : 0);
                    this->pb.val(worker_chosen_label_prediction[j][i][k+1][1]) = this->pb.val(worker_chosen_label_prediction[j][i][k][1]) + (label == k ? task_label_prediction_value[i][k].scale : 0);                   
                }

                Float multiplied_value = task_label_prediction_value[i][label] * Float(1.0f / task_number);
                this->pb.val(multiply_result_variables_2[j][i][0]) = multiplied_value.value;
                this->pb.val(multiply_result_variables_2[j][i][1]) = multiplied_value.scale;
                multiply_gates_2[j][i]->generate_r1cs_witness(task_label_prediction_value[i][label], Float(1.0f / task_number), multiplied_value);

                Float added_value = worker_quality_value[j] + multiplied_value;
                this->pb.val(updated_worker_quality[j][i+1][0]) = added_value.value;
                this->pb.val(updated_worker_quality[j][i+1][1]) = added_value.scale;
                add_gates_2[j][i]->generate_r1cs_witness(worker_quality_value[j], multiplied_value, added_value);
                worker_quality_value[j] = added_value;
            }            
        }
    }
};

#endif