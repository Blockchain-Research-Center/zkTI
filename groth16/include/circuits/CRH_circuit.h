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
    std::vector<pb_variable<FieldT>> constant_variable;

    // (auxiliary input) secret witness
    std::vector<std::vector<pb_variable<FieldT>>> answer_variables;    std::vector<std::vector<pb_variable<FieldT>>> worker_quality;
    std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>> task_label_prediction;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> updated_task_label_prediction;
    std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>> add_result_variables;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> equality_result_variables;
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> distance_sum_variables; 
    std::vector<std::vector<pb_variable<FieldT>>> weight_sum_variables; 
    std::vector<std::vector<pb_variable<FieldT>>> updated_worker_quality;

    // (auxiliary input) auxiliary witness and gadgets
    std::vector<std::vector<std::vector<pb_variable<FieldT>>>> sum;
    // label equality gadgets
    std::vector<std::vector<std::vector<EqualityCheckGadget<FieldT>*>>> equality_check_gadgets; // task * worker * label
    // add gates 1
    std::vector<std::vector<std::vector<FloatAddGadget<FieldT>*>>> add_gates_1;
    // add gates 2
    std::vector<std::vector<FloatAddGadget<FieldT>*>> add_gates_2;
    // divide gates 1
    std::vector<std::vector<FloatDivideGadget<FieldT>*>> divide_gates_1;
    // add gates 3
    std::vector<std::vector<FloatAddGadget<FieldT>*>> add_gates_3;
    // add gates 4
    std::vector<FloatAddGadget<FieldT>*> add_gates_4;
    // divide gates 2
    std::vector<FloatDivideGadget<FieldT>*> divide_gates_2;

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
        init_one_dimension_vec(this->pb, constant_variable, 2, std::string("constant_variable"));

        init_two_dimension_vec(this->pb, worker_quality, worker_number, 2, "worker_quality");
        init_two_dimension_vec(this->pb, updated_worker_quality, worker_number, 2, "updated_worker_quality");
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
        // distance sum variables
        init_three_dimension_vec(this->pb, distance_sum_variables, worker_number, task_number + 1, 2, std::string("distance_sum_variables"));
        // weight sum variables
        init_two_dimension_vec(this->pb, weight_sum_variables, worker_number + 1, 2, std::string("weight_sum_variables"));

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
        }

        // task_label_prediction
        for(int i = 0; i < task_number; i++) {
            for(int k = 0; k < label_class_number; k++) {
                // 1
                this->pb.val(task_label_prediction[i][k][0][0]) = Float::zero().value;
                this->pb.val(task_label_prediction[i][k][0][1]) = Float::zero().scale;
            }
        }

        // label class variables
        for (int i = 0; i < label_class_number; i++)
        {
            this->pb.val(label_class_variables[i]) = i;
        }

        this->pb.val(constant_variable[0]) = Float::one().value;
        this->pb.val(constant_variable[1]) = Float::one().scale;
    }

public:
    CRHCircuit(protoboard<FieldT> &pb, CRH &_crh, std::vector<unsigned> &_predicted_data, int _task_number, int _worker_number,  std::string &annotation = "") : gadget<FieldT>(pb, annotation), crh(_crh), predicted_data(_predicted_data)
    {
        assert(_task_number <= crh.answer_data.size());
        assert(_worker_number <= crh.answer_data[0].size());
        task_number = _task_number;
        worker_number = _worker_number;
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

        // initialize add gadgets 1
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<std::vector<FloatAddGadget<FieldT>*>> row_x;
            for (int j = 0; j < worker_number; ++j)
            {
                std::vector<FloatAddGadget<FieldT>*> row_y;
                for (int k = 0; k < label_class_number; k++) {
                    auto gadget_name = annotation + std::string("add_gates_1_") + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k);   
                    row_y.push_back(new FloatAddGadget<FieldT>(this->pb, task_label_prediction[i][k][j][0], task_label_prediction[i][k][j][1], worker_quality[j][0], worker_quality[j][1], add_result_variables[i][k][j][0], add_result_variables[i][k][j][1], multiplies, random_point, random_coefficients, gadget_name));
                }
                row_x.push_back(row_y);
            }
            add_gates_1.push_back(row_x);
        }

        // initialize add gadgets 2
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<FloatAddGadget<FieldT>*> row_x;
            for (int k = 0; k < label_class_number; ++k)
            {
                auto gadget_name = annotation + std::string("add_gates_2_") + std::to_string(i) + "_" + std::to_string(k);   
                // we leave the implementation of permutation check inside gadget
                row_x.push_back(new FloatAddGadget<FieldT>(this->pb, sum[i][k][0], sum[i][k][1], task_label_prediction[i][k][worker_number][0], task_label_prediction[i][k][worker_number][1], sum[i][k+1][0], sum[i][k+1][1], multiplies, random_point, random_coefficients, gadget_name));      
            }
            add_gates_2.push_back(row_x);
        }

        // initialize divide gadgets 1
        for (int i = 0; i < task_number; ++i)
        {
            std::vector<FloatDivideGadget<FieldT>*> row_x;
            for (int k = 0; k < label_class_number; ++k)
            {
                auto gadget_name = annotation + std::string("divide_gates_1_") + std::to_string(i) + "_" + std::to_string(k);   
                row_x.push_back(new FloatDivideGadget<FieldT>(this->pb, task_label_prediction[i][k][worker_number][0], task_label_prediction[i][k][worker_number][1], sum[i][label_class_number][0], sum[i][label_class_number][1], updated_task_label_prediction[i][k][0], updated_task_label_prediction[i][k][1], gadget_name));      
            }
            divide_gates_1.push_back(row_x);
        }

        // initialize add gadgets 3
        for (int j = 0; j < worker_number; ++j)
        {
            std::vector<FloatAddGadget<FieldT>*> row_x;
            for (int i = 0; i < task_number; ++i)
            {
                auto gadget_name = annotation + std::string("add_gates_3_") + std::to_string(j) + "_" + std::to_string(i);   
                // we leave the implementation of permutation check inside gadget
                row_x.push_back(new FloatAddGadget<FieldT>(this->pb, distance_sum_variables[j][i][0], distance_sum_variables[j][i][1], constant_variable[0], constant_variable[1], distance_sum_variables[j][i+1][0], distance_sum_variables[j][i+1][1], multiplies, random_point, random_coefficients, gadget_name));      
            }
            add_gates_3.push_back(row_x);
        }

        // initialize add gadgets 4
        for (int j = 0; j < worker_number; ++j)
        {
            auto gadget_name = annotation + std::string("add_gates_4_") + std::to_string(j);
            add_gates_4.push_back(new FloatAddGadget<FieldT>(this->pb, weight_sum_variables[j][0], weight_sum_variables[j][1], distance_sum_variables[j][task_number][0],  distance_sum_variables[j][task_number][1], weight_sum_variables[j+1][0], weight_sum_variables[j+1][1], multiplies, random_point, random_coefficients, gadget_name));
        }

        // initialize divide gadgets 2
        for (int j = 0; j < worker_number; ++j)
        {
            auto gadget_name = annotation + std::string("divide_gates_2_") + std::to_string(j);   
            divide_gates_2.push_back(new FloatDivideGadget<FieldT>(this->pb, weight_sum_variables[worker_number][0], weight_sum_variables[worker_number][1], distance_sum_variables[j][task_number][0], distance_sum_variables[j][task_number][1], updated_worker_quality[j][0], updated_worker_quality[j][1], gadget_name));
        }
    }

    ~CRHCircuit() {}

    void generate_r1cs_constraints()
    {
        for(int i = 0; i < task_number; i++) {
            for(int j = 0; j < worker_number; j++) {
                for(int k = 0; k < label_class_number; k++) {
                    add_gates_1[i][j][k]->generate_r1cs_constraints();
                    equality_check_gadgets[i][j][k]->generate_r1cs_constraints();
                    // r*(s_i+1-a) = (1-r)*(s_i+1-s_i)
                    add_r1cs(equality_result_variables[i][j][k], 2 * task_label_prediction[i][k][j+1][0] - add_result_variables[i][k][j][0] - task_label_prediction[i][k][j][0], task_label_prediction[i][k][j+1][0] - task_label_prediction[i][k][j][0]);
                    add_r1cs(equality_result_variables[i][j][k], 2 * task_label_prediction[i][k][j+1][1] - add_result_variables[i][k][j][1] - task_label_prediction[i][k][j][1], task_label_prediction[i][k][j+1][1] - task_label_prediction[i][k][j][1]);
                }
            }

            for(int k = 0; k < label_class_number; k++) {
                add_gates_2[i][k]->generate_r1cs_constraints();
                divide_gates_1[i][k]->generate_r1cs_constraints();
            }
        }

        for(int j = 0; j < worker_number; j++) {
            for(int i = 0; i < task_number; i++) {
                add_gates_3[j][i]->generate_r1cs_constraints();
            }
            // Here we omit the distance function for simplity.
            add_gates_4[j]->generate_r1cs_constraints();
            divide_gates_2[j]->generate_r1cs_constraints();
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
                task_label_prediction_value[i][k] = Float::zero();
            }

            for(int j = 0; j < worker_number; j++) {
                unsigned label = crh.answer_data[i][j];
                for(int k = 0; k < label_class_number; k++) {
                    Float a = worker_quality_value[j];
                    Float add_value = task_label_prediction_value[i][k] + a;
                    
                    add_gates_1[i][j][k]->generate_r1cs_witness(task_label_prediction_value[i][k], a, add_value);

                    this->pb.val(add_result_variables[i][k][j][0]) = add_value.value;
                    this->pb.val(add_result_variables[i][k][j][1]) = add_value.scale;

                    if(label == k) {
                        task_label_prediction_value[i][k] = add_value;
                    } else {
                        task_label_prediction_value[i][k] = task_label_prediction_value[i][k];
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
                divide_gates_1[i][k]->generate_r1cs_witness(task_label_prediction_value[i][k], sums_value, divided_value);
                task_label_prediction_value[i][k] = divided_value;
            }
        }

        Float weight_sum_value = Float::zero();
        std::vector<Float*> distance_sum_values;
        this->pb.val(weight_sum_variables[0][0]) = weight_sum_value.value;
        this->pb.val(weight_sum_variables[0][1]) = weight_sum_value.scale;
        for(int j = 0; j < worker_number; j++) {
            Float distance_sum_value = Float::zero();
            this->pb.val(distance_sum_variables[j][0][0]) = distance_sum_value.value;
            this->pb.val(distance_sum_variables[j][0][1]) = distance_sum_value.scale;
            for(int i = 0; i < task_number; i++) {
                Float constant_value = Float::one();
                Float add_value = distance_sum_value + constant_value;
                this->pb.val(distance_sum_variables[j][i+1][0]) = add_value.value;
                this->pb.val(distance_sum_variables[j][i+1][1]) = add_value.scale;
                add_gates_3[j][i]->generate_r1cs_witness(distance_sum_value, constant_value,add_value);
                distance_sum_value = add_value;
            }    
            distance_sum_values.push_back(new Float());
            *distance_sum_values[j] = distance_sum_value;

            Float added_value =  weight_sum_value + distance_sum_value;
            this->pb.val(weight_sum_variables[j+1][0]) = added_value.value;
            this->pb.val(weight_sum_variables[j+1][1]) = added_value.scale;
            add_gates_4[j]->generate_r1cs_witness(weight_sum_value, distance_sum_value, added_value);
            weight_sum_value = added_value;
        }

        for(int j = 0; j < worker_number; j++) {
            Float divide_value = weight_sum_value / *distance_sum_values[j];
            this->pb.val(updated_worker_quality[j][0]) = divide_value.value;
            this->pb.val(updated_worker_quality[j][1]) = divide_value.scale;
            divide_gates_2[j]->generate_r1cs_witness(weight_sum_value, *distance_sum_values[j], divide_value);
        }
    }
};

#endif