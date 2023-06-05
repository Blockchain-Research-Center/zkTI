#ifndef __zkTI_zencrowd_circuit_h
#define __zkTI_zencrowd_circuit_h

#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#include "../methods/Float.h"
#include "../methods/ZenCrowd.h"
#include "./gadgets.h"
#include "./gadgets2.h"
#include "../common.h"

using namespace libsnark;

#ifndef __BASE
#define __BASE 2
#endif

template <typename FieldT>
class ZenCrowdCircuit : public gadget<FieldT>
{

public:
    ZenCrowd &zc;
    std::vector<unsigned> &predicted_data;
    int task_number, worker_number, label_number;
    int label_class_number;

private:
    // (primary input) public input
    pb_variable<FieldT> base;
    std::vector<pb_variable<FieldT>> multiplies; // base**base**i

    // (auxiliary input) secret witness
    Float x,y;

    // (auxiliary input) auxiliary witness and gadgets
    FloatMultiplyGadget<FieldT> *test_multiplication;
    FloatDivideGadget<FieldT> *test_divide;
    FloatAddGadget<FieldT> *test_add;
    ExponentialGadget<FieldT> *test_exp;

    // helper function
    // initialize all the variables pb needed
    void init_pb_vars()
    {
        auto &prefix = this->annotation_prefix;

        base.allocate(this->pb, std::string("base"));
        init_one_dimension_vec(multiplies, 8, std::string("base_multiplies"));

        // set public number
        this->pb.set_input_sizes(8);
    }

    // assign values for public inputs
    void assign_public_inputs()
    {
        this->pb.val(base) = __BASE;
        for(int i = 0; i < 8; i++) {
            mpz_t base, tmp;
            mpz_init(tmp);
            mpz_init_set_ui(base, __BASE);
            mpz_pow_ui(tmp, base, (u64)pow(__BASE, i));
            this->pb.val(multiplies[i]) = libff::bigint<__limbs>(tmp);
        }

        // // predicted variables
        // for(int i = 0; i < task_number; i++) {
        //     this->pb.val(predicted_class_variables[i]) = predicted_data[i];
        // }   

        // // label class variables
        // for (int i = 0; i < label_class_number; i++)
        // {
        //     this->pb.val(label_class_variables[i]) = i;
        // }

        // // correct
        // for (int i = 0; i < task_number; i++)
        // {
        //     this->pb.val(correct[i]) = 1;
        // }
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
    ZenCrowdCircuit(protoboard<FieldT> &pb, ZenCrowd &_zc, std::vector<unsigned> &_predicted_data, const std::string &annotation = "") : gadget<FieldT>(pb, annotation), zc(_zc), predicted_data(_predicted_data)
    {
        task_number = zc.answer_data.size();
        worker_number = zc.answer_data[0].size();
        label_number = task_number * worker_number;
        label_class_number = zc.label_class_number;

        // initialize variables
        init_pb_vars();

        // assign public inputs
        assign_public_inputs();

        // initialize gadgets
        pb_variable<FieldT> x_v, y_v, x_s, y_s, r_s, r_v;
        x_v.allocate(this->pb, ""); 
        y_v.allocate(this->pb, ""); 
        x_s.allocate(this->pb, ""); 
        y_s.allocate(this->pb, ""); 
        r_s.allocate(this->pb, ""); 
        r_v.allocate(this->pb, ""); 

        // tc1
        // x = Float(3.213124124f);
        // y = Float(1.3f);
        // tc2
        // x = Float(0.0123f);
        // y = Float(0.00041323f);
        // tc3
        // x = Float(0.1f);
        // y = Float(0.000000002f);
        // tc4
        x = Float(0.2f);
        y = Float::zero();

        Float xy = x * y;
        // Float xdy = x / y;
        Float x1y = x + y;

        this->pb.val(x_v) = x.value;
        this->pb.val(x_s) = x.scale;
        this->pb.val(y_v) = y.value;
        this->pb.val(y_s) = y.scale;
        // add
        this->pb.val(r_s) = x1y.scale;
        this->pb.val(r_v) = x1y.value;
        test_add = new FloatAddGadget<FieldT>(this->pb, x_v, x_s, y_v, y_s, r_v, r_s, multiplies, "");
        // mul
        // this->pb.val(r_s) = xy.scale;
        // this->pb.val(r_v) = xy.value;
        // test_multiplication = new FloatMultiplyGadget<FieldT>(this->pb, x_v, x_s, y_v, y_s, r_v, r_s, "");
        // div
        // this->pb.val(r_s) = xdy.scale;
        // this->pb.val(r_v) = xdy.value;
        // test_divide = new FloatDivideGadget<FieldT>(this->pb, x_v, x_s, y_v, y_s, r_v, r_s, "");

        // pb_variable<FieldT> a, b;
        // a.allocate(this->pb, "");
        // b.allocate(this->pb, "");
        // this->pb.val(a) = 120;
        // mpz_t base, tmp;
        // mpz_init(tmp);
        // mpz_init_set_ui(base, __BASE);
        // mpz_pow_ui(tmp, base, 120);
        // this->pb.val(b) = libff::bigint<__limbs>(tmp);
        // test_exp = new ExponentialGadget<FieldT>(this->pb, a, b, multiplies);
    }

    ~ZenCrowdCircuit() {}

    void generate_r1cs_constraints()
    {
        // test_divide->generate_r1cs_constraints();
        // test_multiplication->generate_r1cs_constraints();
        test_add->generate_r1cs_constraints();
        // test_exp->generate_r1cs_constraints();
    }

    void generate_r1cs_witness()
    {
        Float xy = x * y;
        // Float x_y = x / y;
        Float x1y = x + y;
        std::cout << x << std::endl;
        std::cout << y << std::endl;
        std::cout << xy << std::endl;
        // std::cout << x_y << std::endl;
        std::cout << x1y << std::endl;

        // test_divide->generate_r1cs_witness(x, y, x_y);
        // test_multiplication->generate_r1cs_witness(x, y, xy);
        // test_exp->generate_r1cs_witness(120);
        test_add->generate_r1cs_witness(x, y, x1y);
    }
};

#endif