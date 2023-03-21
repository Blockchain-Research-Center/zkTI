#ifndef __zkTI_gadgets_h
#define __zkTI_gadgets_h

#include <common.h>

template <typename FieldT>
class ComparisonGadget: public gadget<FieldT> {

public:
    const pb_variable<FieldT> &comparison_result;
    const pb_variable<FieldT> &x, &y;
    const pb_variable<FieldT> &diff;

    ComparisonGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x_, pb_variable <FieldT> &y_,
                     pb_variable <FieldT> &comparison_result_, pb_variable <FieldT> &diff_,
                     const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation), x(x_), y(y_), comparison_result(comparison_result_), diff(diff_) {
    }

    ~ComparisonGadget() {}

    void generate_r1cs_constraints() {
        // 1 y=x+d 0 x=
        add_r1cs(2 * comparison_result - 1, y - x, diff);
    }

    void generate_r1cs_witness() {}

};

// assume big-endian
template<typename FieldT>
class DecompositionCheckGadget : public gadget<FieldT> {
private:
    pb_variable <FieldT> *vars, *decompositions;
    int n_vars, n_bit_per_var;
public:
    DecompositionCheckGadget(protoboard <FieldT> &pb, pb_variable <FieldT> *vars_,
                             pb_variable <FieldT> *decompositions_,
                             int n_vars_, int n_bit_per_var_, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation) {
        vars = vars_;
        decompositions = decompositions_;
        n_vars = n_vars_;
        n_bit_per_var = n_bit_per_var_;
    }

    void generate_r1cs_constraints() {
        // check summation
        for (int i = 0; i < n_vars; ++i) {
            int base = i * n_bit_per_var;
            auto sum = linear_combination<FieldT>(decompositions[base + n_bit_per_var - 1]);
            for (int j = n_bit_per_var - 2; j >= 0; --j) {
                sum = sum + decompositions[base + j] * FieldT(1ULL << (n_bit_per_var - j - 1));
            }
            add_r1cs(sum, 1, vars[i]);
        }

        // check 0 or 1
        for (int i = 0; i < n_vars * n_bit_per_var; ++i) {
            add_r1cs(decompositions[i], 1 - decompositions[i], 0);
        }
    }

    void generate_r1cs_witness() {
        // nothing to do here assuming decomposition is calculated outside, maybe it's better to put it here
        return;
    }
};

template<typename FieldT>
class ArgmaxGadget : public gadget<FieldT> {
private:
    int n_bits_per_number;
    int n_vars;
    pb_variable <FieldT> *values;
    pb_variable <FieldT> *max_index; // one-hot

    pb_variable <FieldT> max_value;
    pb_variable <FieldT> *max_value_decompositions;
    pb_variable <FieldT> *values_times_index;
    pb_variable <FieldT> *value_decompositions;

    pb_variable <FieldT> *place_holder;

    DecompositionCheckGadget<FieldT> *decompositionCheckGadget;
    DecompositionCheckGadget<FieldT> *decompositionCheckGadget1;

public:
    ArgmaxGadget(protoboard <FieldT> &pb, pb_variable <FieldT> *values_, pb_variable <FieldT> *max_index_,
                 int n_vars_, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation) {
        n_vars = n_vars_;
        values = values_;
        max_index = max_index_;
        n_bits_per_number = log(n_vars) + 1;

        _init_pb_array(pb, value_decompositions, n_vars * n_bits_per_number,
                       annotation + std::string("value_decompositions"));
        _init_pb_array(pb, values_times_index, n_vars, annotation + std::string("values_times_index"));

        max_value.allocate(pb, annotation + std::string("max_value"));
        _init_pb_array(pb, max_value_decompositions, n_bits_per_number, std::string("max_value_decompositions"));

        _init_pb_array(pb, place_holder, n_vars * n_bits_per_number * 2, std::string("placeholder"));

        decompositionCheckGadget = new DecompositionCheckGadget<FieldT>(pb, values_, value_decompositions, n_vars,
                                                                        n_bits_per_number,
                                                                        annotation + std::string("argmax_gadget"));
        decompositionCheckGadget1 = new DecompositionCheckGadget<FieldT>(pb, &max_value, max_value_decompositions, 1,
                                                                         n_bits_per_number,
                                                                         annotation + std::string("argmax_gadget1"));
    }

    void generate_r1cs_constraints() {
        decompositionCheckGadget->generate_r1cs_constraints();
        decompositionCheckGadget1->generate_r1cs_constraints();

        auto sum = linear_combination<FieldT>();
        for (int i = 0; i < n_vars; ++i) {
            add_r1cs(max_index[i], 1 - max_index[i], 0);
            sum = sum + max_index[i];
        }
        add_r1cs(sum, 1, 1);

        sum = linear_combination<FieldT>();
        for (int i = 0; i < n_vars; ++i) {
            add_r1cs(values[i], max_index[i], values_times_index[i]);
            sum = sum + values_times_index[i];
        }
        add_r1cs(max_value, 1, sum);

        // check less than or equal
        for (int i = 0; i < n_vars; ++i) {
            // less than or equal(values[i], max_value);
            for (int j = 0; j < n_bits_per_number; ++j) {
                auto &x = place_holder[i * n_bits_per_number + j], &y = place_holder[(i + n_vars) * n_bits_per_number +
                                                                                     j];
                add_r1cs(x, 1, x); // placeholder
                add_r1cs(y, 1, y);
            }
        }

    }

    void generate_r1cs_witness() {
        decompositionCheckGadget->generate_r1cs_witness();
        decompositionCheckGadget1->generate_r1cs_witness();
        return;
    }
};

template<typename FieldT>
class EqualityCheckGadget : public gadget<FieldT> {
private:
    int n_bits;
    pb_variable <FieldT> &x, &y, &result;
    pb_variable <FieldT> *x_dec, *y_dec; // bit decomposition
    pb_variable <FieldT> *bit_equal; // equality check of each bit
    pb_variable <FieldT> *aggr_equal; // aggregate the result from bit equal

    DecompositionCheckGadget<FieldT> *decompositionCheckGadgets;

public:
    EqualityCheckGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x_, pb_variable <FieldT> &y_,
                        pb_variable <FieldT> &result_,
                        int n_bits_ = 32, const std::string &annotation = "") :
            gadget<FieldT>(pb, annotation), n_bits(n_bits_), x(x_), y(y_), result(result_) {

        decompositionCheckGadgets = (DecompositionCheckGadget<FieldT> *) malloc(sizeof(DecompositionCheckGadget<FieldT>) * 2);

        _init_pb_array(this->pb, x_dec, n_bits, annotation + "x_dec");
        _init_pb_array(this->pb, y_dec, n_bits, annotation + "y_dec");
        _init_pb_array(this->pb, bit_equal, n_bits, annotation + "bit_equal");
        _init_pb_array(this->pb, aggr_equal, n_bits - 1, annotation + "aggr_equal");

        new(decompositionCheckGadgets + 0) DecompositionCheckGadget<FieldT>(pb, &x, x_dec, 1, n_bits, annotation + "decompositionCheckGadgets0");
        new(decompositionCheckGadgets + 1) DecompositionCheckGadget<FieldT>(pb, &y, y_dec, 1, n_bits, annotation + "decompositionCheckGadgets1");
    }

    void generate_r1cs_constraints() {
        decompositionCheckGadgets[0].generate_r1cs_constraints();
        decompositionCheckGadgets[1].generate_r1cs_constraints();
        for (int i = 0; i < n_bits; ++i) {
            add_r1cs(2 * x_dec[i], y_dec[i], x_dec[i] + y_dec[i] + bit_equal[i] - 1); // z_i = x_i == y_i
        }
        add_r1cs(bit_equal[0], bit_equal[1], aggr_equal[0]);
        for (int i = 1; i < n_bits - 1; ++i) {
            add_r1cs(aggr_equal[i - 1], bit_equal[i + 1], aggr_equal[i]);
        }
        add_r1cs(aggr_equal[n_bits - 2], 1, result);
    }

    void generate_r1cs_witness(unsigned x, unsigned y) {
        for (int i = 0; i < n_bits; ++i) {
            pb_eval(x_dec[i]) = (x >> (n_bits - 1 - i)) & 1U;
            pb_eval(y_dec[i]) = (y >> (n_bits - 1 - i)) & 1U;
            pb_eval(bit_equal[i]) = pb_eval(x_dec[i]) == pb_eval(y_dec[i]);
        }
        pb_eval(aggr_equal[0]) = pb_eval(bit_equal[0]) * pb_eval(bit_equal[1]);
        for (int i = 1; i < n_bits - 1; ++i) {
            pb_eval(aggr_equal[i]) = pb_eval(aggr_equal[i - 1]) * pb_eval(bit_equal[i + 1]);
        }
    }

    ~EqualityCheckGadget() {

    }
};

template<typename FieldT>
class MajorityGadget : public gadget<FieldT> {
private:
    int n_vars;
    pb_variable <FieldT> *values, value_majority;

    pb_variable <FieldT> *equal; // equal[i][j] = 1 iff values[i] == values[j]
    pb_variable <FieldT> *counts; // counts how many times a class appear
    pb_variable <FieldT> *max_index; // one-hot representation of the position of the max count

    ArgmaxGadget<FieldT> *argmaxGadget;

public:
    MajorityGadget(protoboard <FieldT> &pb, pb_variable <FieldT> *values_, pb_variable <FieldT> value_majority_,
                   int n_vars_, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation) {
        n_vars = n_vars_;
        values = values_;
        value_majority = value_majority_;

        _init_pb_array(pb, counts, n_vars, annotation + std::string("counts"));
        _init_pb_array(pb, max_index, n_vars, annotation + std::string("max_index"));
        _init_pb_array(pb, equal, n_vars * n_vars / 2, annotation + std::string("equal"));

        argmaxGadget = new ArgmaxGadget<FieldT>(pb, values, max_index, n_vars,
                                                annotation + std::string("argmaxGadget"));
    }

    void generate_r1cs_constraints() {
        for (int i = 0; i < n_vars; ++i) {
            auto sum = linear_combination<FieldT>();
            for (int j = 0; j < n_vars; ++j) {
                if (i < j) {
                    add_r1cs(equal[j], values[i] - values[j], 0);
                }
                sum = sum + equal[j];
            }
            add_r1cs(counts[i], 1, sum);
        }

        argmaxGadget->generate_r1cs_constraints();

    }

    void generate_r1cs_witness() {
        // nothing to do here
        return;
    }

};

#endif