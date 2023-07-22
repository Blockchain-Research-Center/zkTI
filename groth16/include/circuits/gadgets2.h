#ifndef __zkTI_gadgets2_h
#define __zkTI_gadgets2_h

#include <cmath>

#include "../methods/Float.h"
#include "./gadgets.h"
#include "../common.h"

#ifndef __float_constant_
#define __BASE 2
#define __MAXIMUM_BIT 23
#define __MINIMUM_BASE_POW 256
#endif

#define __scale_n_bit 6
#define __accuracy (__MAXIMUM_BIT - 1)
#define __accuracy_bound (__accuracy + 2)
#define __max_delta_bit 8
// #define __max_delta_bit 7
#define __omega __MAXIMUM_BIT // omega is the bit-n of each number

// This gadget is used when x, y is not constrainted.
// and we compare x, y in the field directly.
template <typename FieldT>
class LessThanGadget: public gadget<FieldT> {

public:
    const pb_variable<FieldT> x, y, result;

    pb_variable<FieldT> num;
    pb_variable<FieldT> *num_decomposed;

    DecompositionCheckGadget<FieldT> *num_decomposition;

    int n_bit;
    libff::bigint<__limbs> ceil;
    mpz_t value;

    LessThanGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x_, pb_variable <FieldT> &y_, pb_variable<FieldT> &result_, int n_bit, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation), x(x_), y(y_), result(result_), n_bit(n_bit) {
        _init_pb_array(pb, num_decomposed, n_bit + 1, annotation + std::string("num_decomposed"));
        num.allocate(pb, annotation + std::string("num"));

        num_decomposition = new DecompositionCheckGadget<FieldT>(pb, &num, num_decomposed, 1, n_bit + 1, annotation + std::string("num_decomposition"));

        mpz_init_set_d(value, pow(2, n_bit));
        ceil = libff::bigint<__limbs>(value);
    }

    ~LessThanGadget() {}

    void generate_r1cs_constraints() {
        // calculate x + 2^n -y 
        add_r1cs(x - y + FieldT(ceil) - num, 1, 0);
        // decompose
        num_decomposition->generate_r1cs_constraints();
        // equal
        add_r1cs(result - num_decomposed[0], 1, 0);
    }

    void generate_r1cs_witness(u64 x, u64 y) {
        // num
        mpz_t x_value, y_value, num_value;
        mpz_init(num_value);
        mpz_init_set_ui(x_value, x);
        mpz_init_set_ui(y_value, y);
        mpz_add(x_value, x_value, value);
        mpz_sub(num_value, x_value, y_value);
        this->pb.val(num) = libff::bigint<__limbs>(num_value);

        // num decomposition
        mpz_t tmp;
        mpz_init(tmp);
        for (int i = 0; i < n_bit + 1; ++i) {
            mpz_fdiv_q_2exp(tmp, num_value, n_bit - i);
            u64 tmp_value = mpz_get_ui(tmp);
            this->pb.val(num_decomposed[i]) = tmp_value & 1U;
        }
        num_decomposition->generate_r1cs_witness();
    }

    void generate_r1cs_witness(mpz_t x, mpz_t y) {
        // num
        mpz_t x_value, y_value, num_value;
        mpz_init(num_value);
        mpz_init(x_value);
        mpz_init(y_value);
        mpz_set(x_value, x);
        mpz_set(y_value, y);
        mpz_add(x_value, x_value, value);
        mpz_sub(num_value, x_value, y_value);
        this->pb.val(num) = libff::bigint<__limbs>(num_value);

        // num decomposition
        mpz_t tmp;
        mpz_init(tmp);
        for (int i = 0; i < n_bit + 1; ++i) {
            mpz_fdiv_q_2exp(tmp, num_value, n_bit - i);
            u64 tmp_value = mpz_get_ui(tmp);
            this->pb.val(num_decomposed[i]) = tmp_value & 1U;
        }
        num_decomposition->generate_r1cs_witness();
    }

};

// y = base^x 
// x is at most N which means we use logN+1 bit to represent x
template <typename FieldT>
class ExponentialGadget: public gadget<FieldT> {

public:
    int n_bit;
    pb_variable<FieldT> x, y;
    pb_variable<FieldT> *x_decomposed;
    pb_variable<FieldT> *multiplied;
    pb_variable<FieldT> *intermediates;

    DecompositionCheckGadget<FieldT> *x_decomposition;
    
    const std::vector<pb_variable<FieldT>> &multiplies;

    ExponentialGadget() {

    }

    ExponentialGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x, pb_variable <FieldT> &y, std::vector<pb_variable<FieldT>> &multiplies, int n_bit = 8, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation), x(x), y(y), n_bit(n_bit), multiplies(multiplies) {
        // init pb val
        _init_pb_array(pb, x_decomposed, n_bit,
                       annotation + std::string("x_decomposed"));
        _init_pb_array(pb, multiplied, n_bit, annotation + std::string("after_multiply_values"));
        _init_pb_array(pb, intermediates, n_bit + 1, annotation + std::string("intermediates"));

        x_decomposition = new DecompositionCheckGadget<FieldT>(pb, &x, x_decomposed, 1, n_bit, annotation + std::string("decomposition_gadget"));
    }
 
    ~ExponentialGadget() {
         
    }

    void generate_r1cs_constraints() {
        // x decomposition
        x_decomposition->generate_r1cs_constraints();
        // s[0] = 1
        add_r1cs(intermediates[0] - 1, 1, 0);
        for(int i = 0; i < n_bit; i++) {
            // mul[i] = s[i] * BASE**BASE**i
            add_r1cs(intermediates[i], multiplies[i], multiplied[i]);
            // s_i+1 = b * mul + (1-b) * s_i
            add_r1cs(x_decomposed[n_bit - 1 - i], multiplied[i] - intermediates[i], intermediates[i+1] - intermediates[i]);
        }
        // check equal  
        add_r1cs(intermediates[n_bit], 1, y);
    }

    void generate_r1cs_witness(i16 value) {
        for (int i = 0; i < n_bit; ++i) {
            this->pb.val(x_decomposed[i]) = (value >> (n_bit - i - 1)) & 1U;
        }
        x_decomposition->generate_r1cs_witness();

        // intermediate[0]
        this->pb.val(intermediates[0]) = 1;
        mpz_t intermediate_value;
        mpz_init_set_ui(intermediate_value, 1);
        for(int i = 0; i < n_bit; ++i) {
            // multiply values
            mpz_t base, multiply_value, multiplied_value; 
            mpz_init(multiply_value);
            mpz_init(multiplied_value);
            mpz_init_set_ui(base, __BASE);
            mpz_pow_ui(multiply_value, base, (u64)pow(__BASE, i));
            // multiplied values
            mpz_mul(multiplied_value, intermediate_value, multiply_value);
            this->pb.val(multiplied[i]) = libff::bigint<__limbs>(multiplied_value);
            
            // intermediates
            if((value >> i) & 1U) {
                // 1
                this->pb.val(intermediates[i + 1]) = libff::bigint<__limbs>(multiplied_value);
                mpz_set(intermediate_value, multiplied_value);
            } else {
                // 0
                this->pb.val(intermediates[i + 1]) = libff::bigint<__limbs>(intermediate_value);
            }
        }
    }

};

template <typename FieldT>
class FloatMultiplyGadget: public gadget<FieldT> {

public:
    pb_variable<FieldT> result_v, result_s;
    pb_variable<FieldT> x_v, x_s, y_s, y_v;

    pb_variable<FieldT> theta;
    pb_variable<FieldT> lambda; // lambda is the exp of theta, lambda = BASE**theta
    pb_variable<FieldT> a; // a = x_v * y_v
    pb_variable<FieldT> b; // b = a - result_v * lambda
    pb_variable<FieldT> c; // c = sigma * b

    pb_variable<FieldT> *result_v_decomposed;
    DecompositionCheckGadget<FieldT> *result_v_decomposition;

    u64 exp_base_omega_sub_1 = pow(__BASE, __omega - 1); 
    u64 sigma = pow(__BASE, __accuracy); // actually it should be sigma^-1

    pb_variable<FieldT> lessthan_result;
    LessThanGadget<FieldT> *lessthan_checker;

    FloatMultiplyGadget() {

    }

    FloatMultiplyGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x_v, pb_variable <FieldT> &x_s, pb_variable <FieldT> &y_v, pb_variable <FieldT> &y_s, pb_variable <FieldT> &result_v, pb_variable <FieldT> &result_s, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation), x_v(x_v), x_s(x_s), y_v(y_v), y_s(y_s), result_v(result_v), result_s(result_s) {
        // init pb val
        _init_pb_array(pb, result_v_decomposed, __omega, annotation + std::string("result_v_decomposed"));
        theta.allocate(pb, annotation + std::string("theta"));
        lambda.allocate(pb, annotation + std::string("lambda"));
        a.allocate(pb, annotation + std::string("a"));
        b.allocate(pb, annotation + std::string("b"));
        c.allocate(pb, annotation + std::string("c"));
        lessthan_result.allocate(pb, annotation + std::string("lessthan_result"));
        result_v_decomposition = new DecompositionCheckGadget<FieldT>(pb, &result_v, result_v_decomposed, 1, __omega, annotation + std::string("result_v_decomposition"));
        lessthan_checker = new LessThanGadget<FieldT>(pb, c, a, lessthan_result, 2 * __omega, annotation + std::string("lessthan_checker"));
    }

    ~FloatMultiplyGadget() {
         
    }

    void generate_r1cs_constraints() {
        // result_v belong to [2^(w-1), 2^w)
        result_v_decomposition->generate_r1cs_constraints();
        add_r1cs(result_v_decomposed[0], 1, 1); // most-significant bit is 1
        // result_s + theta = x_s + y_s
        add_r1cs(result_s + theta, 1, x_s + y_s);
        // theta is w or w-1
        add_r1cs(theta - __omega + 1, theta - __omega, 0);
        // lambda = (theta - w + 1) * BASE ** theta - (theta - w) * BASE ** (theta - 1) 
        add_r1cs(lambda, 1, exp_base_omega_sub_1 * (theta - __omega + 2));
        // a = x_v * y_v
        add_r1cs(x_v, y_v, a);
        // b = x_v * y_v - result_v * lambda
        add_r1cs(result_v, lambda, a - b);
        // c = sigma * b
        add_r1cs(sigma * b, 1, c);
        // c < a
        lessthan_checker->generate_r1cs_constraints();
        add_r1cs(lessthan_result, 1, 0);
    }

    void generate_r1cs_witness(Float x, Float y, Float result) {
        // result_v decomposed
        for (int i = 0; i < __omega; ++i) {
            this->pb.val(result_v_decomposed[i]) = (result.value >> (__omega - i - 1)) & 1U;
        }
        result_v_decomposition->generate_r1cs_witness();
        // a
        u64 a_value = x.value * y.value;
        this->pb.val(a) = libff::bigint<__limbs>(a_value);
        // theta
        int shift = 0; // shift time
        u64 value = a_value;
        while(value > Float::get_maximum_value()) {
            value /= __BASE;
            shift++;
        }
        this->pb.val(theta) = shift;
        // lambda
        u64 lambda_value = (u64)pow(__BASE, shift);
        this->pb.val(lambda) = libff::bigint<__limbs>(lambda_value);
        // b
        u64 b_value = a_value - result.value * lambda_value;
        this->pb.val(b) = libff::bigint<__limbs>(b_value);
        // c
        u64 c_value = (u64)(sigma) *  b_value;
        this->pb.val(c) = libff::bigint<__limbs>(c_value);
        // less than
        // std::cout << c_value << " " << a_value << " " << a_value - c_value << std::endl;
        lessthan_checker->generate_r1cs_witness(c_value, a_value);
        this->pb.val(lessthan_result) = 0;
    }

};

template <typename FieldT>
class FloatDivideGadget: public gadget<FieldT> {

public:
    pb_variable<FieldT> result_v, result_s;
    pb_variable<FieldT> x_v, x_s, y_s, y_v;

    pb_variable<FieldT> theta;
    pb_variable<FieldT> lambda; // lambda is the exp of theta, lambda = BASE**theta
    pb_variable<FieldT> a; // a = result_v * y_v
    pb_variable<FieldT> b; // b = d - a * lambda
    pb_variable<FieldT> c; // c = sigma * b
    pb_variable<FieldT> d; // d = x_v * BASE**omega

    pb_variable<FieldT> *result_v_decomposed;
    DecompositionCheckGadget<FieldT> *result_v_decomposition;

    u64 exp_base_omega = pow(__BASE, __omega) - 1; // use 1 to avoid overflow
    u64 sigma = pow(__BASE, __accuracy); // actually it should be sigma^-1

    pb_variable<FieldT> lessthan_result;
    LessThanGadget<FieldT> *lessthan_checker;

    FloatDivideGadget() {

    }

    FloatDivideGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x_v, pb_variable <FieldT> &x_s, pb_variable <FieldT> &y_v, pb_variable <FieldT> &y_s, pb_variable <FieldT> &result_v, pb_variable <FieldT> &result_s, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation), x_v(x_v), x_s(x_s), y_v(y_v), y_s(y_s), result_v(result_v), result_s(result_s) {
        // init pb val
        _init_pb_array(pb, result_v_decomposed, __omega, annotation + std::string("result_v_decomposed"));
        theta.allocate(pb, annotation + std::string("theta"));
        lambda.allocate(pb, annotation + std::string("lambda"));
        a.allocate(pb, annotation + std::string("a"));
        b.allocate(pb, annotation + std::string("b"));
        c.allocate(pb, annotation + std::string("c"));
        d.allocate(pb, annotation + std::string("d"));
        lessthan_result.allocate(pb, annotation + std::string("lessthan_result"));
        result_v_decomposition = new DecompositionCheckGadget<FieldT>(pb, &result_v, result_v_decomposed, 1, __omega, annotation + std::string("result_v_decomposition"));
        lessthan_checker = new LessThanGadget<FieldT>(pb, c, d, lessthan_result, 2 * __omega, annotation + std::string("lessthan_checker"));
    }

    ~FloatDivideGadget() {
         
    }

    void generate_r1cs_constraints() {
        // result_v belong to [2^(w-1), 2^w)
        result_v_decomposition->generate_r1cs_constraints();
        add_r1cs(result_v_decomposed[0], 1, 1); // most-significant bit is 1
        // result_s + theta = x_s + omega - y_s
        add_r1cs(result_s + theta, 1, x_s + __omega - y_s);
        // theta is 1 or 0
        add_r1cs(theta - 1, theta, 0);
        // lambda = theta + 1
        add_r1cs(lambda, 1, theta + 1);
        // a = result_v * y_v
        add_r1cs(result_v, y_v, a);
        // b = d - a * lambda
        add_r1cs(a, lambda, d - b);
        // c = sigma * b
        add_r1cs(sigma * b, 1, c);
        // d = x_v * BASE**omega
        add_r1cs(x_v * exp_base_omega, 1, d);
        // c < d
        lessthan_checker->generate_r1cs_constraints();
        add_r1cs(lessthan_result, 1, 0);
    }

    void generate_r1cs_witness(Float x, Float y, Float result) {
        // result_v decomposed
        for (int i = 0; i < __omega; ++i) {
            this->pb.val(result_v_decomposed[i]) = (result.value >> (__omega - i - 1)) & 1U;
        }
        result_v_decomposition->generate_r1cs_witness();
        // a
        u64 a_value = result.value * y.value;
        this->pb.val(a) = libff::bigint<__limbs>(a_value);
        // d
        u64 d_value = x.value * Float::get_maximum_value();
        this->pb.val(d) = libff::bigint<__limbs>(d_value);
        // theta
        int shift = 0; // shift time
        u64 value = d_value / y.value;
        while(value > Float::get_maximum_value()) {
            value /= __BASE;
            shift++;
        }
        this->pb.val(theta) = shift;
        // lambda
        u64 lambda_value = (u64)pow(__BASE, shift);
        this->pb.val(lambda) = libff::bigint<__limbs>(lambda_value);
        // b
        u64 b_value = d_value - a_value * lambda_value;
        this->pb.val(b) = libff::bigint<__limbs>(b_value);
        // c
        u64 c_value = (u64)(sigma) *  b_value;
        this->pb.val(c) = libff::bigint<__limbs>(c_value);
        // less than
        lessthan_checker->generate_r1cs_witness(c_value, d_value);
        // less than result
        this->pb.val(lessthan_result) = 0;
    }

};

template <typename FieldT>
class FloatAddGadget: public gadget<FieldT> {

public:
    pb_variable<FieldT> result_v, result_s;
    pb_variable<FieldT> x_v, x_s, y_s, y_v;
    pb_variable<FieldT> *a_1, *a_2, *b_1, *b_2;

    pb_variable<FieldT> theta; // y_s - result_s
    pb_variable<FieldT> theta_;
    pb_variable<FieldT> delta; // y_s - x_s
    pb_variable<FieldT> delta_;
    pb_variable<FieldT> lambda; // lambda is the exp of theta, lambda = BASE**theta
    pb_variable<FieldT> a; // a = BASE**delta
    pb_variable<FieldT> b; // b = x_v * a + y_v
    pb_variable<FieldT> c; // c = b - result_v * lambda
    pb_variable<FieldT> d; // d = sigma * c
    pb_variable<FieldT> e; // e = lambda + 2^n - bound
    pb_variable<FieldT> f; // f is the most-significant bit of e
    pb_variable<FieldT> *e_decomposed;
    DecompositionCheckGadget<FieldT> *e_decomposition;
    pb_variable<FieldT> *result_v_decomposed;
    DecompositionCheckGadget<FieldT> *result_v_decomposition;

    u64 exp_base_omega = pow(__BASE, __omega) - 1; // use 1 to avoid overflow
    u64 sigma = pow(__BASE, __accuracy); // actually it should be sigma^-1

    pb_variable<FieldT> lessthan_result;
    LessThanGadget<FieldT> *lessthan_checker;
    ExponentialGadget<FieldT> *exponential_checker;
    PairwisePermutationGadget<FieldT> *pairwise_permutation_checker;

    const std::vector<pb_variable<FieldT>> &multiplies;
    const pb_variable<FieldT> &rand_point, &rand_coef;

    FloatAddGadget() {

    }

    FloatAddGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x_v, pb_variable <FieldT> &x_s, pb_variable <FieldT> &y_v, pb_variable <FieldT> &y_s, pb_variable <FieldT> &result_v, pb_variable <FieldT> &result_s, std::vector<pb_variable<FieldT>> &multiplies, const pb_variable <FieldT> &rand_point, const pb_variable <FieldT> &rand_coef, const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation), x_v(x_v), x_s(x_s), y_v(y_v), y_s(y_s), result_v(result_v), result_s(result_s), multiplies(multiplies), rand_point(rand_point), rand_coef(rand_coef) {
        // init pb val
        _init_pb_array(pb, result_v_decomposed, __omega, annotation + std::string("result_v_decomposed"));
        _init_pb_array(pb, e_decomposed, __max_delta_bit + 1, annotation + std::string("e_decomposed"));
        _init_pb_array(pb, a_1, 2, annotation + std::string("a_1"));
        _init_pb_array(pb, a_2, 2, annotation + std::string("a_2"));
        _init_pb_array(pb, b_1, 2, annotation + std::string("b_1"));
        _init_pb_array(pb, b_2, 2, annotation + std::string("b_2"));

        theta.allocate(pb, annotation + std::string("theta"));
        theta_.allocate(pb, annotation + std::string("theta_"));
        lambda.allocate(pb, annotation + std::string("lambda"));
        delta.allocate(pb, annotation + std::string("delta"));
        delta_.allocate(pb, annotation + std::string("delta_"));
        a.allocate(pb, annotation + std::string("a"));
        b.allocate(pb, annotation + std::string("b"));
        c.allocate(pb, annotation + std::string("c"));
        d.allocate(pb, annotation + std::string("d"));
        e.allocate(pb, annotation + std::string("e"));
        f.allocate(pb, annotation + std::string("f"));
        lessthan_result.allocate(pb, annotation + std::string("lessthan_result"));

        result_v_decomposition = new DecompositionCheckGadget<FieldT>(pb, &result_v, result_v_decomposed, 1, __omega, annotation + std::string("result_v_decomposition"));
        lessthan_checker = new LessThanGadget<FieldT>(pb, d, b, lessthan_result, 2 * __omega, annotation + std::string("lessthan_checker"));
        exponential_checker = new ExponentialGadget<FieldT>(pb, delta_, a, multiplies, __scale_n_bit, annotation + std::string("exponential_checker"));
        e_decomposition = new DecompositionCheckGadget<FieldT>(pb, &e, e_decomposed, 1, __max_delta_bit + 1, annotation + std::string("e_decomposition"));
        pairwise_permutation_checker = new PairwisePermutationGadget<FieldT>(pb, a_1, a_2, b_1, b_2, rand_coef, rand_point, 2,  annotation + std::string("pairwise_permutation_checker"));
    }

    ~FloatAddGadget() {

    }

    void generate_r1cs_constraints() {
        // connector
        add_r1cs(a_1[0], 1, x_v);
        add_r1cs(a_2[0], 1, x_s);
        add_r1cs(a_1[1], 1, y_v);
        add_r1cs(a_2[1], 1, y_s);
        // pairwize permutation check
        // after permutation, y_s >= x_s
        pairwise_permutation_checker->generate_r1cs_constraints();
        // result_v belong to [2^(w-1), 2^w)
        result_v_decomposition->generate_r1cs_constraints();
        add_r1cs(result_v_decomposed[0], 1, 1); // most-significant bit is 1
        // delta = y_s - x_s
        add_r1cs(delta, 1, b_2[1] - b_2[0]);
        // calculate delta - bound
        add_r1cs(delta - FieldT(__accuracy_bound) + FieldT((u64)pow(2, __max_delta_bit)) - e, 1, 0);
        // decompose e
        e_decomposition->generate_r1cs_constraints();
        // if most-significant bit is 0, smaller than bound. else bigger or equal.
        add_r1cs(e_decomposed[0], 1, f);
        // f * (x_v - result_v) = 0
        add_r1cs(f, b_1[0] - result_v, 0);
        // f * (x_s - result_s) = 0
        add_r1cs(f, b_2[0] - result_s, 0);
        // theta = y_s - result_s
        add_r1cs(theta, 1, b_2[1] - result_s);
        // theta is delta / delta+1
        add_r1cs(theta - delta, theta - delta - 1, 0);
        // delta => delta'
        add_r1cs(1 - f, delta - delta_, 0);
        // theta => theta'
        add_r1cs(1 - f, theta - theta_, 0);
        // a = BASE**delta 
        exponential_checker->generate_r1cs_constraints();
        // lambda = a * (theta - delta + 1)
        add_r1cs(a, theta_ - delta_ + 1, lambda);
        // b = x_v * a + y_v
        add_r1cs(a, b_1[0], b - b_1[1]);
        // c = b - result_v * lambda
        add_r1cs(result_v, lambda, b - c);
        // d = sigma * c
        add_r1cs(sigma, c, d);
        // d < b
        lessthan_checker->generate_r1cs_constraints();
        add_r1cs(1 - f, lessthan_result, 0);
    }

    void generate_r1cs_witness(Float x, Float y, Float result) {
        // connect
        this->pb.val(a_1[0]) = x.value;
        this->pb.val(a_2[0]) = x.scale;
        this->pb.val(a_1[1]) = y.value;
        this->pb.val(a_2[1]) = y.scale;
        Float _x, _y;
        // permutation check
        if(x.real_value < y.real_value) {
            _x = y;
            _y = x;
        } else {
            _x = x;
            _y = y; 
        } 
        this->pb.val(b_1[0]) = _x.value;
        this->pb.val(b_2[0]) = _x.scale;
        this->pb.val(b_1[1]) = _y.value;
        this->pb.val(b_2[1]) = _y.scale;
        pairwise_permutation_checker->generate_r1cs_witness();
        // result_v decomposed
        for (int i = 0; i < __omega; ++i) {
            this->pb.val(result_v_decomposed[i]) = (result.value >> (__omega - i - 1)) & 1U;
        }
        result_v_decomposition->generate_r1cs_witness();
        // delta
        this->pb.val(delta) = _y.scale - _x.scale;
        // f
        u64 e_value = _y.scale - _x.scale + (u64)pow(2, __max_delta_bit) - __accuracy_bound;
        this->pb.val(e) = e_value;
        for (int i = 0; i < __max_delta_bit + 1; ++i) {
            this->pb.val(e_decomposed[i]) = (e_value >> (__max_delta_bit - i)) & 1U;
        }
        e_decomposition->generate_r1cs_witness();
        this->pb.val(f) =  (e_value >> __max_delta_bit) & 1U;
        // theta
        this->pb.val(theta) = _y.scale - result.scale;

        if((e_value >> __max_delta_bit) & 1U) {
            // bigger or equal
            this->pb.val(delta_) = 0;
            this->pb.val(theta_) = 0;
            // s
            this->pb.val(a) = 1;
            this->pb.val(lambda) = 1;
            exponential_checker->generate_r1cs_witness(0);
            u64 b_value = _x.value + _y.value;
            this->pb.val(b) = b_value;
            u64 c_value = b_value - result.value;
            this->pb.val(c) = c_value;
            u64 d_value = sigma * c_value;
            this->pb.val(d) = libff::bigint<__limbs>(d_value);

            // less than
            lessthan_checker->generate_r1cs_witness(d_value, b_value);
            this->pb.val(lessthan_result) = d_value < b_value ? 0 : 1;
        } else {
            // smaller
            this->pb.val(delta_) = _y.scale - _x.scale;
            this->pb.val(theta_) = _y.scale - result.scale;
            // a
            mpz_t a_value, base;
            mpz_init(a_value);
            mpz_init_set_ui(base, __BASE);
            mpz_pow_ui(a_value, base, _y.scale - _x.scale);
            this->pb.val(a) = libff::bigint<__limbs>(a_value);
            // lambda
            mpz_t lambda_value, tmp; 
            mpz_init(lambda_value);
            mpz_init_set_ui(tmp, _x.scale - result.scale + 1);
            mpz_mul(lambda_value, a_value, tmp);
            this->pb.val(lambda) = libff::bigint<__limbs>(lambda_value);
            exponential_checker->generate_r1cs_witness(_y.scale - _x.scale);
            // b
            mpz_t b_value, tmp2;
            mpz_init(b_value);
            mpz_init(tmp2);
            mpz_mul_ui(tmp2, a_value, _x.value);
            mpz_add_ui(b_value, tmp2, _y.value);
            this->pb.val(b) = libff::bigint<__limbs>(b_value);
            // c
            mpz_t c_value, tmp3;
            mpz_init(c_value);
            mpz_init(tmp3);
            mpz_mul_ui(tmp3, lambda_value, result.value);
            mpz_sub(c_value, b_value, tmp3);
            this->pb.val(c) = libff::bigint<__limbs>(c_value);
            // d
            mpz_t d_value;
            mpz_init(d_value);
            mpz_mul_ui(d_value, c_value, sigma);
            this->pb.val(d) = libff::bigint<__limbs>(d_value);

            // less than
            lessthan_checker->generate_r1cs_witness(d_value, b_value);
            this->pb.val(lessthan_result) = 0;
        }
    }

};

#endif