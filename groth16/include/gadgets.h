#ifndef __zkTI_gadgets_h
#define __zkTI_gadgets_h

#include <common.h>

template <typename FieldT>
class ComparisonGadgets: public gadget<FieldT> {

public:
    const pb_variable<FieldT> &comparison_result;
    const pb_variable<FieldT> &x, &y, 
    const pb_variable<FieldT> &diff;

    ComparisonGadget(protoboard <FieldT> &pb, pb_variable <FieldT> &x_, pb_variable <FieldT> &y_,
                     pb_variable <FieldT> &comparison_result_, pb_variable <FieldT> &diff_,
                     const std::string &annotation = "")
            : gadget<FieldT>(pb, annotation), x(x_), y(y_), comparison_result(comparison_result_), diff(diff_) {
    }

    ~ComparisonGadget() {}

    void generate_r1cs_constraints() {
        add_r1cs(2 * comparison_result - 1, y - x, diff);
    }

    void generate_r1cs_witness() {}

};

#endif