#ifndef ZKIF_ZKINTERFACE_UTILS_HPP
#define ZKIF_ZKINTERFACE_UTILS_HPP

#include "zkinterface_generated.h"

#include "libsnark/gadgetlib1/gadget.hpp"
#include "libff/algebra/curves/alt_bn128/alt_bn128_init.hpp"
#include "libff/common/default_types/ec_pp.hpp"

namespace zkinterface_utils {

    using namespace std;
    using namespace flatbuffers;
    using namespace zkinterface;
    using namespace libsnark;

    using libff::alt_bn128_r_limbs;
    using libff::bigint;

    const mp_size_t r_limbs = alt_bn128_r_limbs;

    // utils
    void into_le(const bigint<r_limbs> &num, uint8_t *out, size_t size);
    
    // reader
    uoffset_t read_size_prefix(void *buffer);

    const CircuitHeader* read_circuit_header(char *buffer);

    const ConstraintSystem* read_constraint_system(char *buffer);

    const Witness* read_witness(char *buffer);

    class MessageNotFoundException : public std::exception {
    public:
        inline const char *what() const throw() {
            return "message of the required type not found";
        }
    };

}
#endif //ZKIF_ZKINTERFACE_UTILS_HPP
