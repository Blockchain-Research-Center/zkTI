// zkInterface - integration helpers.
//
// @author Aurélien Nicolas <info@nau.re> for QED-it.com
// @date 2020

#include "zkinterface_utils.hpp"

namespace zkinterface_utils {

    const uoffset_t UNKNOWN_BUFFER_SIZE = uoffset_t(4) * 1000 * 1000 * 1000; // 4G.

    void into_le(const libff::bigint<r_limbs> &num, uint8_t *out, size_t size) {
        size_t bytes_per_limb = sizeof(num.data[0]);
        assert(size >= bytes_per_limb * r_limbs);

        for (size_t byte = 0; byte < size; byte++)
        {
        size_t limb = byte / bytes_per_limb;
        size_t limb_byte = byte % bytes_per_limb;
        out[byte] = uint8_t(num.data[limb] >> (limb_byte * 8));
        }
    }


    uoffset_t read_size_prefix(void *buffer) {
        uoffset_t message_length = ReadScalar<uoffset_t>(buffer);
        return sizeof(uoffset_t) + message_length;
    }

    const CircuitHeader* read_circuit_header(char *buffer) {
        return GetRoot(buffer)->message_as_CircuitHeader();
    }

    const ConstraintSystem* read_constraint_system(char *buffer) {
        return GetRoot(buffer)->message_as_ConstraintSystem();
    }

    const Witness* read_witness(char *buffer) {
        return GetRoot(buffer)->message_as_Witness();
    }
} 