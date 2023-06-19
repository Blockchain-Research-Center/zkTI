#ifndef __zkTI_libsnark_exporter_h
#define __zkTI_libsnark_exporter_h

#include <vector>
#include <fstream>
#include <iostream>

#include "zkinterface_utils.hpp"
#include "zkinterface_generated.h"

#include "libff/common/default_types/ec_pp.hpp"
#include "libff/algebra/curves/alt_bn128/alt_bn128_init.hpp"
#include "libsnark/gadgetlib1/gadget.hpp"

using namespace std;
using namespace flatbuffers;
using namespace zkinterface_utils;
using namespace libsnark;
using libff::alt_bn128_r_limbs;
using libff::bigint;
using libff::bit_vector;
using libff::alt_bn128_q_bitcount;

const size_t fieldt_size = alt_bn128_q_bitcount;
const mp_size_t r_limbs = alt_bn128_r_limbs;

template <class FieldT>
class VarIdConverter
{
public:
  const flatbuffers::Vector<uint64_t> *input_ids;
  size_t input_count;
  uint64_t first_local_id;

  VarIdConverter(const CircuitHeader *circuit)
  {
    input_ids = circuit->instance_variables()->variable_ids();
    input_count = input_ids->size();
    first_local_id = circuit->free_variable_id();
  }

  uint64_t get_variable_id(const var_index_t &var_index);

  uint64_t get_local_id(size_t local_index);

  pb_variable<FieldT> get_local_variable(size_t local_index);

  uint64_t free_id_after_protoboard(const protoboard<FieldT> &pb);
};

template <typename FieldT>
uint64_t VarIdConverter<FieldT>::get_variable_id(const var_index_t &var_index)
{
  // Constant one?
  if (var_index == 0)
    return 0;

  // An input?
  size_t input_index = var_index - 1;
  if (input_index < input_count)
    return input_ids->Get(input_index);

  // A local variable.
  size_t local_index = input_index - input_count;
  return first_local_id + local_index;
}

template <typename FieldT>
uint64_t VarIdConverter<FieldT>::get_local_id(size_t local_index)
{
  return first_local_id + local_index;
}

template <typename FieldT>
pb_variable<FieldT> VarIdConverter<FieldT>::get_local_variable(size_t local_index)
{
  return 1 + input_count + local_index;
}

template <typename FieldT>
uint64_t VarIdConverter<FieldT>::free_id_after_protoboard(const protoboard<FieldT> &pb)
{
  size_t new_variables = pb.num_variables() - input_count;
  return first_local_id + new_variables;
}

template <class FieldT>
class zkifExporter
{
public:
  string circuit_name;
  protoboard<FieldT> pb;
  const CircuitHeader *ch;

  const string circuit_header_name = "header.zkif";
  const string constraints_name = "constraints.zkif";
  const string witness_name = "witness.zkif";

  zkifExporter(const string &circuit_name, const protoboard<FieldT> &pb)
  {
    this->circuit_name = circuit_name;
    this->pb = pb;
  }

  void export_protoboard()
  {
    std::cout << "Serializing circuit header..." << std::endl;
    auto ch_builder = serialize_circuit_header();
    std::cout << "Export circuit header into file " << circuit_name + "." + circuit_header_name << std::endl;
    std::cout << "Size: " << ch_builder.GetSize() << std::endl;
    export_builder_to_file(circuit_name + "." + circuit_header_name, (char*)ch_builder.GetBufferPointer(), ch_builder.GetSize());

    const CircuitHeader *ch = read_circuit_header((char *)ch_builder.GetBufferPointer());
    this->ch = ch;

    std::cout << "Serializing constraints..." << std::endl;
    auto cstr_builder = serialize_protoboard_constraints();
    std::cout << "Export constraints into file " << circuit_name + "." + constraints_name << std::endl;
    std::cout << "Size: " << cstr_builder.GetSize() << std::endl;
    export_builder_to_file(circuit_name + "." + constraints_name, (char*)cstr_builder.GetBufferPointer(), cstr_builder.GetSize());

    std::cout << "Serializing witnesses..." << std::endl;
    auto wtns_builder = serialize_protoboard_local_assignment();
    std::cout << "Export witnesses into file " << circuit_name + "." + witness_name << std::endl;
    std::cout << "Size: " << wtns_builder.GetSize() << std::endl;
    export_builder_to_file(circuit_name + "." + witness_name, (char*)wtns_builder.GetBufferPointer(), wtns_builder.GetSize());
  }

private:
  void export_builder_to_file(string file_name, char* buffer, size_t size) 
  {
    std::ofstream out(file_name, ios::binary);
    out.write(buffer, size);
    out.close();
  }

  Offset<Vector<Offset<KeyValue>>> make_configuration(FlatBufferBuilder &builder, vector<pair<string, string>> keyvalues)
  {
    vector<Offset<KeyValue>> config;
    // add config key-value pairs
    for (auto kv = keyvalues.begin(); kv != keyvalues.end(); kv++)
    {
      config.emplace_back(CreateKeyValue(builder, builder.CreateString(kv->first),
                                        0, builder.CreateString(kv->second)));
    }
    return builder.CreateVector(config);
  }

  FlatBufferBuilder serialize_circuit_header()
  {
    flatbuffers::FlatBufferBuilder builder;

    // instance_variables
    r1cs_primary_input<FieldT> public_inputs = pb.primary_input();
    vector<uint64_t> variable_ids(public_inputs.size());
    vector<uint8_t> values(fieldt_size * public_inputs.size());

    for (size_t i = 0; i < public_inputs.size(); i++)
    {
      variable_ids[i] = (uint64_t)(i + 1); // 1 for constant variable one
      FieldT value = public_inputs.at(i);
      into_le(value.as_bigint(), values.data() + fieldt_size * i, fieldt_size);
    }

    auto instance_variables = CreateVariables(builder, builder.CreateVector(variable_ids), builder.CreateVector(values));

    // free_variable_id
    auto free_variable_id = pb.num_inputs() + 1;

    // field_maximum
    const char *fieldt_maximum = "4294967295";
    auto fieldt_maximum_value = bigint<4L>(fieldt_maximum);
    vector<uint8_t> field_maximum_variable(fieldt_size);
    into_le(fieldt_maximum_value, field_maximum_variable.data(), fieldt_size);
    auto field_maximum = builder.CreateVector(field_maximum_variable);

    // configuration
    auto config = make_configuration(builder, {{"circuit_name", circuit_name}});

    // circuit header
    auto circuit_header = CreateCircuitHeader(builder, instance_variables, free_variable_id, field_maximum, config);

    auto root = CreateRoot(builder, Message_CircuitHeader, circuit_header.Union());

    builder.Finish(root);
    return builder;
  }

  FlatBufferBuilder serialize_protoboard_constraints()
  {
    VarIdConverter<FieldT> id_converter(ch);
    FlatBufferBuilder builder;

    // Closure: add a row of a matrix
    auto make_lc = [&](const vector<libsnark::linear_term<FieldT>> &terms)
    {
      vector<uint64_t> variable_ids(terms.size());
      vector<uint8_t> coeffs(fieldt_size * terms.size());

      for (size_t i = 0; i < terms.size(); i++)
      {
        variable_ids[i] = id_converter.get_variable_id(terms[i].index);
        into_le(terms[i].coeff.as_bigint(), coeffs.data() + fieldt_size * i, fieldt_size);
      }

      return CreateVariables(
          builder,
          builder.CreateVector(variable_ids),
          builder.CreateVector(coeffs));
    };

    // Send all rows of all three matrices
    auto lib_constraints = pb.get_constraint_system().constraints;

    // bilinear constraints
    vector<flatbuffers::Offset<BilinearConstraint>> fb_constraints;

    for (auto lib_constraint = lib_constraints.begin();
         lib_constraint != lib_constraints.end(); lib_constraint++)
    {
      fb_constraints.push_back(CreateBilinearConstraint(
          builder,
          make_lc(lib_constraint->a.terms),
          make_lc(lib_constraint->b.terms),
          make_lc(lib_constraint->c.terms)));
    }

    // constraints system
    auto constraint_system = CreateConstraintSystem(builder, builder.CreateVector(fb_constraints));

    auto root = CreateRoot(builder, Message_ConstraintSystem, constraint_system.Union());

    builder.Finish(root);
    return builder;
  }

  FlatBufferBuilder serialize_protoboard_local_assignment()
  {
    VarIdConverter<FieldT> id_converter(ch);
    FlatBufferBuilder builder;

    size_t input_count = id_converter.input_count;
    size_t new_count = pb.num_variables() - input_count;

    vector<uint64_t> variable_ids(new_count);
    vector<uint8_t> elements(fieldt_size * new_count);

    for (size_t i = 0; i < new_count; ++i)
    {
      variable_ids[i] = id_converter.get_local_id(i);
      auto pb_var = id_converter.get_local_variable(i);
      into_le(pb.val(pb_var).as_bigint(), elements.data() + fieldt_size * i, fieldt_size);
    }

    auto values = CreateVariables(
        builder,
        builder.CreateVector(variable_ids),
        builder.CreateVector(elements));

    // witnesses
    auto witness = CreateWitness(builder, values);

    auto root = CreateRoot(builder, Message_Witness, witness.Union());

    builder.Finish(root);
    return builder;
  }
};

#endif