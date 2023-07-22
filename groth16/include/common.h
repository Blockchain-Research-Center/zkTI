#ifndef _common_circuit_h
#define _common_circuit_h

#include <libsnark/gadgetlib1/gadget.hpp>
#include <gmp.h>

using namespace libsnark;

#define EMPTY_LABEL 100;
#define add_r1cs(x, y, z) this->pb.add_r1cs_constraint(r1cs_constraint<FieldT>(x, y, z))
#define pb_eval(x) this->pb.val(x)
#define __limbs 4L

template <typename FieldT>
void _init_pb_array(protoboard<FieldT> &pb, pb_variable<FieldT> *&array, int length, std::string &&name)
{
    array = new pb_variable<FieldT>[length];
    for (int i = 0; i < length; ++i)
    {
        array[i].allocate(pb, name + std::string("_") + std::to_string(i));
    }
}

template <typename FieldT>
void init_one_dimension_vec(protoboard<FieldT> &pb, std::vector<pb_variable<FieldT>> &pb_vec, int x_n, const std::string &prefix_name)
{
    for (int i = 0; i < x_n; i++)
    {
        pb_variable<FieldT> var;
        var.allocate(pb, prefix_name + "_" + std::to_string(i));
        pb_vec.push_back(var);
    }
}

template <typename FieldT>
void init_two_dimension_vec(protoboard<FieldT> &pb, std::vector<std::vector<pb_variable<FieldT>>> &pb_vec, int x_n, int y_n, const std::string &prefix_name)
{
    for (int i = 0; i < x_n; i++)
    {
        std::vector<pb_variable<FieldT>> row_x;
        for (int j = 0; j < y_n; j++)
        {
            pb_variable<FieldT> var;
            var.allocate(pb, prefix_name + "_" + std::to_string(i) + "_" + std::to_string(j));
            row_x.push_back(var);
        }
        pb_vec.push_back(row_x);
    }
}

template <typename FieldT>
void init_three_dimension_vec(protoboard<FieldT> &pb, std::vector<std::vector<std::vector<pb_variable<FieldT>>>> &pb_vec, int x_n, int y_n, int z_n, const std::string &prefix_name)
{
    for (int i = 0; i < x_n; ++i)
    {
        std::vector<std::vector<pb_variable<FieldT>>> row_x;
        for (int j = 0; j < y_n; ++j)
        {
            std::vector<pb_variable<FieldT>> row_y;
            for (int k = 0; k < z_n; ++k)
            {
                pb_variable<FieldT> var;
                var.allocate(pb, prefix_name + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k));
                row_y.push_back(var);
            }
            row_x.push_back(row_y);
        }
        pb_vec.push_back(row_x);
    }
}

template <typename FieldT>
void init_four_dimension_vec(protoboard<FieldT> &pb, std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>> &pb_vec, int x_n, int y_n, int z_n, int l_n, const std::string &prefix_name)
{
    for (int i = 0; i < x_n; ++i)
    {
        std::vector<std::vector<std::vector<pb_variable<FieldT>>>> row_x;
        for (int j = 0; j < y_n; ++j)
        {
            std::vector<std::vector<pb_variable<FieldT>>> row_y;
            for (int k = 0; k < z_n; ++k)
            {
                std::vector<pb_variable<FieldT>> row_k;
                for (int l = 0; l < l_n; ++l)
                {
                    pb_variable<FieldT> var;
                    var.allocate(pb, prefix_name + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" + std::to_string(l));
                    row_k.push_back(var);
                }
                row_y.push_back(row_k);
            }
            row_x.push_back(row_y);
        }
        pb_vec.push_back(row_x);
    }
}

template <typename FieldT>
void init_five_dimension_vec(protoboard<FieldT> &pb, std::vector<std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>>> &pb_vec, int x_n, int y_n, int z_n, int l_n, int m_n, const std::string &prefix_name)
{
    for (int i = 0; i < x_n; ++i)
    {
        std::vector<std::vector<std::vector<std::vector<pb_variable<FieldT>>>>> row_x;
        for (int j = 0; j < y_n; ++j)
        {
            std::vector<std::vector<std::vector<pb_variable<FieldT>>>> row_y;
            for (int k = 0; k < z_n; ++k)
            {
                std::vector<std::vector<pb_variable<FieldT>>> row_k;
                for (int l = 0; l < l_n; ++l)
                {
                    std::vector<pb_variable<FieldT>> row_l;
                    for (int m = 0; m < m_n; ++m)
                    {
                        pb_variable<FieldT> var;
                        var.allocate(pb, prefix_name + "_" + std::to_string(i) + "_" + std::to_string(j) + "_" + std::to_string(k) + "_" + std::to_string(l) + "_" + std::to_string(m));
                        row_l.push_back(var);
                    }
                    row_k.push_back(row_l);
                }
                row_y.push_back(row_k);
            }
            row_x.push_back(row_y);
        }
        pb_vec.push_back(row_x);
    }
}

#endif