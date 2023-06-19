#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <boost/optional.hpp>
#include <numeric>

#include "methods/MV.h"
#include "methods/ZenCrowd.h"
#include "circuits/MV_circuit.h"
#include "circuits/ZenCrowd_circuit.h"
#include "common.h"
#include "zk_proof_system/groth16.h"
#include "libsnark_exporter.h"

#include <libsnark/common/default_types/r1cs_gg_ppzksnark_pp.hpp>
#include <libsnark/relations/constraint_satisfaction_problems/r1cs/examples/r1cs_examples.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>

#define get_label_enum(x) (x == "1" ? 1 : x == "2" ? 2 \
                                      : x == "3"   ? 3 \
                                      : x == "4"   ? 4 \
                                                   : 0) // only consider the label in 5 enums

using namespace libsnark;

std::vector<std::vector<unsigned>> read_answer_data(const std::vector<std::vector<std::string>> &raw_answer_data)
{
    std::vector<std::vector<unsigned>> answer_data; // The transposed data to store the result
    std::map<std::string, int> question_map;        // Map question ids to indices
    std::map<std::string, int> worker_map;          // Map worker ids to indices

    // Create maps to map question ids and worker ids to indices
    for (const auto &row : raw_answer_data)
    {
        const std::string &question_id = row[0];
        const std::string &worker_id = row[1];

        if (question_map.count(question_id) == 0)
        {
            question_map[question_id] = question_map.size();
        }

        if (worker_map.count(worker_id) == 0)
        {
            worker_map[worker_id] = worker_map.size();
        }
    }

    // Resize the transposed data vector
    int num_questions = question_map.size();
    int num_workers = worker_map.size();
    answer_data.resize(num_questions, std::vector<unsigned>(num_workers));

    // Fill all elements as unlabeled
    for (int i = 0; i < answer_data.size(); i++)
    {
        for (int j = 0; j < answer_data[i].size(); j++)
        {
            answer_data[i][j] = EMPTY_LABEL;
        }
    }

    // Fill in the transposed data vector
    for (const auto &row : raw_answer_data)
    {
        const std::string &question_id = row[0];
        const std::string &worker_id = row[1];
        const std::string &label = row[2];

        int question_index = question_map[question_id];
        int worker_index = worker_map[worker_id];

        unsigned value = get_label_enum(label);

        answer_data[question_index][worker_index] = value;
    }

    // verify transposed data
    // for (int i = 0; i < transposed_data.size(); i++) {
    //     for(int j = 0; j < transposed_data.at(i).size(); j++) {
    //         std::cerr << transposed_data[i][j];
    //     }
    //     std::cerr << std::endl;
    // }

    return answer_data;
}

std::vector<unsigned> read_truth_data(std::vector<std::vector<std::string>> raw_answer_data, std::vector<std::vector<std::string>> raw_truth_data)
{
    std::vector<unsigned> truth_data;
    std::map<std::string, int> question_map; // Map question ids to indices
    std::map<std::string, int> worker_map;   // Map worker ids to indices

    // Create maps to map question ids and worker ids to indices
    for (const auto &row : raw_answer_data)
    {
        const std::string &question_id = row[0];
        const std::string &worker_id = row[1];

        if (question_map.count(question_id) == 0)
        {
            question_map[question_id] = question_map.size();
        }

        if (worker_map.count(worker_id) == 0)
        {
            worker_map[worker_id] = worker_map.size();
        }
    }

    // Resize the transposed data vector
    int num_questions = question_map.size();
    int num_workers = worker_map.size();
    truth_data.resize(num_questions);

    for (const auto &row : raw_truth_data)
    {
        const std::string &question_id = row[0];
        const std::string &truth = row[1];

        int question_index = question_map[question_id];

        unsigned value = get_label_enum(truth);

        truth_data[question_index] = value;
    }

    // verify transposed data
    // for (int i = 0; i < truth_data.size(); i++) {
    //     std::cerr << truth_data[i] << std::endl;
    // }

    return truth_data;
}

std::vector<std::vector<std::string>> read_csv_dataset(const std::string &filename)
{
    // The two-dimensional vector to store the data
    std::vector<std::vector<std::string>> read_data;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Error: Failed to open file.");
    }

    std::string line;
    bool jump_header = false;
    while (std::getline(file, line))
    {
        std::istringstream ss(line);
        std::vector<std::string> row;
        std::string value;

        // jump for table header
        if (jump_header == false)
        {
            jump_header = true;
            continue;
        }

        while (std::getline(ss, value, ','))
        {
            value.erase(std::remove(value.begin(), value.end(), '\r'), value.end());
            row.push_back(value);
        }

        read_data.push_back(row);
    }

    file.close();

    return read_data;
}

template <typename ppT>
void algo_MV(std::vector<std::vector<unsigned>> &answer_data, std::vector<unsigned> &truth_data)
{
    typedef libff::Fr<ppT> FieldT;
    default_r1cs_gg_ppzksnark_pp::init_public_params();

    std::cerr << "Task number: " << answer_data.size() << std::endl;
    std::cerr << "Worker number: " << answer_data[0].size() << std::endl;
    std::cerr << "Label number: " << answer_data.size() * answer_data[0].size() << std::endl;

    MV mv = MV(answer_data, truth_data);
    string circuit_name = "MV";
    std::cerr << "Run the Truth Inference algorithm: " << std::endl;
    std::vector<unsigned> result = mv.run();
    std::cerr << "Accuracy: " << mv.get_accuracy(result) << std::endl;

    protoboard<FieldT> pb;
    MajorityVoteCircuit<FieldT> majorityVoteCircuit = MajorityVoteCircuit<FieldT>(pb, mv, result, circuit_name);
    majorityVoteCircuit.generate_r1cs_constraints();
    majorityVoteCircuit.generate_r1cs_witness();

    std::cerr << "Constraints number: " << pb.num_constraints() << std::endl;
    std::cerr << "Public input number: " << pb.num_inputs() << std::endl;
    std::cerr << "Witness number: " << pb.auxiliary_input().size() << std::endl;
    std::cerr << "Variable number: " << pb.num_variables() << std::endl;
    std::cerr << "Protoboard satisfied: " << pb.is_satisfied() << std::endl;

    // export .zkif file
    std::cerr << "Start exporting circuit into .zkif file: "<< std::endl;
    zkifExporter<FieldT> exporter = zkifExporter<FieldT>(circuit_name, pb);
    exporter.export_protoboard();

    // Groth16 zk-SNARK
    run_r1cs_gg_ppzksnark<ppT>(pb, circuit_name + "_proof");
}

template <typename ppT>
void algo_ZC(std::vector<std::vector<unsigned>> &answer_data, std::vector<unsigned> &truth_data)
{
    typedef libff::Fr<ppT> FieldT;
    default_r1cs_gg_ppzksnark_pp::init_public_params();

    std::cerr << "Task number: " << answer_data.size() << std::endl;
    std::cerr << "Worker number: " << answer_data[0].size() << std::endl;
    std::cerr << "Label number: " << answer_data.size() * answer_data[0].size() << std::endl;

    ZenCrowd zc = ZenCrowd(answer_data, truth_data);
    string circuit_name = "ZC";
    std::cerr << "Run the Truth Inference algorithm: " << std::endl;
    std::vector<unsigned> result = zc.run();
    std::cerr << "Accuracy: " << zc.get_accuracy(result) << std::endl;

    protoboard<FieldT> pb;
    ZenCrowdCircuit<FieldT> zcCircuit = ZenCrowdCircuit<FieldT>(pb, zc, result, circuit_name);
    zcCircuit.generate_r1cs_constraints();
    zcCircuit.generate_r1cs_witness();

    std::cerr << "Constraints number: " << pb.num_constraints() << std::endl;
    std::cerr << "Public input number: " << pb.num_inputs() << std::endl;
    std::cerr << "Witness number: " << pb.auxiliary_input().size() << std::endl;
    std::cerr << "Varia ble number: " << pb.num_variables() << std::endl;
    std::cerr << "Protoboard satisfied: " << pb.is_satisfied() << std::endl;

    // export .zkif file
    std::cerr << "Start exporting circuit into .zkif file: "<< std::endl;
    zkifExporter<FieldT> exporter = zkifExporter<FieldT>(circuit_name, pb);
    exporter.export_protoboard();

    // Groth16 zk-SNARK
    run_r1cs_gg_ppzksnark<ppT>(pb, circuit_name + "_proof");
}

int main(int argc, char **argv)
{
    const std::string answer_file = argv[1];
    const std::string truth_file = argv[2];

    std::vector<std::vector<std::string>> raw_answer_data = read_csv_dataset(answer_file);
    std::vector<std::vector<std::string>> raw_truth_data = read_csv_dataset(truth_file);
    std::vector<std::vector<unsigned>> answer_data = read_answer_data(raw_answer_data);
    std::vector<unsigned> truth_data = read_truth_data(raw_answer_data, raw_truth_data);

    // std::cerr << "Run Groth16 zk-SNARK for zkTI MV algorithm: " << std::endl;
    // algo_MV<default_r1cs_gg_ppzksnark_pp>(answer_data, truth_data);
    // std::cerr << "Finish." << std::endl;

    std::cerr << "Run Groth16 zk-SNARK for zkTI ZC algorithm: " << std::endl;
    algo_ZC<default_r1cs_gg_ppzksnark_pp>(answer_data, truth_data);
    std::cerr << "Finish." << std::endl;
}