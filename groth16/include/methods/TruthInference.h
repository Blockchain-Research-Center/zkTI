#ifndef __truth_inference_h
#define __truth_inference_h

#include <vector>

using namespace std;

class TruthInference {

protected:
    std::vector<std::vector<unsigned>> answer_data;
    std::vector<unsigned> truth_data;

public:
    TruthInference(std::vector<std::vector<unsigned>>& _answer_data, std::vector<unsigned> &_truth_data) {
        answer_data = _answer_data;
        truth_data = _truth_data;
    }

    ~TruthInference() {

    }

};

#endif