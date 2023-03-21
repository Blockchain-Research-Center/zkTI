#ifndef __majority_vote_h
#define __majority_vote_h

#include "TruthInference.h"
#include "../common.h"
#include <iostream>
#include <cassert>
#include <vector>
#include <map>


class MV : public TruthInference {

public:
    std::vector<std::vector<int>> label_counts;
    int label_class_number;

public:
    MV(std::vector<std::vector<unsigned>>& _answer_data, std::vector<unsigned> &_truth_data) : 
    TruthInference(_answer_data, _truth_data) {
        label_class_number = 0;
    }

    std::vector<unsigned> run() {
        assert(answer_data.size() == truth_data.size());

        std::vector<unsigned> result = std::vector<unsigned>(truth_data.size());
        for(int i = 0; i < answer_data.size(); i++) {
            std::map<unsigned, int> counts;
            for(int j = 0; j < answer_data[i].size(); j++) {
                if(counts.find(answer_data[i][j]) == counts.end()) {
                    counts.insert(pair<unsigned, int>(answer_data[i][j], 1));
                } else {
                    counts[answer_data[i][j]]++;
                }
            }

            unsigned major_label = EMPTY_LABEL;
            int num = 0;
            for(auto &p : counts) {
                if(p.second > num) {
                    major_label = p.first;
                    num = p.second;
                }
            }

            result[i] = major_label;

            // update total label class
            if(counts.size() > label_class_number) {
                label_class_number = counts.size();
            }

            // update label count
            std::vector<int> row;
            row.resize(label_class_number);
            for(auto &p : counts) {
                assert(p.first < row.size());
                row[p.first] = p.second;
            }
            label_counts.push_back(row);
        }

        return result;
    }

    float get_accuracy(std::vector<unsigned> &predicted_data) {
        assert(truth_data.size() == predicted_data.size());
        int right = 0;
        for(int i = 0; i < truth_data.size(); i++) {
            if(truth_data[i] == predicted_data[i]) {
                right++;
            }
        }
        return ((float)right) / truth_data.size(); 
    }
};

#endif