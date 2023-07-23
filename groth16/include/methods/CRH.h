#ifndef __crh_h
#define __crh_h

#include <cassert>
#include <iostream>
#include <cassert>
#include <vector>
#include <set>

#include "./Float.h"
#include "../common.h"
#include "./TruthInference.h"

class CRH : public TruthInference
{

public:
    std::vector<std::vector<int>> label_counts;
    int label_class_number;
    Float init_quality;

public:
    CRH(std::vector<std::vector<unsigned>> &_answer_data, std::vector<unsigned> &_truth_data) : TruthInference(_answer_data, _truth_data)
    {
        label_class_number = 0;
    }

    std::vector<unsigned> run(int iter_time = 1, Float init_quality = Float(0.8f))
    {
        assert(answer_data.size() == truth_data.size());
        int task_number = answer_data.size();
        int worker_number = answer_data[0].size();
        this->init_quality = init_quality;

        // count label class
        std::set<unsigned> labels;
        for (int i = 0; i < task_number; i++)
        {
            for (int j = 0; j < worker_number; j++)
            {
                if (labels.find(answer_data[i][j]) == labels.end()) {
                    labels.insert(answer_data[i][j]);
                }
            }
        }
        label_class_number = labels.size();

        // initialize worker quality 
        std::vector<Float> worker_quality = std::vector<Float>(worker_number);
        for(int i = 0; i < worker_number; i++) {
            worker_quality[i] = init_quality;
        }

        // below is the CRH
        std::vector<std::vector<Float>> task_label_prediction = std::vector<std::vector<Float>>(task_number);

        for (int iter = 0; iter < iter_time; iter++) {
            // E Step
            for (int i = 0; i < task_number; i++){
                task_label_prediction[i].resize(label_class_number);
                for(int j = 0; j < label_class_number; j++) {
                    task_label_prediction[i][j] = Float::zero();
                }

                for(int j = 0; j < worker_number; j++) {
                    unsigned label = answer_data[i][j];
                    task_label_prediction[i][label] = task_label_prediction[i][label] + worker_quality[j];
                }

                Float sums = Float::zero();
                for(int k = 0; k < label_class_number; k++) {
                    sums = sums + task_label_prediction[i][k];
                }

                for(int k = 0; k < label_class_number; k++) {
                    task_label_prediction[i][k] = task_label_prediction[i][k] / sums;
                }
            }

            // M Step
            Float weight_sum = Float::zero();
            for(int j = 0; j < worker_number; j++) {
                Float distance_sum = Float::zero();
                for(int i = 0; i < task_number; i++) {
                    unsigned label = answer_data[i][j];
                    // distance calculation
                    Float max_prob = Float::zero();
                    unsigned imax_prob = -1;
                    for(int k = 0; k < label_class_number; k++) {
                        if (task_label_prediction[i][k].real_value > max_prob.real_value) {
                            max_prob = task_label_prediction[i][k];
                            imax_prob = k;
                        }
                    }

                    if(label != imax_prob) {
                        distance_sum = distance_sum + Float::one();                       
                    }
                }
                worker_quality[j] = distance_sum;
                weight_sum = weight_sum + distance_sum;
            }

            for(int j = 0; j < worker_number; j++) {
                worker_quality[j] = weight_sum / worker_quality[j];
                // std::cout << worker_quality[j] << std::endl;
            }    
        }

        // calculate result 
        std::vector<unsigned> result = std::vector<unsigned>(task_number);
        for(int i = 0; i < task_number; i++) {
            float max = 0;
            int imax = 0;
            for(int k = 0; k < label_class_number; k++) {
                if(task_label_prediction[i][k].real_value >= max) {
                    max = task_label_prediction[i][k].real_value;
                    imax = k;
                }

            }
            result[i] = imax;
        }

        return result;
    }

    float get_accuracy(std::vector<unsigned> &predicted_data)
    {
        assert(truth_data.size() == predicted_data.size());
        int right = 0;
        for (int i = 0; i < truth_data.size(); i++)
        {
            if (truth_data[i] == predicted_data[i])
            {
                right++;
            }
        }
        return ((float)right) / truth_data.size();
    }
};

#endif