#ifndef __zencrowd_h
#define __zencrowd_h

#include <cassert>
#include <iostream>
#include <cassert>
#include <vector>
#include <set>

#include "./methods/Float.h"
#include "../common.h"
#include "./TruthInference.h"

class ZenCrowd : public TruthInference
{

public:
    std::vector<std::vector<int>> label_counts;
    int label_class_number;

public:
    ZenCrowd(std::vector<std::vector<unsigned>> &_answer_data, std::vector<unsigned> &_truth_data) : TruthInference(_answer_data, _truth_data)
    {
        label_class_number = 0;
    }

    std::vector<unsigned> run(int iter_time = 5, Float init_quanlity = Float(0.8f))
    {
        assert(answer_data.size() == truth_data.size());
        int task_number = answer_data.size();
        int worker_number = answer_data[0].size();

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
        int label_class_number = labels.size();

        // initialize worker quality 
        std::vector<Float> worker_quanlity = std::vector<Float>(worker_number);
        for(int i = 0; i < worker_number; i++) {
            worker_quanlity[i] = init_quanlity;
        }

        // below is the ZenCrowd
        std::vector<std::vector<Float>> task_label_prediction = std::vector<std::vector<Float>>(task_number);

        // Float min = Float(100.0f);
        // Float max = Float();

        for (int iter = 0; iter < iter_time; iter++) {
            // E Step
            for (int i = 0; i < task_number; i++){
                task_label_prediction[i].resize(label_class_number);
                for(int j = 0; j < label_class_number; j++) {
                    task_label_prediction[i][j] = Float::one();
                }

                for(int j = 0; j < worker_number; j++) {
                    unsigned label = answer_data[i][j];
                    for(int k = 0; k < label_class_number; k++) {
                        if(label == k) {
                            task_label_prediction[i][k] = task_label_prediction[i][k] * worker_quanlity[k];
                        } else {
                            task_label_prediction[i][k] = task_label_prediction[i][k] * (worker_quanlity[k].one_minus());
                        }

                        // if(task_label_prediction[i][k].real_value > max.real_value) {
                        //     std::cout << "max: " << task_label_prediction[i][k] << std::endl;
                        //     max = task_label_prediction[i][k];
                        // }

                        // if(task_label_prediction[i][k].real_value < min.real_value) {
                        //     std::cout << "min: " << task_label_prediction[i][k] << std::endl;
                        //     min = task_label_prediction[i][k];
                        // }
                    }
                }

                Float sums = Float::zero();
                for(int k = 0; k < label_class_number; k++) {
                    // std::cout << task_label_prediction[i][k] << " + " << sums << std::endl;
                    sums = sums + task_label_prediction[i][k];
                    // std::cout << sums << std::endl;
                }

                // if(sums.real_value > max.real_value) {
                //     std::cout << "max: " << sums << std::endl;
                //     max = sums;
                // }

                for(int k = 0; k < label_class_number; k++) {
                    // std::cout << task_label_prediction[i][k] << " // " << sums << std::endl;
                    task_label_prediction[i][k] = task_label_prediction[i][k] / sums;
                    // std::cout << "i: " << i << " k: " << k << " p: " << task_label_prediction[i][k] << std::endl;

                    // if(task_label_prediction[i][k].real_value > max.real_value) {
                    //     std::cout << "max: " << task_label_prediction[i][k] << std::endl;
                    //     max = task_label_prediction[i][k];
                    // }

                    // if(task_label_prediction[i][k].real_value < min.real_value) {
                    //     std::cout << "min: " << task_label_prediction[i][k] << std::endl;
                    //     min = task_label_prediction[i][k];
                    // }
                }
            }

            // M Step
            for(int j = 0; j < worker_number; j++) {
                worker_quanlity[j] = Float::zero();
                for(int i = 0; i < task_number; i++) {
                    unsigned label = answer_data[i][j];
                    // std::cout << worker_quanlity[j] << " + " << task_label_prediction[i][label] / Float(task_number) << std::endl;
                    // std::cout << task_label_prediction[i][label] << std::endl;
                    // std::cout << Float(task_number) << std::endl;
                    worker_quanlity[j] = worker_quanlity[j] + task_label_prediction[i][label] / Float(task_number);
                    // std::cout << worker_quanlity[j] << std::endl;

                    // if(worker_quanlity[j].real_value > max.real_value) {
                    //     std::cout << "max: " << worker_quanlity[j] << std::endl;
                    //     max = worker_quanlity[j];
                    // }

                    // if(worker_quanlity[j].real_value < min.real_value) {
                    //     std::cout << "min: " << worker_quanlity[j] << std::endl;
                    //     min = worker_quanlity[j];
                    // }
                }
                // std::cout << j << ": " << worker_quanlity[j] << std::endl;
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