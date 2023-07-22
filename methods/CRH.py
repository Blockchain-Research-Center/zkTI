import math
import csv
import random
import sys

class CRH:
    def __init__(self, e2wl, w2el, label_set, e2t):
        self.e2wl = e2wl
        self.w2el = w2el
        self.weight = {}
        self.label_set = label_set
        self.e2t = e2t

    def distance_calculation(self, example, label):
        if self.truth[example][label] >= (1.0 / len(self.label_set)):
            return 0.0
        else:
            return 1.0

    def examples_truth_calculation(self):
        self.truth = {}

        for example, worker_label_set in self.e2wl.items():
            self.truth[example] = {}
            
            for label in self.label_set:
                self.truth[example][label]= 0

            for worker, label in worker_label_set:
                self.truth[example][label] = self.truth[example][label] + self.weight[worker]
            
            sums = 0
            for label, weight in self.truth[example].items():
                sums += weight
            
            if sums==0:
                for label in self.truth[example].keys():
                    self.truth[example][label] = 1.0 / len(self.label_set)
            else:
                for label in self.truth[example].keys():
                    self.truth[example][label] = self.truth[example][label] * 1.0 / sums 

    def workers_weight_calculation(self):
        weight_sum = 0.0

        self.weight = {}

        for worker, example_label_set in self.w2el.items():
            distance_sum = 0.0
            for example, label in example_label_set:
                distance_sum = distance_sum + self.distance_calculation(example, label)
            if distance_sum == 0.0:
                distance_sum = 0.00000001

            self.weight[worker] = distance_sum
            weight_sum = weight_sum + distance_sum

        for worker in self.w2el.keys():
            print(weight_sum / self.weight[worker])
            self.weight[worker] = weight_sum / self.weight[worker]

    def init_truth(self):
        self.truth = {}
        self.std = {}

        # using majority voting to obtain initial value
        for example, worker_label_set in self.e2wl.items():
            self.truth[example] = {}
            temp = {}
            for _, label in worker_label_set:
                if (label in temp.keys()):
                    temp[label] = temp[label] + 1
                else:
                    temp[label] = 1

            max = 0
            for label, num in temp.items():
                if num > max:
                    max = num

            candidate = []
            for label, num in temp.items():
                if max == num:
                    candidate.append(label)

            init_truth_label = random.choice(candidate)
            for label in self.label_set:
                if label == init_truth_label:
                    self.truth[example][label] = 1.0
                else:
                    self.truth[example][label] = 0.0
                    
    def init_weight(self):
        for worker in self.w2el.keys():
            self.weight[worker] = 0.8
                
    def get_e2lpd(self):
        e2lpd = {}
        for example, worker_label_set in self.e2wl.items():
            temp = {}
            sum = 0.0
            for worker, label in worker_label_set:
                if (label in temp.keys()):
                    temp[label] = temp[label] + self.weight[worker]
                else:
                    temp[label] = self.weight[worker]
                sum = sum + self.weight[worker]

            for label in temp.keys():
                temp[label] = temp[label] / sum

            e2lpd[example] = temp

        return e2lpd

    def get_workerquality(self):
        sum_worker = sum(self.weight.values())
        norm_worker_weight = {}
        for worker in self.weight.keys():
            norm_worker_weight[worker] = self.weight[worker] / sum_worker
        return norm_worker_weight

    def Run(self, iter):
        self.init_truth()
        self.init_weight()
        
        while iter > 0:
            self.examples_truth_calculation()
            self.workers_weight_calculation()
            iter -= 1

        return self.get_e2lpd(), self.weight


###################################
# The above is the CRH method (a class)
# The following are several external functions
###################################

def getaccuracy(truthfile, predict_truth):
    e2truth = {}
    f = open(truthfile, 'r')
    reader = csv.reader(f)
    next(reader)

    for line in reader:
        example, truth = line
        e2truth[example] = truth

    tcount = 0
    count = 0

    for e, ptruth in predict_truth.items():

        if e not in e2truth:
            continue

        count += 1

        truth = '0'
        if ptruth['0'] < ptruth['1']:
            truth = '1'
        if truth == e2truth[e]:
            tcount += 1

    return tcount*1.0/count


def gete2wlandw2el(datafile):
    e2wl = {}
    w2el = {}
    label_set = []

    f = open(datafile, 'r')
    reader = csv.reader(f)
    next(reader)

    for line in reader:
        example, worker, label = line
        if example not in e2wl:
            e2wl[example] = []
        e2wl[example].append([worker, label])

        if worker not in w2el:
            w2el[worker] = []
        w2el[worker].append([example, label])

        if label not in label_set:
            label_set.append(label)

    return e2wl, w2el, label_set


def gete2t(known_true):
    e2t = {}
    f = open(known_true)
    reader = csv.reader(f)
    next(reader)

    for line in reader:
        example, truth = line
        e2t[example] = truth

    f.close()
    return e2t


if __name__ == "__main__":

    datafile = sys.argv[1]
    e2wl, w2el, label_set = gete2wlandw2el(datafile)

    e2t = gete2t(sys.argv[2])

    e2lpd, weight = CRH(e2wl, w2el, label_set, e2t).Run(5)

    # print(weight)
    # print(e2lpd)

    truthfile = sys.argv[2]
    accuracy = getaccuracy(truthfile, e2lpd)
    print("accuracy: ", accuracy)
