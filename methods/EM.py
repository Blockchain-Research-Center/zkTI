import math
import csv
import random
import sys

class EM:
    def __init__(self,e2wl,worker_to_example_label,label_set, initquality=0.7):
        self.e2wl = e2wl
        self.worker_to_example_label = worker_to_example_label
        self.workers = self.worker_to_example_label.keys()
        self.label_set = label_set
        self.initalquality = initquality
             
    # E-step
    def Update_e2lpd(self):
        self.example_to_lable_probility = {}

        for example, worker_label_set in e2wl.items():
            lpd = {}
            total_weight = 0

            for tlabel, prob in self.label_probility.items():
                weight = prob
                for (worker, label) in worker_label_set:
                    weight *= self.w2cm[worker][tlabel][label]
                
                lpd[tlabel] = weight
                total_weight += weight

            for tlabel in lpd:
                if total_weight == 0:
                    # uniform distribution 
                    lpd[tlabel] = 1.0/len(self.label_set)
                else:
                    lpd[tlabel] = lpd[tlabel]*1.0/total_weight

            self.example_to_lable_probility[example] = lpd

    # M-step
    def Update_label_probility(self):
        for label in self.label_probility:
            self.label_probility[label] = 0

        for _, lpd in self.example_to_lable_probility.items():
            for label in lpd:
                self.label_probility[label] += lpd[label]

        for label in self.label_probility:
            self.label_probility[label] *= 1.0/len(self.example_to_lable_probility)


            
    def Update_w2cm(self):

        for worker in self.workers:
            for tlabel in self.label_set:
                for label in self.label_set:
                    self.w2cm[worker][tlabel][label] = 0

        worker_to_label_weight = {}
        for worker in self.worker_to_example_label:
            worker_to_label_weight[worker] = {}
            for label in self.label_set:
                worker_to_label_weight[worker][label] = 0
            for example, _ in self.worker_to_example_label[worker]:
                for label in self.label_set:
                    worker_to_label_weight[worker][label] += self.example_to_lable_probility[example][label]

            
            for tlabel in self.label_set:

                if worker_to_label_weight[worker][tlabel] == 0:
                    for label in self.label_set:
                        if tlabel == label:
                            self.w2cm[worker][tlabel][label] = self.initalquality
                        else:
                            self.w2cm[worker][tlabel][label] = (1-self.initalquality)*1.0/(len(self.label_set)-1)
                    continue

                for example, label in self.worker_to_example_label[worker]:
                    self.w2cm[worker][tlabel][label] += self.example_to_lable_probility[example][tlabel]*1.0/worker_to_label_weight[worker][tlabel]

        return self.w2cm

    #initialization
    def Init_l2pd(self):
        #uniform probability distribution
        label_probility = {}
        for label in self.label_set:
            label_probility[label] = 1.0/len(self.label_set)
        return label_probility

    def Init_w2cm(self):
        w2cm = {}
        for worker in self.workers:
            w2cm[worker] = {}
            for tlabel in self.label_set:
                w2cm[worker][tlabel] = {}
                for label in self.label_set:
                    if tlabel == label:
                        w2cm[worker][tlabel][label] = self.initalquality
                    else:
                        w2cm[worker][tlabel][label] = (1-self.initalquality)/(len(label_set)-1)

        return w2cm

    def Run(self, iterr = 20):
        
        self.label_probility = self.Init_l2pd()
        self.w2cm = self.Init_w2cm()

        while iterr > 0:
            # E-step
            self.Update_e2lpd() 

            # M-step
            # self.Update_label_probility()
            self.Update_w2cm()

            # compute the likelihood
            print(self.computelikelihood())

            iterr -= 1
        
        return self.example_to_lable_probility, self.w2cm


    def computelikelihood(self):
        
        lh = 0

        for _, worker_label_set in self.e2wl.items():
            temp = 0
            for tlabel, prior in self.label_probility.items():
                inner = prior
                for worker, label in worker_label_set:
                    inner *= self.w2cm[worker][tlabel][label]
                temp += inner
            
            lh += math.log(temp)
        
        return lh


###################################
# The above is the EM method (a class)
# The following are several external functions 
###################################

def getaccuracy(truthfile, example_to_lable_probility, label_set):
    e2truth = {}
    f = open(truthfile, 'r')
    reader = csv.reader(f)
    next(reader)

    for line in reader:
        example, truth = line
        e2truth[example] = truth

    tcount = 0
    count = 0

    for e in example_to_lable_probility:

        if e not in e2truth:
            continue

        temp = 0
        for label in example_to_lable_probility[e]:
            if temp < example_to_lable_probility[e][label]:
                temp = example_to_lable_probility[e][label]
        
        candidate = []

        for label in example_to_lable_probility[e]:
            if temp == example_to_lable_probility[e][label]:
                candidate.append(label)

        truth = random.choice(candidate)

        count += 1

        if truth == e2truth[e]:
            tcount += 1

    return tcount*1.0/count


def gete2wlandw2el(datafile):
    e2wl = {}
    worker_to_example_label = {}
    label_set=[]
    
    f = open(datafile, 'r')
    reader = csv.reader(f)
    next(reader)

    for line in reader:
        example, worker, label = line
        if example not in e2wl:
            e2wl[example] = []
        e2wl[example].append([worker,label])

        if worker not in worker_to_example_label:
            worker_to_example_label[worker] = []
        worker_to_example_label[worker].append([example,label])

        if label not in label_set:
            label_set.append(label)

    return e2wl,worker_to_example_label,label_set
    
def get_worker_quality(example_to_lable_probility={},  worker_to_example_label={}):
    worker_quality = {}
    total_count = len(example_to_lable_probility)
    
    for worker, el in worker_to_example_label.items():
        count = 0
        
        for e, w_label in el:
            if e not in example_to_lable_probility:
                continue
            else:
                temp = 0
                
                for label in example_to_lable_probility[e]:
                    if temp < example_to_lable_probility[e][label]:
                        temp = example_to_lable_probility[e][label]

                candidate = []

                for label in example_to_lable_probility[e]:
                    if temp == example_to_lable_probility[e][label]:
                        candidate.append(label)
                        
                if w_label in candidate:
                    count += 1
        
        worker_quality[worker] = count / total_count
                
    return worker_quality

if __name__ == "__main__":
    datafile = sys.argv[1]
    # generate structures to pass into EM
    e2wl,worker_to_example_label,label_set = gete2wlandw2el(datafile) 
    # EM iteration number
    iterations = 20 
    example_to_lable_probility, w2cm = EM(e2wl,worker_to_example_label,label_set).Run(iterations)

    print("w2cm: ", w2cm)
    print("......")
    print("example_to_lable_probility: ", example_to_lable_probility)
    print("......")

    truthfile = sys.argv[2]
    accuracy = getaccuracy(truthfile, example_to_lable_probility, label_set)
    print("accuracy: ", accuracy)
    print("......")
    
    worker_quality = get_worker_quality(example_to_lable_probility, worker_to_example_label)
    print("worker quality: ", worker_quality)
    

