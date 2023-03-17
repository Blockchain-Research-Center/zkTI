import copy
import random
import sys
import csv


class MV:

    def __init__(self,e2wl,w2el,label_set):

        self.e2wl = e2wl
        self.w2el = w2el
        self.workers = self.w2el.keys()
        self.label_set = label_set


    def Run(self):

        e2wl = self.e2wl
        e2lpd={}
        for e in e2wl:
            e2lpd[e]={}

            # multi label
            for label in self.label_set:
                e2lpd[e][label] = 0

            for item in e2wl[e]:
                label=item[1]
                e2lpd[e][label]+= 1
                
            alls = 0
            
            for label in self.label_set:
                alls += e2lpd[e][label]
                
            if alls!=0:
                for label in self.label_set:
                    e2lpd[e][label] = 1.0 * e2lpd[e][label] / alls
            else:
                for label in self.label_set:
                    e2lpd[e][label] = 1.0 / len(self.label_set)

        return e2lpd


def gete2wlandw2el(datafile):
    e2wl = {}
    w2el = {}
    label_set=[]
    
    f = open(datafile, 'r')
    reader = csv.reader(f)
    next(reader)

    for line in reader:
        example, worker, label = line
        if example not in e2wl:
            e2wl[example] = []
        e2wl[example].append([worker,label])

        if worker not in w2el:
            w2el[worker] = []
        w2el[worker].append([example,label])

        if label not in label_set:
            label_set.append(label)

    return e2wl,w2el,label_set

def getaccuracy(truthfile, e2lpd, label_set):
    e2truth = {}
    f = open(truthfile, 'r')
    reader = csv.reader(f)
    next(reader)

    for line in reader:
        example, truth = line
        e2truth[example] = truth

    tcount = 0
    count = 0

    for e in e2lpd:

        if e not in e2truth:
            continue

        temp = 0
        for label in e2lpd[e]:
            if temp < e2lpd[e][label]:
                temp = e2lpd[e][label]

        candidate = []

        for label in e2lpd[e]:
            if temp == e2lpd[e][label]:
                candidate.append(label)

        truth = random.choice(candidate)

        count += 1

        if truth == e2truth[e]:
            tcount += 1

    return tcount*1.0/count

def get_worker_quality(e2lpd={},  w2el={}):
    worker_quality = {}
    total_count = len(e2lpd)
    
    for worker, el in w2el.items():
        count = 0
        
        for e, w_label in el:
            if e not in e2lpd:
                continue
            else:
                temp = 0
                
                for label in e2lpd[e]:
                    if temp < e2lpd[e][label]:
                        temp = e2lpd[e][label]

                candidate = []

                for label in e2lpd[e]:
                    if temp == e2lpd[e][label]:
                        candidate.append(label)
                        
                if w_label in candidate:
                    count += 1
        
        worker_quality[worker] = count / total_count
                
    return worker_quality
        
if __name__ == "__main__":
    datafile = sys.argv[1]
    e2wl,w2el,label_set = gete2wlandw2el(datafile)
    e2lpd = MV(e2wl,w2el,label_set).Run()
    print("e2lpd: ", e2lpd)
    print("......")

    truthfile = sys.argv[2]
    accuracy = getaccuracy(truthfile, e2lpd, label_set)
    print("accuracy: ", accuracy)
    print("......")
    
    worker_quality = get_worker_quality(e2lpd, w2el)
    print("worker quality: ", worker_quality)
    

