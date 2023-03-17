import copy
import sys
import random
import csv


class ZenClowd:

    def __init__(self, e2wl, w2el, label_set):

        self.e2wl = e2wl
        self.w2el = w2el
        self.label_set = label_set

    def InitPj(self):
        l2pd={}
        
        for label in self.label_set:
            l2pd[label]=1.0/len(self.label_set)
            
        return l2pd


    def Initw2cm(self, workers):
        w2cm={}

        if workers=={}:
            workers=self.w2el.keys()
            for worker in workers:
                w2cm[worker]=0.8
        else:
            for worker in workers:
                if worker not in w2cm: # workers --> w2cm
                    w2cm[worker] = 0.8
                else:
                    w2cm[worker]=workers[worker]

        return w2cm

    # E-step
    def ComputeTij(self, e2wl, l2pd, w2cm):
        e2lpd={}
        for e, workerlabels in e2wl.items():
            e2lpd[e]={}
            for label in self.label_set:
                e2lpd[e][label]=1.0  #l2pd[label]

            for worker,label in workerlabels:
                for candlabel in self.label_set:
                    if label==candlabel:
                        e2lpd[e][candlabel]*=w2cm[worker]
                    else:
                        e2lpd[e][candlabel]*=(1-w2cm[worker])*1.0/(len(self.label_set)-1)

            sums=0
            for label in self.label_set:
                sums+=e2lpd[e][label]

            if sums==0:
                for label in self.label_set:
                    e2lpd[e][label]=1.0/self.len(self.label_set)
            else:
                for label in self.label_set:
                    e2lpd[e][label]=e2lpd[e][label]*1.0/sums

        #print e2lpd
        return e2lpd


    # M-step
    def ComputePj(self, e2lpd):
        l2pd = {}

        for label in self.label_set:
            l2pd[label]=0
        for e in e2lpd:
            for label in e2lpd[e]:
                l2pd[label]+=e2lpd[e][label]

        for label in self.label_set:
            l2pd[label]=l2pd[label]*1.0/len(e2lpd)

        return l2pd


    def Computew2cm(self, w2el, e2lpd):
        w2cm={}
        for worker,examplelabels in w2el.items():
            w2cm[worker]=0.0
            for e,label in examplelabels:
                w2cm[worker]+=e2lpd[e][label]*1.0/len(examplelabels)

        return w2cm


    def Run(self, iter = 5, workers={}):
        l2pd = self.InitPj()
        w2cm = self.Initw2cm(workers)
        while iter>0:
            # E-step
            e2lpd = self.ComputeTij(self.e2wl, {}, w2cm)

            # M-step
            #l2pd = self.ComputePj(e2lpd)
            w2cm = self.Computew2cm(self.w2el, e2lpd)

            #print l2pd,w2cm

            iter -= 1

        return e2lpd, w2cm

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

if __name__=='__main__':
    # w2cm   worker_to_confusion_matrix = {}
    # e2pd   example_to_softlabel = {}
    # l2pd   label_to_priority_probability = {}
    datafile = sys.argv[1]
    e2wl,w2el,label_set = gete2wlandw2el(datafile)
    e2lpd, w2cm= ZenClowd(e2wl,w2el,label_set).Run()
    print(w2cm)
    print("e2lpd: ", e2lpd)

    truthfile = sys.argv[2]
    accuracy = getaccuracy(truthfile, e2lpd, label_set)
    print("accuracy: ", accuracy)


    