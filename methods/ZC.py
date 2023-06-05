import copy
import sys
import random
import csv


class ZenCrowd:

    def __init__(self, e2wl, w2el, label_set, scale=1.0):
        self.e2wl = e2wl
        self.w2el = w2el
        self.label_set = label_set
        self.scale = scale

    def Initw2cm(self, workers):
        w2cm = {}

        if workers == {}:
            workers = self.w2el.keys()
            for worker in workers:
                w2cm[worker] = self.scale
        else:
            for worker in workers:
                if worker not in w2cm:  
                    w2cm[worker] = self.scale * 4//5
                else:
                    w2cm[worker] = workers[worker]

        return w2cm

    # E-step
    def ComputeTij(self, e2wl, l2pd, w2cm):
        e2lpd = {}
        for e, workerlabels in e2wl.items():
            e2lpd[e] = {}
            for label in self.label_set:
                e2lpd[e][label] = 1  # l2pd[label]

            for worker, label in workerlabels:
                for candlabel in self.label_set:
                    if label == candlabel:
                        e2lpd[e][candlabel] *= w2cm[worker]
                    else:
                        e2lpd[e][candlabel] *= (self.scale -
                                                w2cm[worker]) // (len(self.label_set)-1)

            sums = 0
            for label in self.label_set:
                sums += e2lpd[e][label]

            if sums == 0:
                for label in self.label_set:
                    e2lpd[e][label] = self.scale // len(self.label_set)
            else:
                for label in self.label_set:
                    e2lpd[e][label] = e2lpd[e][label] * self.scale // sums  # average

        # print e2lpd
        return e2lpd

    # M-step
    def Computew2cm(self, w2el, e2lpd):
        w2cm = {}
        for worker, examplelabels in w2el.items():
            w2cm[worker] = 0
            for e, label in examplelabels:
                w2cm[worker] += e2lpd[e][label] // len(examplelabels)

        return w2cm

    def Run(self, iter=5, workers={}):
        w2cm = self.Initw2cm(workers)

        while iter > 0:
            # E-step
            e2lpd = self.ComputeTij(self.e2wl, {}, w2cm)

            # M-step
            w2cm = self.Computew2cm(self.w2el, e2lpd)

            iter -= 1

        return e2lpd, w2cm


def get_accuracy(truthfile, e2lpd, label_set):
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


def get_data(datafile):
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

    # print(len(e2wl), len(w2el), len(label_set))
    return e2wl, w2el, label_set


if __name__ == '__main__':
    data_file = sys.argv[1]
    e2wl, w2el, label_set = get_data(data_file)

    e2lpd, w2cm = ZenCrowd(e2wl, w2el, label_set, scale = 1).Run()
    # print("w2cm: ", w2cm)
    print(e2lpd)

    # truthfile = sys.argv[2]
    # accuracy = get_accuracy(truthfile, e2lpd, label_set)
    # print("accuracy: ", accuracy)
