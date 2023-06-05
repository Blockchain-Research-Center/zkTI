import csv
import math
import random
import sys
import numpy as np
from scipy.optimize import minimize

class GLAD:

    def __init__(self, e2wl, w2el, label_set):
        self.e2wl = e2wl
        self.w2el = w2el
        self.workers = self.w2el.keys()
        self.examples = self.e2wl.keys()
        self.label_set = label_set

    def sigmoid(self, x):
        if (-x) > math.log(sys.float_info.max):
            return 0
        if (-x) < math.log(sys.float_info.min):
            return 1

        return 1 / (1 + math.exp(-x))

    def logsigmoid(self, x):
        # For large negative x, -log(1 + exp(-x)) = x
        if (-x) > math.log(sys.float_info.max):
            return x
        # For large positive x, -log(1 + exp(-x)) = 0
        if (-x) < math.log(sys.float_info.min):
            return 0

        value = -math.log(1 + math.exp(-x))
        # if (math.isinf(value)):
        #    return x

        return value

    def logoneminussigmoid(self, x):
        # For large positive x, -log(1 + exp(x)) = -x
        if (x) > math.log(sys.float_info.max):
            return -x
        # For large negative x, -log(1 + exp(x)) = 0
        if (x) < math.log(sys.float_info.min):
            return 0

        value = -math.log(1 + math.exp(x))
        # if (math.isinf(value)):
        #    return -x

        return value

    def kronecker_delta(self, answer, label):
        if answer == label:
            return 1
        else:
            return 0

    def expbeta(self, beta):
        if beta >= math.log(sys.float_info.max):
            return sys.float_info.max
        else:
            return math.exp(beta)

    # E step
    def Update_e2lpd(self):
        self.e2lpd = {}
        for example, worker_label_set in self.e2wl.items():
            lpd = {}
            total_weight = 0

            for tlabel, prob in self.prior.items():
                weight = math.log(prob)
                for (worker, label) in worker_label_set:
                    # log[p(Lij=Zj|alpha,beta)]
                    logsigma = self.logsigmoid(
                        self.alpha[worker] * self.expbeta(self.beta[example]))
                    # log[1-p(Lij=Zj|alpha,beta)]
                    logoneminussigma = self.logoneminussigmoid(
                        self.alpha[worker] * self.expbeta(self.beta[example]))
                    delta = self.kronecker_delta(label, tlabel)
                    weight = weight + delta * logsigma + \
                        (1 - delta) * (logoneminussigma -
                                       math.log(len(label_set) - 1))
                    # weight = weight + delta * logsigma + (1 - delta) * logoneminussigma

                if weight < math.log(sys.float_info.min):
                    lpd[tlabel] = 0
                else:
                    lpd[tlabel] = math.exp(weight)
                total_weight = total_weight + lpd[tlabel]

            for tlabel in lpd:
                if total_weight == 0:
                    lpd[tlabel] = 1.0 / len(self.label_set)
                else:
                    lpd[tlabel] = lpd[tlabel] * 1.0 / total_weight

            self.e2lpd[example] = lpd

    # M step
    def gradientQ(self):

        self.dQalpha = {}
        self.dQbeta = {}

        for example, worker_label_set in self.e2wl.items():
            dQb = 0
            for (worker, label) in worker_label_set:
                for tlabel in self.prior.keys():
                    sigma = self.sigmoid(
                        self.alpha[worker] * self.expbeta(self.beta[example]))
                    delta = self.kronecker_delta(label, tlabel)
                    dQb = dQb + self.e2lpd[example][tlabel] * (
                        (delta - sigma) * self.alpha[worker] + (1 - delta) * math.log(len(label_set)-1))
            # self.dQbeta[example] = dQb - (self.beta[example] - self.priorbeta[example])
            self.dQbeta[example] = dQb

        for worker, example_label_set in self.w2el.items():
            dQa = 0
            for (example, label) in example_label_set:
                for tlabel in self.prior.keys():
                    sigma = self.sigmoid(
                        self.alpha[worker] * self.expbeta(self.beta[example]))
                    delta = self.kronecker_delta(label, tlabel)
                    dQa = dQa + self.e2lpd[example][tlabel] * ((delta - sigma) * self.expbeta(
                        self.beta[example]) + (1 - delta) * math.log(len(label_set) - 1))
            # self.dQalpha[worker] = dQa - (self.alpha[worker] - self.prioralpha[worker])
            self.dQalpha[worker] = dQa

    def computeQ(self):

        Q = 0
        # the expectation of examples given priors, alpha and beta
        for worker, example_label_set in self.w2el.items():
            for (example, label) in example_label_set:
                logsigma = self.logsigmoid(
                    self.alpha[worker]*self.expbeta(self.beta[example]))
                logoneminussigma = self.logoneminussigmoid(
                    self.alpha[worker]*self.expbeta(self.beta[example]))
                for tlabel in self.prior.keys():
                    delta = self.kronecker_delta(label, tlabel)
                    Q = Q + self.e2lpd[example][tlabel]*(delta*logsigma+(
                        1-delta)*(logoneminussigma-math.log(len(label_set)-1)))

        # the expectation of the sum of priors over all examples
        for example in self.e2wl.keys():
            for tlabel, prob in self.prior.items():
                Q = Q + self.e2lpd[example][tlabel] * math.log(prob)

        # Gaussian (standard normal) prior for alpha
        # for worker in self.w2el.keys():
        #     Q = Q + math.log(
        #         (pow(2 * math.pi, -0.5)) * math.exp(-pow((self.alpha[worker] - self.prioralpha[worker]), 2) / 2))
        # # Gaussian (standard normal) prior for beta
        # for example in self.e2wl.keys():
        #     Q = Q + math.log(
        #         (pow(2 * math.pi, -0.5)) * math.exp(-pow((self.beta[example] - self.priorbeta[example]), 2) / 2))

        return Q

    def optimize_f(self, x):
        # unpack x
        i = 0
        for worker in self.workers:
            self.alpha[worker] = x[i]
            i = i + 1
        for example in self.examples:
            self.beta[example] = x[i]
            i = i + 1

        return -self.computeQ()  # Flip the sign since we want to minimize

    def optimize_df(self, x):
        # unpack x
        i = 0
        for worker in self.workers:
            self.alpha[worker] = x[i]
            i = i + 1
        for example in self.examples:
            self.beta[example] = x[i]
            i = i + 1

        self.gradientQ()

        # pack x
        der = np.zeros_like(x)
        i = 0
        for worker in self.workers:
            # Flip the sign since we want to minimize
            der[i] = -self.dQalpha[worker]
            i = i + 1
        for example in self.examples:
            # Flip the sign since we want to minimize
            der[i] = -self.dQbeta[example]
            i = i + 1

        return der

    def Update_alpha_beta(self):

        x0 = []
        for worker in self.workers:
            x0.append(self.alpha[worker])
        for example in self.examples:
            x0.append(self.beta[example])

        minimize(self.optimize_f, x0, method='CG', jac=self.optimize_df, tol=0.0001,
                 options={'disp': True, 'maxiter':  100})
        
        # test gradient
        x = []
        for worker in self.workers:
            x.append(self.alpha[worker])
        for example in self.examples:
            x.append(self.beta[example])
        
        print(self.optimize_df(x))

    # initialization
    def Init_prior(self):
        # uniform probability distribution
        prior = {}
        for label in self.label_set:
            prior[label] = 1.0 / len(self.label_set)
        return prior

    def Init_alpha_beta(self):
        prioralpha = {}
        priorbeta = {}
        for worker in self.w2el.keys():
            prioralpha[worker] = 1
        for example in self.e2wl.keys():
            priorbeta[example] = 1
        return prioralpha, priorbeta

    def Run(self, threshold=1e-5):

        self.prior = self.Init_prior()
        self.prioralpha, self.priorbeta = self.Init_alpha_beta()

        self.alpha = self.prioralpha
        self.beta = self.priorbeta

        Q = 0
        # E-step
        # self.Update_e2lpd()
        # Q = self.computeQ()

        while True:
            # E-step
            self.Update_e2lpd()
            Q = self.computeQ()
            print(Q)

            lastQ = Q

            # M-step
            self.Update_alpha_beta()
            Q = self.computeQ()
            print(Q)

            if (math.fabs((Q - lastQ) / lastQ)) < threshold:
                break

        return self.e2lpd, self.alpha


###################################
# The above is the EM method (a class)
# The following are several external functions
###################################

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

    return tcount * 1.0 / count


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


if __name__ == "__main__":
    datafile = sys.argv[1]
    e2wl, w2el, label_set = gete2wlandw2el(datafile)
    e2lpd, weight = GLAD(e2wl, w2el, label_set).Run(1e-4)

    print("weight:", weight)
    print("e2lpd:", e2lpd)

    truthfile = sys.argv[2]
    accuracy = getaccuracy(truthfile, e2lpd, label_set)
    print("accuracy", accuracy)
