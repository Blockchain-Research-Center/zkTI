import os
import subprocess

def exe(cmd):
    result = os.popen(cmd).read()
    return result

r1 = exe("python3 ZC.py ../datasets/d_duck_identification/answer.csv ../datasets/d_duck_identification/truth.csv")

r1 = eval(r1)


r2 = exe("python3 ZenCrowd.py ../datasets/d_duck_identification/answer.csv ../datasets/d_duck_identification/truth.csv")

r2 = eval(r2)

for key in r1.keys():
    print(r1[key]['0']/100000000-r2[key]['0'])