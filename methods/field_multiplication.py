from Crypto.Util.number import *
from tqdm import tqdm

p = getPrime(96)  # bit prime
sa = 3450065664
sb = 2791728640
sc = 2242542599

c = 32
c_1 = 2**32
c_2 = 2**31
mul = sa * sb
delta = 2**31

print("start: ")
# for i in tqdm(range(2**31, 2**34)): 
i = 11271428967
if (((delta * ((mul - i*c_1) % p)) % p) < mul):
    print("c: {}, i: {}".format(c_1, i))
elif (((delta * ((mul - i*c_2) % p)) % p) < mul):
    print("c: {}, i: {}".format(c_2, i))
print("end: ")