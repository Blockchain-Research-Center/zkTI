from Crypto.Util.number import *
from tqdm import tqdm

p = getPrime(128)  # 24 bit prime
sa = 3450065664
sb = 2791728640
sc = 2242542599

c = 32
delta = 2**31

# assert(2**32*(sa*sb-sc*2**32)<sa*sb)

for i in tqdm(range(2**31, 2**32)):
    if ((delta * ((sa*sb - i*2**c) % p) % p) < sa*sb and 2**31 < i < 2**32):
        print("c: {}, i: {}".format(c, i))
    elif ((delta * ((sa*sb - i*2**(c-1)) % p) % p) < sa*sb and 2**31 < i < 2**32):
        print("c: {}, i: {}".format(c-1, i))

# for i in tqdm(range(0, p)):
#     if((delta * ((sa*sb - i*2**c)%p)%p)<sa*sb and 2**31<i<2**32):
#         print("c: %d, i: %d", c, i)
#     elif((delta * ((sa*sb - i*2**(c-1))%p)%p)<sa*sb and 2**31<i<2**32):
#         print("c: %d, i: %d", c-1, i)
