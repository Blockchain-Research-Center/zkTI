class Float:
    def __init__(self,v,scale,value):
        self.scale = scale
        self.v = v
        self.value = value

    def __mul__(self,other):
        scale = self.scale + other.scale
        v = self.v * other.v
        while v > 2**32:
            scale -= 1
            v = v // 2
        return Float(v,scale,self.value * other.value)

    def __add__(self,other):
        v = 0
        scale = 0
        if(self.scale == other.scale):
            scale = self.scale
            v = self.v + other.v
            while v > 2**32:
                scale -= 1
                v = v // 2
            return Float(v,scale,self.value + other.value)
        else:
            if abs(self.scale - other.scale)>=10:
                if(self.scale < other.scale):
                    return Float(self.v,self.scale,self.value + other.value)
                else:
                    return Float(other.v,other.scale,self.value + other.value)
            else:
                if(self.scale < other.scale):
                    bigs = self.scale
                    bigv = self.v
                    smalls = other.scale
                    smallv = other.v
                else:
                    smalls = self.scale
                    smallv = self.v
                    bigs = other.scale
                    bigv = other.v

                v = 2 ** (smalls-bigs) * bigv + smallv
                scale = smalls
                while v > 2**32:
                    scale -= 1
                    v = v // 2

            return Float(v,scale,self.value + other.value)

    def __str__(self):
        return "M: {0}*2^-{1}, R: {2}".format(self.v, self.scale, self.value)

    def __truediv__(self,other):
        if isinstance(other,int) == True or isinstance(other,float) == True:
            v = int(self.v / other)
            scale = self.scale
            while v > 2**32:
                scale -= 1
                v = v // 2
            return Float(v,scale,self.value/other)
        else:
            v = self.v * 2**32 / other.v
            scale = self.scale + 32 - other.scale
            while v > 2**32:
                scale -= 1
                v = v // 2
            return Float(v,scale,self.value/other.value)
    
    def one_sub(self):
        if self.scale >= 32+10:
            return Float(2**31,31,1.0-self.value)
        else:
            v = 2**31 * 2 **(self.scale - 31) - self.v
            scale = self.scale
            while v > 2**32:
                scale -= 1
                v = v // 2
            return Float(v,scale,1.0 - self.value)


# f1 = Float(2570830560,24,2.570830560230192e-15)
# f2 = Float(9999999999,10,0.9999999999999973)
# f = Float(8000000000,10,0.8)
# f3 = f1 + f2
# f4 = f1/(f3)
# print(f.one_sub())
# f1 = Float(2**31,32, 0.5)
# f2 = Float(2**31,45,0.25)
# print(f2.one_sub())
# print(f1+f2)