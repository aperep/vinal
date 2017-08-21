# taken from https://gist.github.com/sirodoht/ee2abe82eca70f5b1869
from operator import mul, mod

def egcd(a, b):
    if a == 0:
        return (b, 0, 1)
    else:
        g, y, x = egcd(b % a, a)
        return (g, x - (b // a) * y, y)


def crt(m, a):
    M = reduce(mul, m) # the product of m elements
    m_i = [M / item for item in m]
    b = map(mod, m_i, m)
    g, k, l = map(egcd, b, m)
    g, k, l = zip(g, k, l) # transpose g, k and l arrays
    t = map(mod, k, m)
    e = map(mul, m_i, t)
    
    x_sum = sum(map(mul, a, e))
    x = x_sum % M
    return x
