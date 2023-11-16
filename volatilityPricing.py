from numpy import exp, sqrt
from math import comb
import sys

def euroCallOption(S0, r, sigma, K, T, N, rq):
    sum = 0
    dt = T / N
    u = exp(sigma * sqrt(dt))
    d = 1 / u
    p = (exp((r - rq)* dt) - d) / (u - d)   
    q = 1 - p
    for j in range(0, N + 1):
        c = max(0, S0 * u**j * d**(N - j)- K)
        sum += comb(N, j) * p**(j) * q**(N - j) * c
    return exp(-r * T) * sum, p, u, d

def euroPutOption(S0, r, sigma, K, T, N, rq):
    sum = 0
    dt = T / N
    u = exp(sigma * sqrt(dt))
    d = 1 / u
    p = (exp((r - rq) * dt) - d) / (u - d)
    q = 1 - p
    for j in range(0, N + 1):
        c = max(K - S0 * u**j * d**(N - j), 0)
        sum += comb(N, j) * p**(j) * q**(N - j) * c
    return exp(-r * T) * sum, p, u, d


def americanCallOption(S0, r, sigma, K, T, N, rq):
    dt = T / N
    u = exp(sigma * sqrt(dt))
    d = 1 / u
    p = (exp((r - rq) * dt) - d) / (u - d)
    q = 1 - p
    fSave = []
    discount = exp(-(r - rq) * dt)
    for i in range(N, -1, -1):
        f = []
        for j in range(0, i + 1):
            if i == N:
                f.append(max(S0 * d**j * u**(i - j) - K, 0))
            else:
                fu = fSave[j]
                fd = fSave[j + 1]
                f.append(max(S0 * d**j * u**(i - j) - K, discount * (p * fu + q * fd)))
        fSave = f
    return fSave[0]

def americanPutOption(S0, r, sigma, K, T, N, rq):
    dt = T / N
    print(T, N, sigma)
    u = exp(sigma * sqrt(dt))
    d = 1 / u
    p = (exp((r - rq) * dt) - d) / (u - d)
    print(p, d, u, exp(r * dt))
    q = 1 - p
    fSave = []
    discount = exp(-(r - rq) * dt)
    for i in range(N, -1, -1):
        f = []
        for j in range(0, i + 1):
            if i == N:
                f.append(max(K - S0 * d**j * u**(i - j), 0))
            else:
                fu = fSave[j]
                fd = fSave[j + 1]
                f.append(max(K - S0 * d**j * u**(i - j), discount * (p * fu + q * fd)))
        fSave = f
    return fSave[0]


if __name__ == "__main__":
    try:
        method = sys.argv[1]
        S0 = float(sys.argv[2])
        r = float(sys.argv[3])
        sigma = float(sys.argv[4])
        K = float(sys.argv[5])
        T = float(sys.argv[6])
        N = int(sys.argv[7])
        q = 0
        if len(sys.argv) > 8:
            q = float(sys.argv[8])

    except:
        print("Please call the function as follows: >> python3 Pricing.py 'method' 'S0' 'r' 'sigma' 'K' 'T' 'N' 'q-(optional)")
        print("Example: >> python3 Pricing.py call 35 0.05 0.22 40 0.5 2")

    if method == 'call':
        price, p, u, d = euroCallOption(S0, r, sigma, K, T, N, q)
        APrice = americanCallOption(S0, r, sigma, K, T, N, q)
        print('The risk-neutral probability is: ', p)
        print('The upturn u is: ', u)
        print('The downturn d is: ', d)
        print('The price of the European call option is: ', price)
        print('The price of the American call option is: ', APrice)
    elif method == 'put':
        price, p, u, d = euroPutOption(S0, r, sigma, K, T, N, q)
        APrice = americanPutOption(S0, r, sigma, K, T, N, q)
        print('The risk-neutral probability is: ', p)
        print('The upturn u is: ', u)
        print('The downturn d is: ', d)
        print('The price of the European put option is: ', price)
        print('The price of the American put option is: ', APrice)
    else:
        raise Exception("Please put either 'call' or 'put' for method")




