import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

np.seterr(all='ignore')

def square_frac(n):
    x = 1.0;
    for i in range(2, n + 1):
        x *= i

    return x * x

def nearzero(func, n, q):
    assert np.min(n) == np.max(n), 'invalid n'

    n = int(n[0])

    h = q * q
    
    if int == 0:
        mean = h * (-1 / 2 + h * (7 / 128 + h * (-29 / 2304 + h * (68687 / 18874368))))
        sub = 0
    elif n == 1:
        mean = h * (-1 / 8 + h * (-1 / 1536 + h * (49 / 589824 + h * (-83 / 35389440))))
        sub = q * (1 + h * (-1 / 64 + h * (11 / 36864 + h * (55 / 9437184))))
    elif n == 2:
        mean = h * (1 / 6 + h * (-379 / 13824 + h * (7829 / 1244160)))
        sub = q**2 * (1 / 4 + h * (-1 / 36 + h * (11141 / 1769472)))
    elif n == 3:
        mean = h * (1 / 16 + h * (13 / 20480 + h * (-1961 / 23592960)))
        sub = q**3 * (1 / 64)
    elif n == 4:
        mean = h * (1 / 30 + h * (29 / 432000 + h * (1087 / 1360800000)))
        sub = 0
    else:
        if n == 5:
            mean = h * (1 / 48 + h * (11 / 774144 + h * (37 / 891813888)))
        elif n == 6:
            mean = h * (1 / 70 + h * (187 / 43904000 + h * (13781 / 2904249600000)))
        else:
            n = float(n)

            n_sq = n * n
            n_sq_m1 = n_sq - 1
            n_sq_m4 = n_sq - 4
            n_sq_m9 = n_sq - 9

            c2 = 1 / (2 * n_sq_m1)
            c4 = (7 + n_sq * 5) / (32 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4)
            c6 = (29 + n_sq * (58 + n_sq * 9)) / (64 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4 * n_sq_m9)

            mean = h * (c2 + h * (c4 + h * c6))
        sub = 0

    y = mean + (sub if func == 'A' else -sub)

    return y

def limit(func, n, q, c6=0):
    u = np.sqrt(q)
    v = 1 / u
    s = (2 * n + 1) if func == 'A' else (2 * n - 1)
    s_sq = s * s

    c0 = -np.ldexp(1 + s_sq, -3)
    c1 = -np.ldexp(s * (3 + s_sq), -7)
    c2 = -np.ldexp(9 + s_sq * (34 + s_sq * 5), -12)
    c3 = -np.ldexp(s * (405 + s_sq * (410 + s_sq * 33)), -17)
    c4 = -np.ldexp(486 + s_sq * (2943 + s_sq * (1260 + s_sq * 63)), -20)
    c5 = -np.ldexp(s * (41607 + s_sq * (69001 + s_sq * (15617 + s_sq * 527))), -25)

    y = (2 * (-q + s * u) - n * n) + c0 + v * (c1 + v * (c2 + v * (c3 + v * (c4 + v * (c5 + v * c6)))))

    return y

def bump(x, s, t, a, b):
    c = (x - s) / (t - s);

    w = 1 / (np.exp(1 / c - 1 / (1 - c)) + 1)
    y = a * (1 - w) + b * w

    y = np.where(c < 0.001, a, np.where(c > 0.998, b, y))

    return y

def ffuncA(x, s, t):
    n, q, c6 = x
    m = np.maximum(1, n*n)
    
    y = bump(q / m, s, t, nearzero('A', n, q), limit('A', n, q, c6 * m * 1e6))

    return y

def ffuncB(x, s, t):
    n, q, c6 = x
    m = np.maximum(1, n*n)
    
    y = bump(q / m, s, t, nearzero('B', n, q), limit('B', n, q, c6 * m * 1e6))

    return y

def flimitA(x, c6):
    n, q = x
    m = np.maximum(1, n*n)
    
    y = limit('A', n, q, c6 * m * 1e6)

    return y

def flimitB(x, c6):
    n, q = x
    m = np.maximum(1, n*n)

    y = limit('B', n, q, c6 * m * 1e6)

    return y

for func in ['A', 'B']:
    ffunc = ffuncA if func == 'A' else ffuncB
    flimit = flimitA if func == 'A' else flimitB

    c6s, nss = [], []

    c6_init = 0

    for n in range(4, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, raw, expected = data['q'].astype(float), data['approx'], data['convergence']

        peak = np.argmax(expected)
        thr = peak * 3 // 4

        ns = np.broadcast_to(n, len(q)).astype(float)

        #nz = nearzero(ns, q)
        #asymp = limit(func, ns, q)

        popt, _ = curve_fit(flimit, (ns[thr:], q[thr:]), expected[thr:], p0=(c6_init))
        c6, = popt
        c6_init = c6

        print(c6)
        
        asymp = flimit((ns, q), c6)
        
        c6 *= np.maximum(1, n*n) * 1e6
        
        nss.append(n)
        c6s.append(c6)

        plt.clf()
        plt.plot(q[:len(q)//4], expected[:len(q)//4], label='expected', linewidth = 1, alpha=0.5)
        plt.plot(q[:len(q)//4], raw[:len(q)//4], label='raw', linewidth = 1, alpha=0.5)
        plt.plot(q[thr:len(q)//4], asymp[thr:len(q)//4], label='asymptotic\n c6=%.4f' % c6, linewidth = 1, alpha=0.5)
        #plt.plot(q, nz, label='nearzero')
        #plt.plot(q, asymp, label='asymptotic')
        plt.legend(loc='lower left')

        plt.xlabel('q')

        plt.savefig('../sandbox/eigen_%s_%d_approx_asymp_c6.png' % (func, n))

    approx_res = pd.DataFrame([nss, c6s], index=['n', 'c6']).T
    approx_res.to_csv('../sandbox/asymp_c6_%s.csv' % (func), index=False)