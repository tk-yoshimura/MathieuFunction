import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

np.seterr(all='ignore')

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

def limit(func, n, q):
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

    if func == 'A':
        c6 = -5.682576740891 * (n**4.7) - 0.005055889635 * (n**7.9)
    else:
        c6 = 0.517862332643 * (n**5.5) - 0.002906986648 * (n**8.0)

    y = (2 * (-q + s * u) - n * n) + c0 + v * (c1 + v * (c2 + v * (c3 + v * (c4 + v * (c5 + v * c6)))))

    return y

def bump(x, s, t, a, b):
    c = (x - s) / (t - s);

    w = 1 / (np.exp(1 / c - 1 / (1 - c)) + 1)
    y = a * (1 - w) + b * w

    y = np.where(c < 0.001, a, np.where(c > 0.998, b, y))

    return y

def ffuncA(x, s, t):
    n, q = x
    m = np.maximum(1, n*n)
    
    y = bump(q / m, s, t, nearzero('A', n, q), limit('A', n, q))

    return y

def ffuncB(x, s, t):
    n, q = x
    m = np.maximum(1, n*n)
    
    y = bump(q / m, s, t, nearzero('B', n, q), limit('B', n, q))

    return y

for func in ['A', 'B']:
    ffunc = ffuncA if func == 'A' else ffuncB

    pss, pts, ss = [], [], []

    for n in range(4, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, raw, expected = data['q'].astype(float), data['approx'], data['convergence']

        ns = np.broadcast_to(n, len(q)).astype(float)
        
        if func == 'A':
            s = max(0.25, 0.500 - 8.492 * (n**-1.307))
            t = min(1.00, 0.763 + 4.734 * (n**-1.173)) 
        else:
            s = min(0.60, 0.442 + 9.609e+3 * (n**-3.299))
            t = max(0.80, 0.817 - 9.282e+2 * (n**-2.573))

        approx = ffunc((ns, q), s, t)

        sn = (2 * n + 1) if func == 'A' else (2 * n - 1)

        pss.append(s)
        pts.append(t)
        ss.append(sn)

        plt.clf()
        plt.plot(q[:len(q)//4], expected[:len(q)//4], label='expected', linewidth = 1, alpha=0.5)
        plt.plot(q[:len(q)//4], raw[:len(q)//4], label='raw', linewidth = 1, alpha=0.5)
        #plt.plot(q, nz, label='nearzero')
        #plt.plot(q, asymp, label='asymptotic')
        plt.plot(q[:len(q)//4], approx[:len(q)//4], label='approx\n bump(%.4f, %.4f)' % (s, t), linewidth = 1, alpha=0.5)
        
        plt.legend(loc='lower left')

        plt.xlabel('q')

        plt.savefig('../sandbox/eigen_%s_%d_approx_asymp.png' % (func, n))

    approx_res = pd.DataFrame([ss, pss, pts], index=['s', 'ps', 'pt']).T
    approx_res.to_csv('../sandbox/asymp_st_%s_r2.csv' % (func), index=False)