import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

np.seterr(all='ignore')

def nearzero(n, q):
    h = q * q

    n_sq = n * n
    n_sq_m1 = n_sq - 1
    n_sq_m4 = n_sq - 4
    n_sq_m9 = n_sq - 9

    c2 = 1 / (2 * n_sq_m1)
    c4 = (7 + n_sq * 5) / (32 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4)
    c6 = (29 + n_sq * (58 + n_sq * 9)) / (64 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4 * n_sq_m9)

    y = h * (c2 + h * (c4 + h * c6))

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

def ffuncA(x, s, t, c6):
    n, q = x

    y = bump(q / (n * n), s, t, nearzero(n, q), limit('A', n, q, c6 * n*n * 1e6))

    return y

def ffuncB(x, s, t, c6):
    n, q = x

    y = bump(q / (n * n), s, t, nearzero(n, q), limit('B', n, q, c6 * n*n * 1e6))

    return y

for func in ['A', 'B']:
    ffunc = ffuncA if func == 'A' else ffuncB

    pss, pts, c6s, ss = [], [], [], []

    for n in range(4, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../results/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, raw, expected = data['q'].astype(float), data['approx'], data['convergence']

        ns = np.broadcast_to(n, len(q)).astype(float)

        #nz = nearzero(ns, q)
        #asymp = limit(func, ns, q)

        min_error = float('inf')
        best_s, best_t = 0.5, 0.75

        for s_init, t_init in itertools.product(np.linspace(0.1, 1, 10), np.linspace(0.1, 1, 10)):
            if s_init >= t_init:
                continue
            
            popt, _ = curve_fit(ffunc, (ns, q), expected, p0=(s_init, t_init, 0.0))
            s, t, c6 = popt
            approx = ffunc((ns, q), s, t, c6)

            error = np.sum(np.abs(expected - approx))

            if error < min_error:
                min_error = error
                best_s, best_t = s, t

        s, t = best_s, best_t
        approx = ffunc((ns, q), s, t, c6)

        sn = (2 * n + 1) if func == 'A' else (2 * n - 1)
        c6 *= n*n * 1e6

        pss.append(s)
        pts.append(t)
        ss.append(sn)
        c6s.append(c6)

        plt.clf()
        plt.plot(q[:len(q)//4], expected[:len(q)//4], label='expected', linewidth = 1, alpha=0.5)
        plt.plot(q[:len(q)//4], raw[:len(q)//4], label='raw', linewidth = 1, alpha=0.5)
        #plt.plot(q, nz, label='nearzero')
        #plt.plot(q, asymp, label='asymptotic')
        plt.plot(q[:len(q)//4], approx[:len(q)//4], label='approx\n bump(%.4f, %.4f, %.4f)' % (s, t, c6), linewidth = 1, alpha=0.5)
        
        plt.legend(loc='lower left')

        plt.xlabel('q')

        plt.savefig('../figures/eigen_%s_%d_approx_asymp.png' % (func, n))

    approx_res = pd.DataFrame([ss, pss, pts, c6s], index=['s', 'ps', 'pt', 'c6']).T
    approx_res.to_csv('../results/asymp_stc6_%s.csv' % (func), index=False)