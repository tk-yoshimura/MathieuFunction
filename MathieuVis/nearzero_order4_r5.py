import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

def approx_order3(n, q):
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

def approx_term4(n, q, r0, r1, r2, r3, r4, r5, r6, r7, r8):
    h = q * q

    n_sq = n * n
    n_sq_m1 = n_sq - 1
    n_sq_m4 = n_sq - 4
    n_sq_m9 = n_sq - 9

    c2 = (r0 + n_sq * r1) / (2 * n_sq_m1)
    c4 = (r2 + n_sq * (r3 + n_sq * r4)) / (32 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4)
    c6 = (r5 + n_sq * (r6 + n_sq * (r7 + n_sq * r8))) / (64 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4 * n_sq_m9)
    
    y = h * (c2 + h * (c4 + h * c6)) / n_sq

    return y

def ffunc(x, r0, r1, r2, r3, r4, r5, r6, r7, r8):
    n, q = x

    return approx_term4(n, q, r0, r1, r2, r3, r4, r5, r6, r7, r8)

for func in ['A', 'B']:
    qs, ns, residuals = [], [], []

    for n in range(5, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../results/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, expected = data['q'].astype(float), data['convergence']
        q, expected = q[:len(q)//8], expected[:len(expected)//8]

        approx = approx_order3(float(n), q) 
        residual = (approx - expected) / (n * n)

        plt.clf()
        plt.plot(q, expected, label='expected')
        plt.plot(q, approx, label='approx')

        plt.xlabel('q')

        plt.savefig('../figures/eigen_%s_%d_approx_term3.png' % (func, n))

        qs.append(q)
        ns.append(n)
        residuals.append(residual)    

    qs = np.stack(qs)
    ns = np.stack(ns).astype(float)[:, np.newaxis]
    ns = np.broadcast_to(ns, qs.shape)
    residuals = np.stack(residuals)

    with open('../results/eigen_%s_term4.csv' % (func), 'w') as f:

        popt, _ = curve_fit(ffunc, (ns.flatten(), qs.flatten()), residuals.flatten(), p0=(0, 0, 0, 0, 0, 0, 0, 0, 0))

        r0, r1, r2, r3, r4, r5, r6, r7, r8 = popt

        actuals = approx_term4(ns, qs, r0, r1, r2, r3, r4, r5, r6, r7, r8)

        error = np.sum(np.abs(actuals - residuals))

        print(popt)
        f.write(str(popt))
        f.write(str(error))

    for n in range(5, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('plot %s %d' % (func, n))

        data = pd.read_csv('../results/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, expected = data['q'].astype(float), data['convergence']
        q, expected = q[:len(q)//8], expected[:len(expected)//8]

        approx = approx_order3(float(n), q) 
        residual = (approx - expected) / (n * n)
        actual = approx_term4(float(n), q, r0, r1, r2, r3, r4, r5, r6, r7, r8)

        plt.clf()
        plt.plot(q, residual, label='expected')
        plt.plot(q, actual, label='actual')

        plt.xlabel('q')

        plt.savefig('../figures/eigen_%s_%d_approx_term4.png' % (func, n))