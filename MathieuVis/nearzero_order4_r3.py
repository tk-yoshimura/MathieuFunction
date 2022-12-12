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

def approx_term4(n, q, r0, r1, r2, r3):
    h = q * q

    n_sq = n * n
    n_sq_m1 = n_sq - 1
    n_sq_m4 = n_sq - 4
    n_sq_m9 = n_sq - 9
    n_sq_m16 = n_sq - 16

    c8 = (r0 + n_sq * (r1 + n_sq * (r2 + n_sq * r3))) / ((n_sq_m1**7) * n_sq_m4 * n_sq_m9 * n_sq_m16)

    if not np.all(np.isfinite(c8)):
        return np.nan

    y = h**4 * c8 / n_sq

    return y

def ffunc(x, r0, r1, r2, r3):
    n, q = x

    return approx_term4(n, q, r0, r1, r2, r3)

for func in ['A', 'B']:
    qs, ns, residuals = [], [], []

    for n in range(5, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../results/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, expected = data['q'], data['convergence']
        q, expected = q[len(q)//16:len(q)//8], expected[len(expected)//16:len(expected)//8]

        approx = approx_order3(n, q) 
        residual = (approx - expected) / (n * n)

        qs.append(q)
        ns.append(n)
        residuals.append(residual)    

    qs = np.stack(qs)
    ns = np.stack(ns).astype(float)[:, np.newaxis]
    ns = np.broadcast_to(ns, qs.shape)
    residuals = np.stack(residuals)

    with open('../results/eigen_%s_term4.csv' % (func), 'w') as f:

        popt, _ = curve_fit(ffunc, (ns.flatten(), qs.flatten()), residuals.flatten(), p0=(0, 0, 0, 1))

        r0, r1, r2, r3 = popt

        actuals = approx_term4(ns, qs, r0, r1, r2, r3)

        error = np.sum(np.abs(actuals - residuals))

        print(popt)
        f.write(popt)
        f.write(error)