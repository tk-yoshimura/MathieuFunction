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

def approx_term4(q, r0, r1, r2, r3):
    h = q * q

    y = q * (r0 + q * (r1 + q * (r2 + (q * r3))))

    return y

for func in ['A', 'B']:
    qs, ns, residuals = [], [], []

    for n in range(4, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../results/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, expected = data['q'], data['convergence']
        q, expected = q[len(q)//16:len(q)//8], expected[len(expected)//16:len(expected)//8]


        approx = approx_order3(float(n), q) 
        residual = (approx - expected) / (n * n)

        q /= (n * n)
        q -= 0.5

        popt, _ = curve_fit(approx_term4, q, residual, p0=(0, 0, 0, 0))
        r0, r1, r2, r3 = popt

        actual = approx_term4(q, r0, r1, r2, r3)

        plt.clf()
        plt.plot(q, residual, label='expected')
        plt.plot(q, actual, label='actual')

        plt.xlabel('q')

        plt.savefig('../figures/eigen_%s_%d_approx_term4.png' % (func, n))