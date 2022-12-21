import numpy as np
import pandas as pd
import itertools
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

def ffunc(x, s0, p0, s1, p1):
    y = s0 * (x**p0) + s1 * (x**p1)
    
    return y

def ffunc_log(x, s0, p0, s1, p1):
    y = ffunc(x, s0, p0, s1, p1)
    y = np.log10(-y)
    
    return y

for func in ['A', 'B']:
    data = pd.read_csv('../sandbox/asymp_c6_%s.csv' % func, delimiter=',')

    n, expected = data['n'].astype(float), data['c6'].astype(float)
    best_s0, best_p0, best_s1, best_p1 = -8, 4.5, -1e-3, 8.0
    min_error = float('inf')

    n, expected = n[expected < 0], expected[expected < 0]

    s0_init = -1e-3 if func == 'A' else 1e-3
    for p0, p1 in itertools.product(np.linspace(-2, 6, 81), np.linspace(6, 10, 41)):
        if p0 >= p1:
            continue

        def ffunc_s(x, s0, s1):
            y = ffunc(x, s0, p0, s1, p1)
            y = np.log10(-y)
    
            return y

        try:
            popt, _ = curve_fit(ffunc_s, n, np.log10(-expected), p0=(s0_init, -1e-3))
        except:
            continue

        s0, s1 = popt

        approx = ffunc(n, s0, p0, s1, p1)

        error = 10**np.max(np.abs(np.log10(expected / approx)))

        if min_error > error:
            min_error = error
            best_s0, best_p0, best_s1, best_p1 = s0, p0, s1, p1

            print(min_error)
            print('s0=%.12f, p0=%.12f s1=%.12f, p1=%.12f' % (best_s0, best_p0, best_s1, best_p1))

    try:
        ds0, dp0, ds1, dp1 = np.abs(best_s0 * 0.01), np.abs(best_p0 * 0.01), np.abs(best_s1 * 0.01), np.abs(best_p1 * 0.01)
 
        bounds=((best_s0 - ds0, best_p0 - dp0, best_s1 - ds1, best_p1 - dp1),  
                (best_s0 + ds0, best_p0 + dp0, best_s1 + ds1, best_p1 + dp1))
        popt, _ = curve_fit(ffunc_log, n, np.log10(-expected), p0=(best_s0, best_p0, best_s1, best_p1), bounds=bounds)
        s0, p0, s1, p1 = popt

        approx = ffunc(n, s0, p0, s1, p1)

        error = 10**np.max(np.abs(np.log10(expected / approx)))

        if min_error > error:
            min_error = error
            best_s0, best_p0, best_s1, best_p1 = s0, p0, s1, p1
    except:
        continue

    approx = ffunc(n, best_s0, best_p0, best_s1, best_p1)

    error = 10**np.max(np.abs(np.log10(expected / approx)))
    print(error)
    print('s0=%.12f, p0=%.12f s1=%.12f, p1=%.12f' % (best_s0, best_p0, best_s1, best_p1))

    plt.clf()
    plt.plot(n, -expected, label='expected', linewidth = 1, alpha=0.5)
    plt.plot(n, -approx, label='approx\n s0=%.12f, p0=%.12f\n s1=%.12f, p1=%.12f' % (best_s0, best_p0, best_s1, best_p1), linewidth = 1, alpha=0.5)
    plt.xlabel('n')
    plt.xscale('log')
    plt.yscale('log')

    plt.legend(loc='upper left')

    plt.savefig('../sandbox/approx_error_c6_%s.png' % (func))
