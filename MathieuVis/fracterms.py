import numpy as np
import pandas as pd
import scipy.stats as scstat
import matplotlib.pyplot as plt

medianwidth = 8

for func in ['A', 'B']:

    ns, intercepts, slopes = [], [], [] 

    for n in range(0, 128+1):
        if func == 'B' and n == 0:
            continue

        print('plot %s %d' % (func, n))
        
        plt.clf()

        data = np.loadtxt('../results/needs_frac_%s_%d.csv' % (func, n), delimiter=',')

        q, terms = data[:, 0], data[:, 2]

        terms = np.pad(terms, (medianwidth, medianwidth), mode='edge')
        terms_sliced = []
        for i in range(0, medianwidth * 2 + 1):
            terms_sliced.append(terms[i:][:len(q)])

        terms = np.sort(np.stack(terms_sliced), axis=0)[medianwidth]

        x, y = np.log2(q), np.log2(terms)
        res = scstat.linregress(x[len(x)//4:], y[len(y)//4:])
        intercept, slope = res.intercept, res.slope
        u = 2**(intercept + slope * x)

        intercept_actual = 2.13317 + np.sqrt(n + 1) * 0.104619
        slope_actual = 0.239519 * (n + 1) ** 0.0192171
        v = 2**(intercept_actual + slope_actual * x)

        ns.append(n)
        intercepts.append(intercept)
        slopes.append(slope)

        plt.plot(q, terms, label='n=' + str(n), linewidth='1')
        plt.plot(q, u, label='approx(q) = pow2(%.4f + %.4f q)' % (intercept, slope), linewidth='1')
        plt.plot(q, v, label='approx(q) = pow2(%.4f + %.4f q)' % (intercept_actual, slope_actual), linewidth='1')
        
        plt.legend(loc='upper left')

        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('q')
        plt.ylabel('c.frac terms')

        plt.savefig('../figures/needs_frac_%s_%d.png' % (func, n))

    plt.clf()
    
    for n in 2**np.arange(0, 8):
        if func == 'B' and n == 0:
            continue

        print('plot %s %d' % (func, n))
    
        data = np.loadtxt('../results/needs_frac_%s_%d.csv' % (func, n), delimiter=',')

        q, terms = data[:, 0], data[:, 2]

        terms = np.pad(terms, (medianwidth, medianwidth), mode='edge')
        terms_sliced = []
        for i in range(0, medianwidth * 2 + 1):
            terms_sliced.append(terms[i:][:len(q)])

        terms = np.sort(np.stack(terms_sliced), axis=0)[medianwidth]

        plt.plot(q, terms, label='n=' + str(n), linewidth='1')
        
    plt.legend(loc='upper left')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('q')
    plt.ylabel('c.frac terms')

    plt.savefig('../figures/needs_frac_%s.png' % (func))

    plt.clf()
    plt.plot(ns, intercepts, label='intercept', linewidth='1')
    plt.legend(loc='upper left')
    
    plt.xlabel('n')

    plt.savefig('../figures/needs_frac_intercept_%s.png' % (func))

    plt.clf()
    plt.plot(ns, slopes, label='slopes', linewidth='1')
    plt.legend(loc='upper left')
    
    plt.xlabel('n')

    plt.savefig('../figures/needs_frac_slope_%s.png' % (func))

    approx_res = pd.DataFrame([ns, slopes, intercepts], index=['n', 'slope', 'intercept']).T
    approx_res.to_csv('../results/needs_frac_%s.csv' % (func), index=False)