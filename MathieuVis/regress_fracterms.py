import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 

def ffunc(x, a, b, c):
    x = np.log2(x)
    y = a + b * x + c * x * x
    y = 2**y
    
    return y

medianwidth = 1

for length in [4, 8, 16]:
    plt.clf()
    plt.figure(figsize=(12, 6))

    for func in ['A', 'B']:
        for j, n in enumerate([0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 16, 32, 64, 128, 256]):
            if func == 'B' and n == 0:
                continue

            print('plot %s %d' % (func, n))
    
            data = np.loadtxt('../sandbox/needs_frac_log2_mp%d_%s_%d.csv' % (length, func, n), delimiter=',')

            q, terms = data[:, 0], data[:, 2]

            terms = np.pad(terms, (medianwidth, medianwidth), mode='edge')
            terms_sliced = []
            for i in range(0, medianwidth * 2 + 1):
                terms_sliced.append(terms[i:][:len(q)])

            terms = np.sort(np.stack(terms_sliced), axis=0)[medianwidth]

            color = (1, j * 0.06, j * 0.06) if func == 'A' else (j * 0.06, j * 0.06, 1)

            if (j % 2) == 0:
                plt.plot(q, terms, label='$%s_{%d}$' % (func.lower(), n), linewidth='1', color=color, zorder=-n)
            else:
                plt.plot(q, terms, label='$%s_{%d}$' % (func.lower(), n), linewidth='1', color=color, zorder=-n, linestyle='dashed')

    qs, termss = [], []

    for func in ['A', 'B']:
        for j, n in enumerate(list(range(0, 32 + 1)) + [64, 128, 256]):
            if func == 'B' and n == 0:
                continue

            print('plot %s %d' % (func, n))
    
            data = np.loadtxt('../sandbox/needs_frac_log2_mp%d_%s_%d.csv' % (length, func, n), delimiter=',')

            q, terms = data[:, 0], data[:, 2]

            qs.append(q)
            termss.append(terms)

    qs = np.concatenate(qs)
    termss = np.concatenate(termss).astype(float)
    qs, termss = qs[qs < 8000], termss[qs < 8000] 
    
    popt, _ = curve_fit(ffunc, qs, termss, p0=(0, 0, 0))
    a, b, c = popt

    for bias in np.arange(1001) / 100:
        approx = ffunc(qs, a+bias, b, c)
        if np.all(approx + 0.5 > termss):
            break

    qs = 10**(np.arange(-350, 701) / 100)
    approx = ffunc(qs, a+bias, b, c)

    plt.plot(qs, approx, 
             linewidth='1', color='black', zorder=1, linestyle='dashed')
    plt.text(qs[-1], approx[-1], '$pow2(%.4lf + %.4e~log_2(q) + %.4e~log_2(q)^2)$' % (a+bias, b, c), horizontalalignment='right')
                        
    plt.legend(loc='upper left', ncol=2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xticks([1e-4, 1e-3, 1e-2, 1e-1, 1, 10, 100, 1e+3, 1e+4, 1e+5, 1e+6, 1e+7])
    plt.xlabel('q')
    plt.ylabel('c.frac terms')

    plt.savefig('../sandbox/needs_frac_log2_mp%d.png' % length)
    plt.savefig('../figures/needs_frac_log2_mp%d.svg' % length, bbox_inches='tight', pad_inches=0.1)