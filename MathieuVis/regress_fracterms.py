import numpy as np
import pandas as pd
import scipy.stats as scstat
import matplotlib.pyplot as plt

medianwidth = 1

for length in [4, 8, 16]:
    plt.clf()
    plt.figure(figsize=(12, 6))

    for func in ['A', 'B']:
        for j, n in enumerate([0, 1, 2, 3, 4, 5, 6, 7, 8, 16, 32, 64, 128, 256]):
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
                plt.plot(q, terms, label='$%s_{%d}$' % (func.lower(), n), linewidth='1', color=color)
            else:
                plt.plot(q, terms, label='$%s_{%d}$' % (func.lower(), n), linewidth='1', color=color, linestyle='dashed')
                        
    plt.legend(loc='lower right', ncol=2)

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('q')
    plt.ylabel('c.frac terms')

    plt.savefig('../sandbox/needs_frac_log2_mp%d.png' % length)
    plt.savefig('../figures/needs_frac_log2_mp%d.svg' % length, bbox_inches='tight', pad_inches=0.1)