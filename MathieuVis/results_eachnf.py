import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for func in ['A', 'B']:
    for n in range(0, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('plot %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, approx, convergence = data['q'], data['approx'], data['convergence']

        plt.clf()
        plt.plot(q, approx, label='approx')
        plt.plot(q, convergence, label='convergence', linewidth='1')
        plt.legend(loc='upper right')

        plt.savefig('../sandbox/eigen_%s_%d_approx.png' % (func, n))

        plt.clf()
        plt.plot(q[:len(q)//4], approx[:len(approx)//4], label='approx')
        plt.plot(q[:len(q)//4], convergence[:len(convergence)//4], label='convergence', linewidth='1')
        plt.legend(loc='upper right')

        plt.savefig('../sandbox/eigen_%s_%d_approx_short.png' % (func, n))