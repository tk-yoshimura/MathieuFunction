import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


for n in range(32 + 1):
    for func in ['m', 'd']:
        if func == 'd' and n == 0:
            continue

        data = pd.read_csv('../results/ddouble/eigen_ddouble_results_%s_n%d.csv' % (func, n), delimiter=',')

        u, error = data['x'].to_numpy().astype(float), data['relative_error'].to_numpy().astype(float)
        
        plt.clf()
        plt.figure(figsize=(12, 6))
        
        plt.plot(u, error)

        plt.grid()

        plt.xlabel('$u$')
        plt.ylabel('relative_error')

        #plt.xscale('log')

        plt.xlim([0, 1024])

        plt.savefig('../sandbox/eigen_ddouble_error_%s_n%d.png' % (func, n), bbox_inches='tight', pad_inches=0.1)