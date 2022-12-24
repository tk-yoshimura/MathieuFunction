import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for func in ['A', 'B']:
    for n in range(0, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('plot %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % (func, n), delimiter=',')

        q, y = data['q'].to_numpy(), data['convergence'].to_numpy()

        curvature = np.empty(len(y))
        curvature[1:-1] = np.abs(y[:-2] + y[2:] - 2 * y[1:-1])
        curvature[0] = curvature[1]
        curvature[-1] = curvature[-2]

        plt.clf()
        plt.plot(q, curvature, label='curvature', linewidth='1')
        plt.legend(loc='upper right')

        plt.xlabel('q')
        plt.yscale('log')

        plt.savefig('../sandbox/eigen_%s_%d_curvature.png' % (func, n))