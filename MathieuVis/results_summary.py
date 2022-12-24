import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.clf()
plt.figure(figsize=(12, 6))

for n in range(0, 4 + 1):
    for func in ['B', 'A']:
        if func == 'B' and n == 0:
            continue

        print('plot %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % (func, n), delimiter=',')

        q, y = data['q'], data['convergence'] + n * n

        plt.plot(q, y, label='${}_{}$'.format(func.lower(), str(n)))

plt.legend(loc='lower left', ncol=2)

plt.xlabel('q')
plt.xlim([0, 16])
plt.ylim([-10, 20])

plt.savefig('../figures/eigen_plot.svg', bbox_inches='tight', pad_inches=0.1)