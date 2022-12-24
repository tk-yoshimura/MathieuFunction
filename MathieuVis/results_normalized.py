import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.clf()
plt.figure(figsize=(12, 6))

data = pd.read_csv('../sandbox/eigen_A_0_approx_ver2.csv', delimiter=',')
q, y = data['q'].to_numpy().astype(float), data['convergence'].to_numpy().astype(float)        
plt.plot(q, y, label='$a_0$', color=(0.2, 0.2, 0.2), zorder = 0)

for func in ['A', 'B']:
    for i, n in enumerate(list(range(1, 4 + 1)) + [8, 16, 32, 64, 128, 256]):
        print('plot %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % (func, n), delimiter=',')

        q, y = data['q'].to_numpy().astype(float), data['convergence'].to_numpy().astype(float)

        q /= max(1, n * n)
        y /= max(1, n * n)
        
        color = (1, i * 0.08, i * 0.08) if func == 'A' else (i * 0.08, i * 0.08, 1)

        if (i % 2) == 0:
            plt.plot(q, y, label='$%s_{%d}$' % (func.lower(), n), color=color, zorder = -n)
        else:
            plt.plot(q, y, label='$%s_{%d}$' % (func.lower(), n), color=color, zorder = -n, linestyle='dashed')

plt.legend(loc='lower left', ncol=2)

plt.xlabel('$q / max(1, n^2)$')
plt.ylabel('$(a,b - n^2) / max(1, n^2)$')
plt.xlim([0, 12])
plt.ylim([-16, 2])

plt.savefig('../figures/eigen_plot_normalized.svg', bbox_inches='tight', pad_inches=0.1)