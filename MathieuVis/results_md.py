import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.clf()
plt.figure(figsize=(12, 6))

for i, n in enumerate([1, 2, 3, 4, 5, 6, 7, 8, 16, 32]):
    data = pd.read_csv('../results/eigen_precision40_n%d.csv' % (n), delimiter=',', skiprows=2)

    u, m = data['u'].to_numpy().astype(float), data['mean'].to_numpy().astype(float)

    m /= n * n
        
    color = (i * 0.08, i * 0.08, i * 0.08)

    if (i % 2) == 0:
        plt.plot(u, m, label='$n=%d$' % (n), color=color, zorder = -n)
    else:
        plt.plot(u, m, label='$n=%d$' % (n), color=color, zorder = -n, linestyle='dashed')        
        
plt.grid()

plt.xlabel('$u$')
plt.ylabel('$m$')

plt.xlim([0, 16])
plt.xticks(np.arange(1, 16+1))
plt.ylim([-2, 0.5])
plt.yticks(np.arange(-8, 2+1) / 4)

plt.legend(loc='lower left')

plt.savefig('../figures/eigen_plot_m.svg', bbox_inches='tight', pad_inches=0.1)

plt.clf()
plt.figure(figsize=(12, 6))

for i, n in enumerate([1, 2, 3, 4, 5, 6, 7, 8]):
    data = pd.read_csv('../results/eigen_precision40_n%d.csv' % (n), delimiter=',', skiprows=2)

    u, d = data['u'].to_numpy().astype(float), data['scaled_diff'].to_numpy().astype(float)
    
    d = 1 / d - 1
        
    color = (i * 0.08, i * 0.08, i * 0.08)

    if (i % 2) == 0:
        plt.plot(u, d, label='$n=%d$' % (n), color=color, zorder = -n)
    else:
        plt.plot(u, d, label='$n=%d$' % (n), color=color, zorder = -n, linestyle='dashed')        
        
plt.grid()

plt.xlabel('$u$')
plt.ylabel('$d$')

plt.xlim([0, 16])
plt.xticks(np.arange(1, 16+1))
plt.ylim([0, 32])
plt.yticks(np.arange(1, 8+1) * 4)

plt.legend(loc='upper right')

plt.savefig('../figures/eigen_plot_d.svg', bbox_inches='tight', pad_inches=0.1)