import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.clf()
fig, ax1 = plt.subplots(figsize=(12, 6))
ax2 = ax1.twiny()

for i, n in enumerate(list(range(1, 4 + 1)) + [8, 16, 32, 64]):
    data_a = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % ('A', n), delimiter=',')
    data_b = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % ('B', n), delimiter=',')

    q, a, b = data_a['q'].to_numpy().astype(float), data_a['convergence'].to_numpy().astype(float), data_b['convergence'].to_numpy().astype(float)

    y = (a + b)/2

    q /= max(1, n * n)
    y /= max(1, n * n)
        
    color = (i * 0.08, i * 0.08, i * 0.08)

    if (i % 2) == 0:
        ax1.plot(q*q, y, label='$n=%d$' % (n), color=color, zorder = -n)
    else:
        ax1.plot(q*q, y, label='$n=%d$' % (n), color=color, zorder = -n, linestyle='dashed')        
        
ax2.set_xticks(np.arange(0, 16 + 1), [])
ax2.grid()

ax1.set_xlabel('$q / n^2~~~(square~scale)$')
ax1.set_ylabel('$((a+b)/2 - n^2) / n^2$')
ax1.set_xticks(np.arange(0, 4 + 1)**2, np.arange(0, 4 + 1))

ax1.set_xlim([0, 16])
ax1.set_ylim([-2, 0.5])


ax1.legend(loc='lower left')

plt.savefig('../figures/eigen_plot_mean.svg', bbox_inches='tight', pad_inches=0.1)