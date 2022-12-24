import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.clf()
fig = plt.figure(figsize=(12, 6))

for i, n in enumerate(list(range(1, 4 + 1)) + [8, 16, 32, 64]):
    data_a = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % ('A', n), delimiter=',')
    data_b = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % ('B', n), delimiter=',')

    q, a, b = data_a['q'].to_numpy().astype(float), data_a['convergence'].to_numpy().astype(float), data_b['convergence'].to_numpy().astype(float)

    y = (a + b)/2

    q /= max(1, n * n)
    y /= max(1, n * n)
        
    color = (i * 0.08, i * 0.08, i * 0.08)

    if (i % 2) == 0:
        plt.plot(q, y, label='$n=%d$' % (n), color=color, zorder = -n)
    else:
        plt.plot(q, y, label='$n=%d$' % (n), color=color, zorder = -n, linestyle='dashed')        

plt.xlabel('$q / n^2$')
plt.ylabel('$((a+b)/2 - n^2) / n^2$')

plt.xlim([0, 4])
plt.ylim([-2, 0.5])

plt.legend(loc='lower left')

plt.savefig('../figures/eigen_plot_mean_raw.svg', bbox_inches='tight', pad_inches=0.1)