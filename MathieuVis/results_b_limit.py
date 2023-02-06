import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def limit_b(n, q):
    s = 2 * n - 1    
    y = 2 * (s * np.sqrt(q) - q) - ((6 * n * n - 2 * n + 1) / 4)

    return y

for n in range(1, 16 + 1):
    data = pd.read_csv('../results/eigen_limit_precision64_n%d.csv' % (n), delimiter=',', skiprows=2)

    u, b = data['u'].to_numpy().astype(float), data['b'].to_numpy().astype(float)

    b_limit = limit_b(n, np.sqrt(u) * max(1, n * n))
            
    plt.clf()
    plt.figure(figsize=(12, 6))

    plt.plot(u, -b, label='$n=%d$' % (n), color='black')
    plt.plot(u, -b_limit, label='$n=%d (limit)$' % (n), color='black')
        
    plt.grid()

    plt.xlabel('$u$')
    plt.ylabel('$-b$')

    plt.xscale('log')
    plt.yscale('log')

    plt.xlim([1024, 1.099511627776e12])

    plt.legend(loc='upper right')

    plt.savefig('../sandbox/eigen_limit_plot_b_n%d.png' % (n), bbox_inches='tight', pad_inches=0.1)

for n in range(1, 16 + 1):
    data = pd.read_csv('../results/eigen_limit_precision64_n%d.csv' % (n), delimiter=',', skiprows=2)

    u, b = data['u'].to_numpy().astype(float), data['b'].to_numpy().astype(float)

    b_limit = limit_b(n, np.sqrt(u) * max(1, n * n))

    b_diff = b_limit - b

    plt.clf()
    plt.figure(figsize=(12, 6))

    plt.plot(1 / u, b_diff**4, label='$n=%d$' % (n), color='black')
        
    plt.grid()

    plt.xlabel('$1 / u$')
    plt.ylabel('$\Delta b^4$')

    plt.xlim([0, 0.001])
    plt.xticks(np.arange(0, 10+1) / 10000)

    plt.legend(loc='upper right')

    plt.savefig('../sandbox/eigen_limitdiff_plot_b_n%d.png' % (n), bbox_inches='tight', pad_inches=0.1)