import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def limit_a(n, q, level):
    h = np.sqrt(q)

    s = 2 * n + 1    
    y = 2 * (s * h - q) - ((6 * n * n + 2 * n + 1) / 4)
    
    if level > 0:
        y -= (s * (3 + s * s)) / (128 * h)
    if level > 1:
        y -= (9 + s * s * (34 + s * s * 5)) / (4096 * h * h)
    if level > 2:
        y -= (s * (405 + s * s * (410 + s * s * 33))) / (131072 * h * h * h)
    if level > 3:
        y -= (486 + s * s * (2943 + s * s * (1260 + s * s * 63))) / (1048576 * h * h * h * h)
    if level > 4:
        y -= (s * (41607 + s * s * (69001 + s * s * (15617 + s * s * 527)))) / (33554432 * h * h * h * h * h)

    return y

#for n in range(14 + 1):
#    data = pd.read_csv('../results/eigen_limit_precision64_n%d.csv' % (n), delimiter=',', skiprows=2)
#
#    u, a = data['u'].to_numpy().astype(float), data['a'].to_numpy().astype(float)
#
#    a_limit = limit_a(n, np.sqrt(u) * max(1, n * n))
#            
#    plt.clf()
#    plt.figure(figsize=(12, 6))
#
#    plt.plot(u, -a, label='$n=%d$' % (n), color='black')
#    plt.plot(u, -a_limit, label='$n=%d (limit)$' % (n), color='black')
#        
#    plt.grid()
#
#    plt.xlabel('$u$')
#    plt.ylabel('$-a$')
#
#    plt.xscale('log')
#    plt.yscale('log')
#
#    plt.xlim([1024, 1.099511627776e12])
#
#    plt.legend(loc='upper right')
#
#    plt.savefig('../sandbox/eigen_limit_plot_a_n%d.png' % (n), bbox_inches='tight', pad_inches=0.1)

for n in range(16 + 1):
    data = pd.read_csv('../results/eigen_limit_r2_precision64_n%d.csv' % (n), delimiter=',', skiprows=2)

    invu, a = data['1/u'].to_numpy().astype(float)[1:], data['a'].to_numpy().astype(float)[1:]
    u = 1 / invu

    plt.clf()
    plt.figure(figsize=(12, 6))

    for level in range(5):
        a_limit = limit_a(n, np.sqrt(u) * max(1, n * n), level)

        a_diff = a_limit - a
    
        plt.plot(np.sqrt(np.sqrt(invu)), a_diff / np.max(a_diff), label='$level=%d$' % (level))
        
    plt.grid()

    plt.xlabel('$sqrt(1 / u)$')
    plt.ylabel('$\Delta a$')

    plt.xlim([0, np.sqrt(1 / 32)])
    #plt.xticks(np.arange(0, 10+1) / 10000)

    #plt.yscale('log')

    plt.legend(loc='upper right')

    plt.savefig('../sandbox/eigen_limitdiff_plot_a_n%d.png' % (n), bbox_inches='tight', pad_inches=0.1)