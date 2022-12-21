import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

for func in ['A', 'B']:

    ns, max_errors = [], []

    for n in range(0, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('plot %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        approx, convergence = data['approx'], data['convergence']

        error = np.abs(approx - convergence)
        error = error[error==error]
        max_error = np.max(error)

        ns.append(n)
        max_errors.append(max_error)

    approx_res = pd.DataFrame([ns, max_errors], index=['n', 'max_error']).T
    approx_res.to_csv('../sandbox/approx_error_%s.csv' % (func), index=False)

    plt.clf()
    plt.plot(ns, max_errors)
    plt.xlabel('n')
    plt.xlabel('max_error')

    plt.savefig('../sandbox/approx_error_%s.png' % (func))
