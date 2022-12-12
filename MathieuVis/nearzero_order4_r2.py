import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import itertools

def approx_order3(n, q):
    h = q * q

    n_sq = n * n
    n_sq_m1 = n_sq - 1
    n_sq_m4 = n_sq - 4
    n_sq_m9 = n_sq - 9

    c2 = 1 / (2 * n_sq_m1)
    c4 = (7 + n_sq * 5) / (32 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4)
    c6 = (29 + n_sq * (58 + n_sq * 9)) / (64 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m1 * n_sq_m4 * n_sq_m9)

    y = h * (c2 + h * (c4 + h * c6))

    return y

def approx_term4(n, q, rn, r0, r1, r2, r3):
    h = q * q

    n_sq = n * n
    n_sq_m1 = n_sq - 1
    n_sq_m4 = n_sq - 4
    n_sq_m9 = n_sq - 9
    n_sq_m16 = n_sq - 16

    c8 = (r0 + n_sq * (r1 + n_sq * (r2 + n_sq * r3))) / (rn * (n_sq_m1**7) * n_sq_m4 * n_sq_m9 * n_sq_m16)

    if not np.all(np.isfinite(c8)):
        return np.nan

    y = h**4 * c8[:, np.newaxis]

    return y

for func in ['A', 'B']:
    qs, ns, residuals = [], [], []

    for n in range(5, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../results/eigen_%s_%d_approx.csv' % (func, n), delimiter=',')

        q, expected = data['q'], data['convergence']
        q, expected = q[len(q)//16:len(q)//8], expected[len(expected)//16:len(expected)//8]

        approx = approx_order3(n, q) 

        qs.append(q)
        ns.append(n)
        residuals.append(approx - expected)    

    qs = np.stack(qs)
    ns = np.stack(ns).astype(float)
    residuals = np.stack(residuals)

    with open('../results/eigen_%s_term4.csv' % (func), 'w') as f:

        for rn in 2**np.arange(3, 16):
            min_error = float('inf')
            min_error_r0, min_error_r1, min_error_r2, min_error_r3 = 0, 0, 0, 0
            
            r0_min, r0_max = -16384, 16384
            r1_min, r1_max = -16384, 16384
            r2_min, r2_max = -16384, 16384
            r3_min, r3_max = -16384, 16384
            
            for r3 in range(r3_min, r3_max+1):
                actuals = approx_term4(ns, qs, rn, min_error_r0, min_error_r1, min_error_r2, r3)

                if np.any(np.isnan(actuals)):
                    continue

                error = residuals - actuals
                                                                                
                sum_abserror = np.sum(np.abs(error) / (ns * ns)[:, np.newaxis])

                if not min_error > sum_abserror:
                    continue

                min_error = sum_abserror

                print(min_error)
                print('{},{},{},{},{}'.format(rn, min_error_r0, min_error_r1, min_error_r2, r3))

                min_error_r3 = r3

            r3_min, r3_max = min_error_r3-8, min_error_r3+8

            for r3 in range(r3_min, r3_max+1):
                for r2 in range(r2_min, r2_max+1):
                    actuals = approx_term4(ns, qs, rn, min_error_r0, min_error_r1, r2, r3)

                    if np.any(np.isnan(actuals)):
                        continue

                    error = residuals - actuals
                                                                                
                    sum_abserror = np.sum(np.abs(error))

                    if not min_error > sum_abserror:
                        continue

                    min_error = sum_abserror

                    print(min_error)
                    print('{},{},{},{},{}'.format(rn, min_error_r0, min_error_r1, r2, r3))

                    min_error_r3 = r3
                    min_error_r2 = r2

            r3_min, r3_max = min_error_r3-8, min_error_r3+8
            r2_min, r2_max = min_error_r2-8, min_error_r2+8

            for r3 in range(r3_min, r3_max+1):
                for r2 in range(r2_min, r2_max+1):
                    for r1 in range(r1_min, r1_max+1):
                        actuals = approx_term4(ns, qs, rn, min_error_r0, r1, r2, r3)

                        if np.any(np.isnan(actuals)):
                            continue

                        error = residuals - actuals
                                                                                
                        sum_abserror = np.sum(np.abs(error))

                        if not min_error > sum_abserror:
                            continue

                        min_error = sum_abserror

                        print(min_error)
                        print('{},{},{},{},{}'.format(rn, min_error_r0, r1, r2, r3))

                        min_error_r3 = r3
                        min_error_r2 = r2
                        min_error_r1 = r1

            r3_min, r3_max = min_error_r3-8, min_error_r3+8
            r2_min, r2_max = min_error_r2-8, min_error_r2+8
            r1_min, r1_max = min_error_r1-8, min_error_r1+8

            for r3 in range(r3_min, r3_max+1):
                for r2 in range(r2_min, r2_max+1):
                    for r1 in range(r1_min, r1_max+1):
                        for r0 in range(r0_min, r0_max+1):
                            actuals = approx_term4(ns, qs, rn, r0, r1, r2, r3)

                            if np.any(np.isnan(actuals)):
                                continue

                            error = residuals - actuals
                                                                                
                            sum_abserror = np.sum(np.abs(error))

                            if not min_error > sum_abserror:
                                continue

                            min_error = sum_abserror

                            print(min_error)
                            print('{},{},{},{},{}'.format(rn, r0, r1, r2, r3))

                            min_error_r3 = r3
                            min_error_r2 = r2
                            min_error_r1 = r1
                            min_error_r0 = r0

            f.write('{},{},{},{},{},{}\n'.format(min_error, rn, min_error_r0, min_error_r1, min_error_r2, min_error_r3))