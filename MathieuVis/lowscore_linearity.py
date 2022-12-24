import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.clf()

for func in ['A', 'B']:

    ns, max_errors = [], []

    for n in range(0, 256 + 1):
        if func == 'B' and n == 0:
            continue

        print('read %s %d' % (func, n))

        data = pd.read_csv('../sandbox/eigen_%s_%d_approx_ver2.csv' % (func, n), delimiter=',')
        data = data[len(data)//24:]

        q, y, score = data['q'].to_numpy().astype(float), data['convergence'].to_numpy().astype(float), data['score'].to_numpy().astype(float)

        peak_index = np.argmax(y)

        lowscore_indexes = np.where(score < 0.999)[0]
        lowscore_indexes = lowscore_indexes[lowscore_indexes > 4]
        lowscore_indexes = lowscore_indexes[lowscore_indexes < len(score) - 4]
        #lowscore_indexes = lowscore_indexes[np.abs(lowscore_indexes - peak_index) >= 8]

        if len(lowscore_indexes) == 0:
            continue
        
        yn4, yn3, yn2, yn1 = y[lowscore_indexes-4], y[lowscore_indexes-3], y[lowscore_indexes-2], y[lowscore_indexes-1]
        y0 = y[lowscore_indexes]
        yp4, yp3, yp2, yp1 = y[lowscore_indexes+4], y[lowscore_indexes+3], y[lowscore_indexes+2], y[lowscore_indexes+1]

        lowscores = score[lowscore_indexes]

        linearities1 = np.abs((yn1 + yp1 - y0 * 2) / np.maximum(1e-30, np.abs(yn1 - yp1)))
        linearities2 = np.abs((yn2 + yp2 - y0 * 2) / np.maximum(1e-30, np.abs(yn2 - yp2)))
        linearities3 = np.abs((yn3 + yp3 - y0 * 2) / np.maximum(1e-30, np.abs(yn3 - yp3)))
        linearities4 = np.abs((yn4 + yp4 - y0 * 2) / np.maximum(1e-30, np.abs(yn4 - yp4)))

        linearities = np.minimum(np.minimum(linearities1, linearities2), np.minimum(linearities3, linearities4))
         
        plt.scatter(lowscores, linearities, s=np.full(len(lowscores), 1), color='gray', alpha=0.5)

plt.xlabel('score')
plt.ylabel('linearity')
plt.yscale('log')
plt.savefig('../sandbox/eigen_lowscores_withpeak.png')