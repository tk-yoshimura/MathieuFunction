import numpy as np
import pandas as pd

for n in range(16 + 1):
    data = pd.read_csv('../results/eigen_precision40_n%d.csv' % (n), delimiter=',', skiprows=2, dtype=str)

    a, b = data['a'][1:], data['b'][1:]
            
    a.to_csv('../sandbox/eigen_precision40_a_n%d.csv' % (n), index=False, header=False, line_terminator='",\n"')
    b.to_csv('../sandbox/eigen_precision40_b_n%d.csv' % (n), index=False, header=False, line_terminator='",\n"')