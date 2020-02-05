import scipy.io as sio
from scipy.stats import sem
import numpy as np
import sys

auprs = sio.loadmat(sys.argv[1])
aupr_mat = auprs['all_aups']
means = np.mean(aupr_mat, axis=1)
std_err = sem(aupr_mat, axis=1)
print(means)
print(std_err)
