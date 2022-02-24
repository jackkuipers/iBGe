

from causaldag import unknown_target_igsp
import causaldag as cd
import random
from causaldag.utils.ci_tests import gauss_ci_suffstat, gauss_ci_test, MemoizedCI_Tester
from causaldag.utils.invariance_tests import gauss_invariance_suffstat, gauss_invariance_test, MemoizedInvarianceTester
import pandas
import numpy as np
import time

#kk = 2
#seed_number = 101

scaleN = 1

n = 100 # number of nodes
N = 4*n*scaleN # number of observations

exp_ints = scaleN*[5, 10, 20][kk - 1] # number of observations per intervention
exp_parents = 2 # expected number of parents

ni = 10 # number of possible interventions


# store values for later dataframe
setup_vec = [n, N, exp_ints, exp_parents, ni, seed_number]
names_setup_vec = ["n", "N", "ints", "parents", "ni", "seed"]


mystring = [""] * len(setup_vec)
for ii in range(len(setup_vec)):
    mystring[ii] = "_".join([names_setup_vec[ii], str(setup_vec[ii])])

f_name = "_".join(mystring)
subdir_name = "_".join(mystring[0:5])
#dir_name = "/".join(["./simresults_s", subdir_name])
dir_name = "/".join(["./simresults_g", subdir_name])


Ts = pandas.read_csv("".join([dir_name, "/", f_name, "_Ts.csv"]))["x"]
data = pandas.read_csv("".join([dir_name, "/", f_name, "_data.csv"]))


obs_samples = np.array(data.iloc[np.where(Ts == 0)])


diffTs = np.array(list(set(Ts)))
diffTs = diffTs[np.where(diffTs > 0)]

iv_samples_list = [np.array(data.iloc[np.where(Ts == ii)]) for ii in diffTs]


# Form sufficient statistics
obs_suffstat = gauss_ci_suffstat(obs_samples)
invariance_suffstat = gauss_invariance_suffstat(obs_samples, iv_samples_list)
setting_list = [dict(known_interventions=[]) for ii in diffTs]


SP_alphas = np.array([2.48e-08, 1.11e-07, 4.98e-07, 2.23e-06, 1.00e-05, 1.65e-05, 2.72e-05, 4.48e-05,7.39e-05])
DAGstring = [""] * len(SP_alphas)

print(seed_number)

start_time = time.time()
# Run UT-IGSP
# Create conditional independence tester and invariance tester
for ii in range(len(SP_alphas)):
    alpha = SP_alphas[ii]
    print(alpha)
    alpha_inv = alpha
    ci_tester = MemoizedCI_Tester(gauss_ci_test, obs_suffstat, alpha=alpha)
    invariance_tester = MemoizedInvarianceTester(gauss_invariance_test, invariance_suffstat, alpha=alpha_inv)
    try:
        est_dag, est_targets_list = unknown_target_igsp(setting_list, set(range(n)), ci_tester, invariance_tester)
    except:
        est_dag = "error"
        print("SVD error")
    DAGstring[ii] = est_dag

time_taken = time.time() - start_time
print(time_taken)

DAGs = pandas.DataFrame(DAGstring)

DAGs.to_csv("".join([dir_name, "/", f_name, "_UT-IGSP_DAGs.csv"]), header=False, index=False)

