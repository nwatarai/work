import numpy as np
import pandas as pd
import sys
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from matplotlib import cm

def mean_multiple_env(env, number_of_env, index, nmr_of_sample):
	out = []
	if nmr_of_sample * number_of_env > 500:
		quit("nmr_of_sample * number_of_env too many")
	for i in range(nmr_of_sample):
		multi = env[index[i*number_of_env:(i+1)*number_of_env],:]
		out.append(np.mean(multi, axis=0))
	return np.array(out)

def fetch_multi(before, after, number_of_env, nmr_of_sample):
	index = np.arange(500)
	np.random.shuffle(index)
	multibefore = mean_multiple_env(before, number_of_env, index, nmr_of_sample)
	multiafter = mean_multiple_env(after, number_of_env, index, nmr_of_sample)
	return multibefore, multiafter

def mean_std_correlation(before, after, nmr_of_sample):
	corr = []
	for i in range(nmr_of_sample):
		pearson = pearsonr(before[i],after[i])[0]
		corr.append(1.0 - pearson)
	return np.mean(corr), np.std(corr)

matrix = pd.read_csv(sys.argv[1], sep=",", header=0, index_col=0)

nmr_of_sample = 5
before = matrix.ix[:500,:].values
after = matrix.ix[500:1000,:].values

for i in range(1,101,3):
	nmr_of_sample = 500 / i
	b, a = fetch_multi(before, after, i, nmr_of_sample)
	mean, std = mean_std_correlation(b, a, nmr_of_sample)
	plt.errorbar(i, mean, yerr=std, marker="*", color="k")
plt.xlabel("# environments", fontsize=15, fontname='serif')
plt.ylabel("1 - Pearson (before and after purterbation)", fontsize=15, fontname='serif')
plt.ylim([0.04,0.16])

plt.show()