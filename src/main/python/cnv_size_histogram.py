import os
import pandas as pd
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys

cnv_file = sys.argv[1]
output_file = sys.argv[2]

print "Plotting size of cnvs from %s" % str(cnv_file)

all_truth = pd.read_table(cnv_file)
all_truth['span'] = np.log10(all_truth.end - all_truth.start)

# Little bit of a hack: not sure why targeted bp is null for some cases?
all_truth = all_truth[all_truth.tgbp != "null"]
all_truth['tgbp'] = np.log10(pd.to_numeric(all_truth.tgbp))


fig = plt.figure(figsize=(4,7))

plt.subplot(311)
plt.tick_params(axis='both', which='major', labelsize=7)
dp = sns.distplot(all_truth.span, kde=False, bins=10)
dp.set_title("Distribution of CNV Sizes")
dp.set_xlabel(r"Spanned Base Pairs ($log_{10}(bp)$)")

plt.subplot(312)
plt.tick_params(axis='both', which='major', labelsize=7)
dp = sns.distplot(all_truth.tgbp, kde=False, bins=20)
# dp.set_title("Distribution of CNV Sizes")
dp.set_xlabel(r"Targeted Base Pairs ($log_{10}(bp)$)")

plt.subplot(313)
plt.tick_params(axis='both', which='major', labelsize=7)
dp = sns.distplot(all_truth.tgbp, kde=False)
# dp.set_title("Distribution of CNV Sizes")
dp.set_xlabel(r"Target Regions")

fig.subplots_adjust(hspace=0.1)

plt.savefig(output_file)
