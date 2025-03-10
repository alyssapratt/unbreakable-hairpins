import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import re
from scipy import stats

#############
#Subroutines#
#############

def read_id_list(structure_file):
    bprna_ids = []
    with open(structure_file, 'r') as F:
        for line in F:
            bprna_ids.append(line)
    return(bprna_ids)

def get_dbn_files()

######
#Main#
######

usage = "usage: " + sys.argv[0] + " <bulge file>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

bulge_file = sys.argv[1]

pk_gcContent,pkf_gcContent = read_bulge_file(bulge_file)

print "PK bulges found: ", len(pk_gcContent)
print "PKF bulges found: ", len(pkf_gcContent)
ks_results = stats.ks_2samp(pk_gcContent, pkf_gcContent)
pval = ks_results[1]
print "K-S p-value: " , pval
all_data = pk_gcContent + pkf_gcContent
binSize = 0.1
bins = list(np.arange(min(all_data),max(all_data),binSize)) + [max(all_data)]

plt.figure()
plt.xlabel("GC Content")
plt.ylabel("Density")
plt.hist(pk_gcContent, bins=bins, density=True, alpha=0.5, label="PK bulges")
plt.hist(pkf_gcContent, bins=bins, density=True, alpha=0.5, label="PKF bulges")
plt.annotate("KS test p-val: " + str(pval), xy=(-0.021,3.55))
plt.legend()
plt.savefig("gcContentHist_overOneBaseLong.pdf")
