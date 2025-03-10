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

def read_bulge_file(structure_file):
    pk_gcContent = []
    pkf_gcContent = []
    with open(structure_file, 'r') as F:
        for line in F:
            if 'rna_ID' not in line:
                rna_id,sub_id,seq, start_pos,stop_pos, closing_pair_5p,closing_pos1_5p,closing_pos2_5p, closing_pair_3p, closing_pos1_3p, closing_pos2_3p, length, ispk=line.strip().split('\t')
                seq = seq.upper()
                #if  re.search('[NRYKMBDHV]', seq) is None and int(length) != 0:
                if re.search('[NRYKMBDHV]', seq) is None and int(length) > 1:
                    gcContent = (seq.count("G") + seq.count("C") + seq.count("S"))/float(length)
                    if ispk == "1":
                        pk_gcContent.append(gcContent)
                    else:
                        pkf_gcContent.append(gcContent)
    return(pk_gcContent,pkf_gcContent)

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





