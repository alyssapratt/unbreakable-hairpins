import scipy.stats as stats
import sys,statistics

## MAIN ##
usage = "usage: " + sys.argv[0] + " <endospore tsv file> <output file>"
## EXAMPLE:
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

tsv = sys.argv[1]
output = sys.argv[2]
endospore_genera = []
endosporeless_genera = []
endospore_breakable = 0
endospore_unbreakable = 0
endosporeless_breakable = 0
endosporeless_unbreakable = 0
with open(tsv, 'r') as F:
    for line in F:
        if not line.startswith("Genus"):
            proportion_unbreakable = float(line.split('\t')[5])
            endospore_forming_str = line.split('\t')[6]
            number_unbreakable = int(line.split('\t')[3])
            number_breakable = int(line.split('\t')[4]) - number_unbreakable
            if endospore_forming_str.lower() == 'true':
                endospore_genera.append(proportion_unbreakable)
                endospore_breakable += number_breakable
                endospore_unbreakable += number_unbreakable
            elif endospore_forming_str.lower() == 'false':
                endosporeless_genera.append(proportion_unbreakable)
                endosporeless_breakable += number_breakable
                endosporeless_unbreakable += number_unbreakable
F.close()
mean_endospore_unbreakable = statistics.mean(endospore_genera)
mean_endosporeless_unbreakable = statistics.mean(endosporeless_genera)

tstatistic,tpvalue = stats.ttest_ind(endospore_genera,endosporeless_genera,equal_var=False,alternative='greater')
print("t-test p-value: " + str(tpvalue))
ks_statistic,ks_pvalue = stats.ks_2samp(endosporeless_genera,endospore_genera,alternative='greater')
print("ks-test p-value: " + str(ks_pvalue))
oddsratio,fisher_pvalue = stats.fisher_exact([[endospore_breakable,endospore_unbreakable], [endosporeless_breakable,endosporeless_unbreakable]])
print("Fisher contingency table p-value: " + str(fisher_pvalue))
