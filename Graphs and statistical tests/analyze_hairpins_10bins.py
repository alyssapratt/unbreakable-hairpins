import sys
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
#############
#Subroutines#
#############

def read_hairpin_file(hairpin_file):
    hairpin_data = {}
    with open(hairpin_file, 'r') as F:
        for line in F:
            if 'rna_ID' not in line:                
                rna_id,sub_id,start_pos,stop_pos,seq,length,closing_pair,closing_pos1,closing_pos2,ispk=line.strip().split('\t')
                if "N" not in seq:
                    if rna_id not in hairpin_data:
                        hairpin_data[rna_id] = []
                    hairpin_data[rna_id].append((seq,length,ispk))
    return hairpin_data


def sort_plot_data(hairpin_data):
    pkf_gc = []
    pk_gc = []
    pkf_len = []
    pk_len = []
    for rna_id in hairpin_data:
        for hairpin in hairpin_data[rna_id]:
            seq,length,ispk = hairpin
            if int(length) != 0:
                gccount = seq.count("G") + seq.count("C")
                gcpercent = (float(gccount) / float(length))*100
                if int(ispk) == 1:
                    pk_gc.append(gcpercent)
                    pk_len.append(length)
                else:
                    pkf_gc.append(gcpercent)
                    pkf_len.append(length)
    
    plt.hist(pk_gc, normed=1, bins=range(0,101,10),cumulative=False,fill=False,histtype='step',label='hairpins with PKs')
    plt.hist(pkf_gc,normed=1,bins=range(0,101,10),cumulative=False,fill=False,histtype='step',label='hairpins without PKS')
    plt.xlabel('hairpin GC content')
    plt.ylabel('distribution')
    plt.legend(loc='upper right')
    plt.title('histogram of hairpin GC content')
    plt.savefig('../../plots/hairpins/hist_gccontent_10bins.pdf')


######
#Main#
######


usage = "usage: " + sys.argv[0] + " <hairpin file>"

hairpin_file = sys.argv[1]

hairpin_data = read_hairpin_file(hairpin_file)
sort_plot_data(hairpin_data)
print "done"
