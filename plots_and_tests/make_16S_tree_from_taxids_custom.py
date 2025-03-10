import sys,re, time
from ete3 import Tree, NodeStyle, TreeStyle, NCBITaxa, AttrFace, TextFace
from Bio import SeqIO, Entrez
import matplotlib.cm as cm
from urllib.request import urlopen
from urllib.error import HTTPError

def collapsed_leaf(node):
    if len(node2labels[node]) == 1:
       return True
    else:
       return False

## MAIN ##
usage = "usage: " + sys.argv[0] + " <lineages file> <T1000 file> <ids to taxids file> <output tree file> <taxonomy filter> <rank filter>"
## EXAMPLE: make_16S_tree_from_taxids_custom.py all_16S_lineages.txt T1000.fasta ids_to_taxids_all_16S.txt tree.pdf Bacilli family
if len(sys.argv) != 7:
    print(usage)
    sys.exit()

lineages_file = sys.argv[1]
T1000_file = sys.argv[2]
ids_to_taxids_file = sys.argv[3]
output_file = sys.argv[4]
name_filter = sys.argv[5]
rank_filter = sys.argv[6]

T1000s = set()
taxid_to_ids = {}

ncbi = NCBITaxa()
Entrez.email = "prattal@oregonstate.edu"


with open(ids_to_taxids_file, 'r') as F:
    for line in F:
        if not line.startswith("downloaded"):
            try:
                id,taxid = line.strip().split(":")
                if int(taxid) in taxid_to_ids.keys():
                    taxid_to_ids[int(taxid)].append(int(id))
                else:
                    taxid_to_ids[int(taxid)]=[int(id)]
            except ValueError:
                continue
    F.close()

for record in SeqIO.parse(T1000_file, 'fasta'):
    id = int(record.id.split('_')[0])
    seq = str(record.seq)
    T1000s.add(id)

rank_taxids = {}
for taxid in taxid_to_ids.keys():
    try:
        rank = ncbi.get_rank([taxid])[taxid].lower()
        lineage = ncbi.get_lineage(taxid)
        rank_names = list(ncbi.get_taxid_translator(lineage).values())
        if rank == rank_filter and name_filter in rank_names:
            if taxid in rank_taxids.keys():
                rank_taxids[taxid].extend(taxid_to_ids[taxid])
            else:
                rank_taxids[taxid] = taxid_to_ids[taxid]
        elif rank is not None:
            lineage_ranks = ncbi.get_rank(lineage)
            rank_found = 0
            if name_filter in rank_names:
                for lineage_taxid in lineage_ranks.keys():
                    if lineage_ranks[lineage_taxid].lower() == rank_filter:
                        rank_found = 1
                        if lineage_taxid in rank_taxids.keys():
                            rank_taxids[lineage_taxid].extend(taxid_to_ids[taxid])
                        else:
                            rank_taxids[lineage_taxid] = taxid_to_ids[taxid]
        else:
            print("rank for " + taxid + " is None")
    except KeyError:
        print("no rank found: " + str(taxid))


t = ncbi.get_topology(rank_taxids.keys())

total_samples_represented = 0
with open("cur_taxa_info.txt", 'w') as myfile:
    for node in t.traverse("postorder"):
        if node.name != '':
            n = NodeStyle()
            n["hz_line_width"] = 2
            n["vt_line_width"] = 2
            n["shape"] = "square"
            n["size"] = 1
            n["fgcolor"] = "black"
            if int(node.name) in rank_taxids.keys():
                num_ids_for_taxid = len(rank_taxids[int(node.name)])
                total_samples_represented += num_ids_for_taxid
                num_unbreakable_ids = 0
                for id in rank_taxids[int(node.name)]:
                    if id in T1000s:
                        num_unbreakable_ids += 1
                proportion_unbreakable = num_unbreakable_ids/num_ids_for_taxid
                if num_ids_for_taxid < 1:
                    print("fewer than 4 examples, taxid excluded: " + node.name + ", ids: " + str(rank_taxids[int(node.name)]))
                else:
                    b = round(proportion_unbreakable * 255.0)
                    n["bgcolor"] = '#%02x%02x%02x' % (255, 255-b, 255-b)
                    if node.is_leaf():
                        sciname_face = AttrFace("sci_name", fsize=9)
                        sciname = getattr(node, "sci_name", None)
                        node.add_face(sciname_face, column=0, position="branch-right")
                        node.add_face(TextFace("("+str(num_unbreakable_ids) +"/" +str(num_ids_for_taxid)+")", fsize=9), column=0)
                        myfile.write(sciname +"\t" + str(num_unbreakable_ids) + "\t" + str(num_ids_for_taxid) + "\t" + str(proportion_unbreakable)+ "\n")
            node.set_style(n)
myfile.close()

def my_layout(node):
    if getattr(node, "rank", None) == 'phylum':
        rank_face = AttrFace("sci_name", fsize=7)
        node.add_face(rank_face, column=0, position="branch-top")



ts = TreeStyle()
ts.mode = "c"
ts.arc_start = -360 # 0 degrees = 3 o'clock
ts.arc_span = 360
ts.layout_fn = my_layout
ts.show_leaf_name = False
ts.title.add_face(TextFace("n = " + str(total_samples_represented), fsize=20), column=0)

t.render(output_file, tree_style=ts, w=500, units="mm")
t.write(outfile="16S_tree.nw")
