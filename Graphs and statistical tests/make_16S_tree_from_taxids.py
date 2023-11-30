import sys,re, time
from ete3 import Tree, NodeStyle, TreeStyle, NCBITaxa, AttrFace
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
usage = "usage: " + sys.argv[0] + " <lineages file> <fasta file> <T1000 file> <ids to taxids file> <output tree file>"
if len(sys.argv) != 6:
    print(usage)
    sys.exit()

lineages_file = sys.argv[1]
fasta_file = sys.argv[2]
T1000_file = sys.argv[3]
ids_to_taxids_file = sys.argv[4]
output_file = sys.argv[5]

sequence_dict = {}
T1000s = set()
T1000_nodes = set()
sample_taxids = set()
taxid_to_ids = {}
known_ids = set()

ncbi = NCBITaxa()
Entrez.email = "prattal@oregonstate.edu"


with open(ids_to_taxids_file, 'r') as F:
    for line in F:
        if not line.startswith("downloaded"):
            try:
                id,taxid = line.strip().split(":")
                known_ids.add(int(id))
                sample_taxids.add(taxid)
                if int(taxid) in taxid_to_ids.keys():
                    taxid_to_ids[int(taxid)].append(int(id))
                else:
                    taxid_to_ids[int(taxid)]=[int(id)]
            except ValueError:
                continue
    F.close()

for record in SeqIO.parse(fasta_file, 'fasta'):
    id = int(record.id.split("|")[0].split("_")[2])
    seq = str(record.seq)
    if seq in sequence_dict:
        sequence_dict[seq].append(id)
    else:
        sequence_dict[seq] = [id]

for record in SeqIO.parse(T1000_file, 'fasta'):
    id = int(record.id.split('_')[0])
    seq = str(record.seq)
    T1000s.add(id)

nodes_to_delete = set()
nodes_to_modify = set()
nodes_deleted = {}
genus_taxids = {}
for taxid in taxid_to_ids.keys():
    try:
        rank = ncbi.get_rank([taxid])[taxid].lower()
        lineage = ncbi.get_lineage(taxid)
        rank_names = list(ncbi.get_taxid_translator(lineage).values())
        if rank == 'class' and "Bacteria" in rank_names:
            if taxid in genus_taxids.keys():
                genus_taxids[taxid].extend(taxid_to_ids[taxid])
            else:
                genus_taxids[taxid] = taxid_to_ids[taxid]
        elif rank is not None:
            lineage_ranks = ncbi.get_rank(lineage)
            genus_found = 0
            if "Bacteria" in rank_names:
                for lineage_taxid in lineage_ranks.keys():
                    if lineage_ranks[lineage_taxid].lower() == 'class':
                        genus_found = 1
                        if lineage_taxid in genus_taxids.keys():
                            genus_taxids[lineage_taxid].extend(taxid_to_ids[taxid])
                        else:
                            genus_taxids[lineage_taxid] = taxid_to_ids[taxid]
            #if genus_found == 0:
            #    print("no valid rank in lineage: " + str(taxid))
        else:
            print("rank for " + taxid + " is None")
    except KeyError:
        print("no rank found: " + str(taxid))


nstyle = NodeStyle()
nstyle["shape"] = "square"
nstyle["fgcolor"] = "darkred"
nstyle["size"] = 10
nstyle["bgcolor"] = "white"

nstylet1000 = NodeStyle()
nstylet1000["fgcolor"] = "green"
nstylet1000["size"] = 10
nstylet1000["bgcolor"] = "LightSteelBlue"

nstylet1000_dup = NodeStyle()
nstylet1000_dup["shape"] = "square"
nstylet1000_dup["fgcolor"] = "darkred"
nstylet1000_dup["size"] = 10
nstylet1000_dup["bgcolor"] = "LightSteelBlue"



t = ncbi.get_topology(genus_taxids.keys())

#for node in t.traverse("postorder"):
#    if node.name != '':
#        if accession_to_id[node.name] in nodes_to_delete:
#            node.delete()
#        if accession_to_id[node.name] in nodes_to_modify:
#            node.set_style(nstyle)
#            if accession_to_id[node.name] in T1000s:
#                node.set_style(nstylet1000_dup)
#        elif accession_to_id[node.name] in T1000s:
#            node.set_style(nstylet1000)
#            T1000_nodes.add(node)

                #num_identical_nodes = len(nodes_deleted[node.name])
                #print(node.name + " has identical nodes")
            #print(node.name + " has identical nodes")

for node in t.traverse("postorder"):
    if node.name != '':
        n = NodeStyle()
        n["hz_line_width"] = 2
        n["vt_line_width"] = 2
        n["shape"] = "square"
        n["size"] = 1
        n["fgcolor"] = "black"
        if int(node.name) in genus_taxids.keys():
            num_ids_for_taxid = len(genus_taxids[int(node.name)])
            num_unbreakable_ids = 0
            for id in genus_taxids[int(node.name)]:
                if id in T1000s:
                    num_unbreakable_ids += 1
            proportion_unbreakable = num_unbreakable_ids/num_ids_for_taxid
            if proportion_unbreakable > 0.95:
                print(node.name + " ids: " + str(genus_taxids[int(node.name)]) + " " + str(proportion_unbreakable))
            b = round(proportion_unbreakable * 255.0)
            n["bgcolor"] = '#%02x%02x%02x' % (255, 255-b, 255-b)
        node.set_style(n)

def my_layout(node):
    if node.is_leaf():
        sciname_face = AttrFace("sci_name", fsize=9)
        node.add_face(sciname_face, column=0, position="branch-right")


ts = TreeStyle()
ts.mode = "c"
ts.arc_start = -360 # 0 degrees = 3 o'clock
ts.arc_span = 360
ts.layout_fn = my_layout
ts.show_leaf_name = False

#print("there are a total of " + str(len(nodes_to_delete)) +" deleted nodes")
#print("There are " + str(len(T1000s)) + " T1000s")

#ca = t.get_common_ancestor(T1000_nodes)
#ca.set_style(ca_style)

t.render(output_file, tree_style=ts, w=500, units="mm")
t.write(outfile="16S_tree.nw")
