import sys, re, math, random
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

#############
#Subroutines#
#############

def make_deletions(fasta_file, fasta_output):
    modified_sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        randpos = random.randrange(0,len(seq))
        mod_seq = seq[:randpos] + seq[randpos+1:]
        modified_sequences.append(SeqRecord(Seq(mod_seq), id=rna_id, description=''))
    SeqIO.write(modified_sequences,fasta_output,'fasta')

def make_insertions(fasta_file, fasta_output):
    modified_sequences = []
    for record in SeqIO.parse(fasta_file, 'fasta'):
        rna_id = record.id
        seq = str(record.seq)
        randpos = random.randrange(0,len(seq))
        randnucl = random.choice('UCGA')
        mod_seq = seq[:randpos] + randnucl+ seq[randpos:]
        modified_sequences.append(SeqRecord(Seq(mod_seq), id=rna_id, description=''))
    SeqIO.write(modified_sequences,fasta_output,'fasta')

######
#Main#
######

usage = "usage: " + sys.argv[0] + " <input fasta file> <output fasta file>"
if len(sys.argv) != 3:
    print(usage)
    sys.exit()

fasta_file = sys.argv[1]
output_file = sys.argv[2]

make_insertions(fasta_file, output_file)
