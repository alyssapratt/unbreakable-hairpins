This repository contains code used in characterizing unbreakable hairpins: hairpins that retain structure despite repeated Altschul-Erickson dinucleotide shuffles. There is also a script to calculate the number of unique sequences generatable by dinucleotide shuffling.

Building Stem Loop Dataset:

example:

$ python processHairpins_v3.py allHairpins.txt allSegments.txt





Finding Unbreakable Hairpins:

$ python filterHairpins_T1000.py allHairpins.fasta




Plotting and Statistical Tests:


Instructions for running count_unique_paths.py (tested on Python 3.6+)

Requires Biopython and Numpy

Run:

python3 count_unique_paths.py [RNA sequence]


