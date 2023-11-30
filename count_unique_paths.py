from Bio import SeqIO
import numpy as np
import math
import sys

DEBUG = False # set to True to see contents of matrix representations

#####################################################################
# 2023 Alyssa Pratt, David Hendrix                                  #
# http://hendrixlab.cgrb.oregonstate.edu                            #
# Counts the number of unique sequences generatable by Altschul-    #
# Erickson dinucleotide shuffling and prints the result             #
# ----------------------------------------------------------------- #              
#####################################################################

usage = "usage: " + sys.argv[0] + " <RNA sequence>"
if len(sys.argv) != 2:
    print(usage)
    sys.exit()

seq = sys.argv[1]

#############
#Subroutines#
#############

def matrix_representation(sequence):
    rows,cols = (4,4)
    matrix = [[0]*cols for r in range(rows)]
    indegrees = [0,0,0,0]
    multiplicities=[[0]*cols for r in range(rows)]
    if DEBUG:
        for col in matrix:
            print("Row:" + str(col))
    nums = RNA_to_nums(sequence)
    indegrees[nums[0]]=1 # give false edge to first element
    multiplicities[nums[len(nums)-1]][nums[0]] = 1 # false edge between last and first element
    for i in range(0, len(nums) - 1):
        cur_node = nums[i]
        next_node = nums[i+1]
        indegrees[next_node]+=1
        multiplicities[cur_node][next_node]+=1
        if cur_node != next_node:
            matrix[cur_node][next_node] -= 1
            matrix[next_node][next_node]+=1
    root_index = nums[0]
    if DEBUG:
        for i in range (0,4):
            print("Row " + str(i) + ":" + str(matrix[i]))
    del matrix[root_index]
    for row in matrix:
        del row[root_index]
        if DEBUG:
            print("Row: " + str(col))
    return multiplicities,indegrees,matrix

def RNA_to_nums(sequence):
    new_seq_list = []
    for char in sequence:
        if (char == 'A'):
            new_seq_list.append(0)
        elif char == 'C':
            new_seq_list.append(1)
        elif char == 'G':
            new_seq_list.append(2)
        elif char == 'U':
            new_seq_list.append(3)
    return new_seq_list

def get_determinant(matrix_representation):
    numpy_matrix = np.asarray(matrix_representation)
    return np.linalg.det(numpy_matrix)

def get_numerator(indegrees):
    product = 1
    for vertex in indegrees:
        if vertex > 0:
            factorial = math.factorial(vertex-1)
            product *= factorial
    return product

def get_denominator(multiplicities):
    product = 1
    for i in range (len(multiplicities)):
        for j in range (len(multiplicities[i])):
            factorial = math.factorial(multiplicities[i][j])
            product *= factorial
    return product

######
#Main#
######

multiplicities,indegrees,matrix = matrix_representation(seq)
numerator = get_numerator(indegrees)
denominator = get_denominator(multiplicities)
det = get_determinant(matrix)
upaths = (det*numerator)/denominator
print(upaths)
