# script to count HTS reads of ribozyme extension reaction intemediates onto expected species and analyse the counts
#
# v1.0 original script - Count in based on expected species, allowing variable errors at 3' and 5' addition sites.
# Suggested rules are internal triplets 100% match, extremem triplets allow 1 missmatch each, except for hexamers where total missmatches <= 1.
# v1.2 add fractional tally counts for multiple mapped species

import numpy as np                      # import numpy
import csv                              # import csv reader / writer
import argparse
import operator
from functools import reduce
from Bio.Seq import Seq                 # import Seq type from biopython
from Bio.Alphabet import IUPAC          # import IUPAC alphabet from biopython
from Bio import SeqIO                   # import seq In/out from biopython
from Bio.SeqRecord import SeqRecord     # import SeqRecord from biopython

parser = argparse.ArgumentParser(description='intermaster1.0: Counts HTS reads of ribozyme extension reaction intemediates onto expected species and analyse the counts')
parser.add_argument('-i', help='Input text file detailing clipped HTS FASTA file names.', type=argparse.FileType('r'), required=True)
parser.add_argument('-s', help='File listing expected speices sequences as DNA', type=argparse.FileType('r'), required=True)
parser.add_argument('-d', help='Hamming distances allowed in the order proximal_triplet distal_triplet hexamer', type=int, nargs=3, required=True)
parser.add_argument('-csvout', help="Specify filename for .csv export")
args = parser.parse_args()

species = []
for specie in args.s:
    if specie[0] != '#':
        species.append(specie.split())
species = reduce(operator.add, species)
nspec = len(species)

filehandles = []
for filehandle in args.i:
    if filehandle[0] != '#':
        filehandles.append(filehandle.split())
filehandles = reduce(operator.add, filehandles)
nfile = len(filehandles)

#hamms = []
#for i in range(3):
#    hamms[i] = int(args.d[i])

def hamming(str1, str2):                                            #function to calculate hamming distance between two sequences
    return sum(map(str.__ne__, str1, str2))

def search_list(element, list_element):                             #function to search list and return index of hits
    try:
        index_element = list_element.index(element)
        return index_element
    except ValueError:
        return nspec

def align(subject, targets, dist5, dist3, hexdist, array, col):          #function to map reads onto targets and count into an array
    results = []
    for target in targets:
        if len(subject) == len(target):
            if len(subject) ==6:
                if hamming(subject, target) <= (hexdist):
                    results.append(target)
            else:    
                if hamming(subject[3:-3], target[3:-3]) == 0:
                    if hamming(subject[:3], target[:3]) <= dist5 and hamming(subject[-3:], target[-3:]) <= dist3:
                        results.append(target)
    if len(results) == 0:
        array[nspec,col] = array[nspec,col] + 1
    elif len(results) == 1:
        ind = search_list(results[0], targets)
        array[ind,col] = array[ind,col] + 1
    else:
        for result in results:
            ind = search_list(result, targets)
            array[ind,col] = array[ind,col] + (1/len(results))
    return() 

count = np.zeros(((nspec+1),(nfile+2)), dtype=object)
for i in range(nspec+1):
    count[i,0] = i
count[:nspec,1] = np.asarray(species)
count[nspec,1] = 'Unmatched'

i = -1
for file in filehandles:
    i = i+1
    for seq_record in SeqIO.parse((file + '.fasta'), "fasta"):
        read = str(seq_record.seq)
        align(read, species, args.d[0], args.d[1], args.d[2], count, (2+i))

totals = np.zeros((1,nfile), dtype = int)
norfacs = np.zeros((1,nfile), dtype = int)
#for k in range(nfile):
#    count[:nspec,(k+2)] = (count[:nspec,(k+2)] / count[0,(k+2)])*20           #into pmol
#    count[1:nspec,(k+2)] = (count[1:nspec,(k+2)] / np.sum(count[1:nspec,(k+2)]))           #as % of rxn molecules

print(count)

printfile = open('count_' + args.csvout + '.csv', 'w')
with printfile: 
    writer = csv.writer(printfile, dialect='excel')
    writer.writerows(count)