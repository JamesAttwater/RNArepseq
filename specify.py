import random as ran                # random number generation
import numpy as np                  # import numpy
import sys                          # import sys for kill command
import Bio                          # import biopython
from Bio.Seq import Seq             # import Seq type from biopython
from Bio.Alphabet import IUPAC      # import IUPAC alphabet from biopython
from Bio import SeqIO               # import seq In/out from biopython
from Bio.SeqRecord import SeqRecord # import SeqRecord from biopython

print('Input sequence')
full = str(input())      
#full = 'GGAGACGGUCGGGUCCAGAUAUUCGUAUCUGUCGAGUAGAGUGUGGGCUCC'
d = len(full)
ntri = int(d/3)
print('file name')
fn = str(input())      

# Generate species

with open(fn + '.txt', 'w') as myfile:
    for j in range(ntri-1):
        for i in range(j+1):
            spec = full[i*3:(3*(ntri-j))+(i*3)]
            myfile.write('>1' + '\n' + spec + '\n')
myfile.close() 
