
import os 
os.chdir(r'C:\Users\Teresa\OneDrive - University Of Cambridge\St Johns Cambridge\1-Part II Biochemistry\Part II Project\Files\Randomer sequencing')

import sys
from Bio.Seq import Seq
# def to read fasta files
def readFasta(fastaFile):
    fh = open(fastaFile, 'r')
    for line in fh:
        if line[0] == '>':
            header = line.rstrip()[1:]
            if sys.version_info[0] < 3:
                seq = fh.next().rstrip()
            else:
                seq = fh.readline().rstrip()
        yield [header, seq]
    fh.close()

dict_rounds = {}

rounds_names = ['RZ5N','F201N', 'F205N', 'N17L', 'N73L', 'N73H', 'FN17L', 'FN73L','FN73H'] 

count = {}
for i in rounds_names:
    dict_rounds[i]={}
    count[i]=0
    for line in readFasta(i+'.fasta'): # looks through each line of each file
        count[i]+=1
        if  line[1] not in dict_rounds[i]: #adds to dict if seq not there
            dict_rounds[i][line[1]] = 1
            
        else:
            dict_rounds[i][line[1]] +=1            

pool_s = ['RZ5N','F201N', 'F205N', 'N17L', 'N73L', 'N73H', 'FN17L', 'FN73L','FN73H'] 

dict_l = {}
for i in pool_s:
    dict_l[i]={}
    for seq in dict_rounds[i]:
        if len(seq) in dict_l[i]:
            dict_l[i][len(seq)]+=dict_rounds[i][seq]
        else:
            dict_l[i][len(seq)]=dict_rounds[i][seq]

#saving counts of length with number of length corresponding to count
import pandas as pd
dict_l_c = {}
dict_l_c['Pool']=[]
dict_l_c['Length']=[]
for i in dict_l:
    for l in dict_l[i]:
        for j in range(dict_l[i][l]):
            dict_l_c['Pool'].append(i)
            dict_l_c['Length'].append(l)
        
df = pd.DataFrame.from_dict(dict_l_c)
df.to_csv('Length_allselec_unfiltered.csv')