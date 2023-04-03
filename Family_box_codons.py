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


rounds_names = ['RZ5N','F201N', 'F205N', 'N17L', 'N73L', 'N73H', 'FN17L', 'FN73L','FN73H'] #paper

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

dict_filt3_range_rounds = {} #only take sequences that are multiples of 3 within right range 
#except pools with PNK are not multiple of 3 but right range and above 8 nt
dict_filt3_rest={} #all sequences that are not multiple of 3

def is_integer_num(n):
    if isinstance(n, int):
        return True
    if isinstance(n, float):
        return n.is_integer()
    return False

count_range = {} #new counting dictionary to include different categories
length_count ={} #counts number of sequences in pool at length all in range and filtered 3


count_categories=['total','greater5len','filtered3_greater5len','in_range','in_range_filtered3', 
                  'fraction_greater5len_total',
                  'fraction_filtered3_greater5len_total',
                  'fraction_in_range_total',
                  'fraction_in_range_filtered3_total']

for x in rounds_names:
    count_range[x] = dict.fromkeys([count_categories[y] for y in range(len(count_categories))],0)
    
for i in rounds_names:
    dict_filt3_range_rounds[i]={}
    length_count[i]={}
    dict_filt3_rest[i]={}
    for seq in dict_rounds[i]:
        count_range[i]['total']+= dict_rounds[i][seq]
        
        if len(seq)>5: #sequences with length 6 or more
            count_range[i]['greater5len']+=dict_rounds[i][seq]
            
            if  seq not in dict_filt3_range_rounds[i]:
                if is_integer_num(len(seq)/3): #sequences with length multiple of 3
                    count_range[i]['filtered3_greater5len']+=dict_rounds[i][seq]
                    
                    if i == 'N17L' or i =='N73L' or i =='FN17L' or i =='FN73L': #range 12-30 nt plus minus 6 (below 6 already filtered before)
                        if len(seq) <36:

                            dict_filt3_range_rounds[i][seq]=dict_rounds[i][seq]
                            count_range[i]['in_range_filtered3'] += dict_rounds[i][seq]
                            
                            if len(seq) not in length_count[i]:
                                length_count[i][len(seq)]=dict_rounds[i][seq]
                            else:
                                length_count[i][len(seq)]+=dict_rounds[i][seq]

                    elif i == 'N73H' or i =='FN73H': #range 33-150 nt plus minus 6 
                        #--> adjusted to lower end as most reads seem to be shorter in pool
                        if len(seq) > 21 and len(seq) < 156:
                            dict_filt3_range_rounds[i][seq]=dict_rounds[i][seq] 
                            count_range[i]['in_range_filtered3'] += dict_rounds[i][seq]
                            
                            if len(seq) not in length_count[i]:
                                length_count[i][len(seq)]=dict_rounds[i][seq]
                            else:
                                length_count[i][len(seq)]+=dict_rounds[i][seq]

                    elif i == 'FN847':
                        if len(seq) > 15 and len(seq) < 33: #around 27 nt plus minus 6 
                            #adjusted to lower end as most reads here seem to be shorter (see plot)
                            dict_filt3_range_rounds[i][seq]=dict_rounds[i][seq]
                            count_range[i]['in_range_filtered3'] += dict_rounds[i][seq]
                            
                            if len(seq) not in length_count[i]:
                                length_count[i][len(seq)]=dict_rounds[i][seq]
                            else:
                                length_count[i][len(seq)]+=dict_rounds[i][seq]

                    elif i == 'RZ1N' or i == 'RZ5N' or i== 'DG1N' or i == 'DG5N' or i == 'F201N' or i =='F205N':
                        if len(seq) < 27: #6-21 nt plus minus 6
                            dict_filt3_range_rounds[i][seq]=dict_rounds[i][seq]
                            count_range[i]['in_range_filtered3'] += dict_rounds[i][seq]
                            
                            if len(seq) not in length_count[i]:
                                length_count[i][len(seq)]=dict_rounds[i][seq]
                            else:
                                length_count[i][len(seq)]+=dict_rounds[i][seq]

                    elif i == 'N145' or i == 'FN145':
                        if len(seq) < 51: #9-45 nt plus minus 6
                            dict_filt3_range_rounds[i][seq]=dict_rounds[i][seq]
                            count_range[i]['in_range_filtered3'] += dict_rounds[i][seq]
                            
                            if len(seq) not in length_count[i]:
                                length_count[i][len(seq)]=dict_rounds[i][seq]
                            else:
                                length_count[i][len(seq)]+=dict_rounds[i][seq]
                                
                    elif i == 'N73A' or i == 'FN73A': #range 12-150 nt plus minus 6
                        if len(seq) <156:
                            dict_filt3_range_rounds[i][seq]=dict_rounds[i][seq] 
                            count_range[i]['in_range_filtered3'] += dict_rounds[i][seq]
                            count_range[i]['in_range'] += dict_rounds[i][seq]
                            
                            if len(seq) not in length_count[i]:
                                length_count[i][len(seq)]=dict_rounds[i][seq]
                            else:
                                length_count[i][len(seq)]+=dict_rounds[i][seq]
                else: #takes additional sequences for these pools without requiring multiple of 3
                    if i == 'N73A' or i == 'FN73A': #range 12-150 nt plus minus 6
                        if len(seq) <156:
                            dict_filt3_range_rounds[i][seq]=dict_rounds[i][seq] 
                            count_range[i]['in_range'] += dict_rounds[i][seq]
                            
                            if len(seq) not in length_count[i]:
                                length_count[i][len(seq)]=dict_rounds[i][seq]
                            else:
                                length_count[i][len(seq)]+=dict_rounds[i][seq]
                    if not is_integer_num(len(seq)/3):
                        
                        if i == 'N17L' or i =='N73L' or i =='FN17L' or i =='FN73L': #range 12-30 nt plus minus 6 (below 6 already filtered before)
                            if len(seq) <36:

                                dict_filt3_rest[i][seq]=dict_rounds[i][seq]

                                

                        elif i == 'N73H' or i =='FN73H': #range 33-150 nt plus minus 6 
                            #--> adjusted to lower end as most reads seem to be shorter in pool
                            if len(seq) > 21 and len(seq) < 156:
                                dict_filt3_rest[i][seq]=dict_rounds[i][seq] 

                                
                        elif i == 'RZ1N' or i == 'RZ5N' or i== 'DG1N' or i == 'DG5N' or i == 'F201N' or i =='F205N':
                            if len(seq) < 27: #6-21 nt plus minus 6
                                dict_filt3_rest[i][seq]=dict_rounds[i][seq]

                            
                                                        
for i in rounds_names:
    count_range[i]['fraction_greater5len_total'] = count_range[i]['greater5len']/count_range[i]['total']
    count_range[i]['fraction_filtered3_greater5len_total']=count_range[i]['filtered3_greater5len']/count_range[i]['total']
    count_range[i]['fraction_in_range_total'] = count_range[i]['in_range']/count_range[i]['total']
    count_range[i]['fraction_in_range_filtered3_total']=count_range[i]['in_range_filtered3']/count_range[i]['total']
    

sub_selection = ['N73L','N73H'] # pool H and I

dict_comb = {}
for i in sub_selection:
    for seq in dict_filt3_range_rounds[i]:
        if seq not in dict_comb:
            dict_comb[seq]=dict_filt3_range_rounds[i][seq]
        else:
            dict_comb[seq]+=dict_filt3_range_rounds[i][seq]

family_codon_start = ['CT','GT','TC','CC','AC','GC','CG','GG']
nt = ['A','C','G','T']
all_codon =[]
family_codon =[]
for n1 in nt:
    for n2 in nt:
        for n3 in nt:
            all_codon.append(n1+n2+n3)
for f in range(len(family_codon_start)):
    for n in nt:
        full = family_codon_start[f]+n
        #replace family_codon[f] in family_codon with full
        family_codon.append(full)
family_codon

#get reverse complement of sequences in family_codon
family_codon_reversecomp = []
for seq in family_codon:
    seq = seq[::-1]
    seq = seq.replace('A','t')
    seq = seq.replace('T','a')
    seq = seq.replace('C','g')
    seq = seq.replace('G','c')
    seq = seq.upper()
    family_codon_reversecomp.append(seq)
family_codon_reversecomp

#make list of codons from all_codons that are not in non_family_codon_reversecomp 
not_family_codon_reversecomp = []
for codon in all_codon:
    if codon not in family_codon_reversecomp:
        not_family_codon_reversecomp.append(codon)

#counts the longest stretch of translated peptide basedon family codons - need to replace dict_random with dict_comb for sample
lengths_peptide = {}

for seq in dict_comb:
    
    bool_family =[]
    one_list =[]
    for i in range(len(seq)):
        if i%3 ==0:
            codon = seq[i:i+3]
            if codon in family_codon:
                one_list.append(1)
            else:
                one_list.append(0)

    max_val = len(max("".join(map(str, one_list)).split("0")))
    if max_val in lengths_peptide:
        lengths_peptide[max_val]+=dict_comb[seq]
    else:
        lengths_peptide[max_val]=dict_comb[seq]


import pandas as pd
#convert lengths_peptide to dataframe
lengths_peptide_df = pd.DataFrame.from_dict(lengths_peptide, orient='index')
lengths_peptide_df = lengths_peptide_df.reset_index()
lengths_peptide_df = lengths_peptide_df.rename(columns={'index':'Length', 0:'Frequency'})
lengths_peptide_df['Frequency'] = lengths_peptide_df['Frequency']/sum(lengths_peptide_df['Frequency'])
#sort by length
lengths_peptide_df = lengths_peptide_df.sort_values(by=['Length'])
#reset index
lengths_peptide_df = lengths_peptide_df.reset_index(drop=True)

#save lengths_peptide_df to csv
lengths_peptide_df.to_csv('Peptide_length_family_box_length_HI.csv', index=False)
                    


#new dict_random only based on length
# random sequences based on length distribution in dataset
import numpy as np
nucleotides = ['G', 'T', 'C', 'A']
nucleotides2 = ['T','C','A','G']
codon ={}
dict_random = {}
len_dict ={}

total = sum(dict_comb.values())
    
for seq in dict_comb:
    len_s =len(seq)
    if len_s in len_dict:
        len_dict[len_s]+=dict_comb[seq]/total
    else:
        len_dict[len_s]=dict_comb[seq]/total

#list of keys in len_dict
len_dict_keys = list(len_dict.keys())
#list of values in len_dict
len_dict_values = list(len_dict.values())

codons = []
for n1 in nucleotides:
    for n2 in nucleotides2:
        for n3 in nucleotides2:
            codons.append(n1+n2+n3)
    
for j in range(len(dict_comb)):
    length = np.random.choice(len_dict_keys, p=len_dict_values, size=1)
    codon_number = length/3
    sequence =''
    for c in range(int(codon_number)):
        codon = np.random.choice(codons, size=1)
        sequence += codon[0]
    
    if  sequence not in dict_random: #adds to dict if seq not there
        dict_random[sequence] = 1
        
    else:
        dict_random[sequence] +=1


#counts the longest stretch of translated peptide basedon family codons - need to replace dict_random with dict_comb for sample
lengths_peptide = {}

for seq in dict_random:
    
    bool_family =[]
    one_list =[]
    for i in range(len(seq)):
        if i%3 ==0:
            codon = seq[i:i+3]
            if codon in family_codon:
                one_list.append(1)
            else:
                one_list.append(0)

    max_val = len(max("".join(map(str, one_list)).split("0")))
    if max_val in lengths_peptide:
        lengths_peptide[max_val]+=dict_random[seq]
    else:
        lengths_peptide[max_val]=dict_random[seq]
                    
import pandas as pd
#convert lengths_peptide to dataframe
lengths_peptide_df = pd.DataFrame.from_dict(lengths_peptide, orient='index')
lengths_peptide_df = lengths_peptide_df.reset_index()
lengths_peptide_df = lengths_peptide_df.rename(columns={'index':'Length', 0:'Frequency'})
lengths_peptide_df['Frequency'] = lengths_peptide_df['Frequency']/sum(lengths_peptide_df['Frequency'])
#sort by length
lengths_peptide_df = lengths_peptide_df.sort_values(by=['Length'])
#reset index
lengths_peptide_df = lengths_peptide_df.reset_index(drop=True)

#save lengths_peptide_df to csv
lengths_peptide_df.to_csv('Peptide_length_family_box_length_HI_codon_random.csv', index=False)