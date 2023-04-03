gam="GCAGAGGCGGCAGCCTTCGGTGGC"
l1 = list(gam)
tally = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

import random
triplets = ("GCA", "GAG", "GCG", "GCC", "TTC", "GGT", "GGC", "CGA", "AGG", "CTG", "CCG", "CCT", "CAT")

j = 1
while j < 20000:
    nhy = []
    times = 0
    while times<8:
        nhy.append(random.choice(triplets))
        times+=1
    else:
        seq = "".join(nhy)
        l2 = list(seq)
        index = 0
        count = 0
        while index<24:
            if l1[index] == l2[index]:
                count+=1
            else:
                count+=0
            index+=1
        else:
            tally[count]+=1
        j+=2
else:
    print(tally)
    print(sum(tally))