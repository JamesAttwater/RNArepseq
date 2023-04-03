gam="GCAGAGGCGGCAGCCTTCGGTGGC"
l1 = list(gam)
tally = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

f=open('C:\\Users\\jamie\\Desktop\\Sequencing\\211221_RyD_replication\\results\\0Dg.fasta', 'r')
lines=f.readlines()
print(len(lines))

j = 1
while j<len(lines):
    seq = lines[j]
    l2 = list(seq)
    index = 0
    count = 0
    if len(l2)<26:
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
        j+=2
else:
    print(tally)
    print(sum(tally))