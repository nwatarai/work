import sys
argvs = sys.argv
if len(argvs)<3:
    print "Usage $python %s fasta nmr_of_splits"
    quit()
fasta = open(argvs[1],'r').read().rstrip().split(">")
del fasta[0]
nmr_of_data = len(fasta)
nmr_of_splits = int(argvs[2])

print nmr_of_data
print nmr_of_splits

if nmr_of_data<nmr_of_splits:
    print "data too small"
    quit()

if nmr_of_data%nmr_of_splits == 0:
    nmr_of_data_in_one_split = int(nmr_of_data/nmr_of_splits)
else:
    nmr_of_data_in_one_split = int(nmr_of_data/nmr_of_splits)+1

i=0
j=1
for f in fasta:
    if i%nmr_of_data_in_one_split == 0:
        file = open("query"+str(j),'w')
        j += 1
    file.writelines(">"+f)
    if (i+1)%nmr_of_data_in_one_split == 0:
        file.close
    i += 1
try:
    file.close
except:
    print ""