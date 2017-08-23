import sys

extract_targets = open(sys.argv[2],'r').read().rstrip().split("\n")

out = open(sys.argv[1]+".extracted", 'w')
lines = open(sys.argv[1], 'r').readlines()
out.write(lines[0])
for i in lines[1:]:
    if i.split("\t")[0] in extract_targets:
        out.write(i)