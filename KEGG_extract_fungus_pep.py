import sys
import gzip

taxonomy = sys.argv[1]
# pep is gunzipped
pep = sys.argv[2]

flag = False
fungi = []
for i in open(taxonomy, "r"):
	if flag:
		if i[0] != "#":
			fungi.append(i.split("\t")[1])
		elif i[:3] != "## ":
			break
	if i.rstrip() == "## Fungi":
		flag = True

fasta = gzip.open(pep, "r").split(">")[1:]
out = open(sys.argv[2].replace(".gz", ".fungi"), "w")

for f in fasta:
	f[:3] in fungi:
		out.write(">" + f)