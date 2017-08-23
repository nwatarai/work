import sys

def main(fasta_file, output):
    fastas = open(fasta_file, "r").read().split(">")[1:]
    out = open(output, "w")
    for f in fastas:
        if not "hypothetical" in f: 
        	out.write(">" + f)

if __name__ =="__main__":
	main(sys.argv[1], sys.argv[1]+".delhyp.faa")