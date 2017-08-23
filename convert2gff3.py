import sys

def parse_gtf(gtf_file, outfile):
    out = open(outfile, "w")
    out.write("##gff-version    3\n")
    ID = None
    all_pos = []
    for i in open(gtf_file, "r"):
        if i[0] != "#":
            values = i.rstrip().split("\t")
            if ID != values[8]:
                out.write("##sequence-region\n")
            ID = values[8]
            new = values[:8]
            new[0] = new[0].split(" ")[0]
            if values[2] == "gene":
                ID = values[8]
                new.append(values[8].replace("gene_id \"","ID=").replace("\"; transcript_id \"", ";Name=").rstrip("\";"))
                out.write("\t".join(new) + "\n")
            elif values[2] == "exon" or values[2] == "CDS":    
                new.append(values[8].replace("gene_id \"","ID=").replace("\"; transcript_id \"", ";Parent=").rstrip("\";"))
                out.write("\t".join(new) + "\n")


if __name__ == "__main__":
    parse_gtf(sys.argv[1], sys.argv[1] + ".gff")