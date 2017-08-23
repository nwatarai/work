import sys

def n2p(seq):
    DNA2Protein = {
            'TTT' : 'F', 'TCT' : 'S', 'TAT' : 'Y', 'TGT' : 'C',
            'TTC' : 'F', 'TCC' : 'S', 'TAC' : 'Y', 'TGC' : 'C',
            'TTA' : 'L', 'TCA' : 'S', 'TAA' : '*', 'TGA' : '*',
            'TTG' : 'L', 'TCG' : 'S', 'TAG' : '*', 'TGG' : 'W',

            'CTT' : 'L', 'CCT' : 'P', 'CAT' : 'H', 'CGT' : 'R',
            'CTC' : 'L', 'CCC' : 'P', 'CAC' : 'H', 'CGC' : 'R',
            'CTA' : 'L', 'CCA' : 'P', 'CAA' : 'Q', 'CGA' : 'R',
            'CTG' : 'L', 'CCG' : 'P', 'CAG' : 'Q', 'CGG' : 'R',

            'ATT' : 'I', 'ACT' : 'T', 'AAT' : 'N', 'AGT' : 'S',
            'ATC' : 'I', 'ACC' : 'T', 'AAC' : 'N', 'AGC' : 'S',
            'ATA' : 'I', 'ACA' : 'T', 'AAA' : 'K', 'AGA' : 'R',
            'ATG' : 'M', 'ACG' : 'T', 'AAG' : 'K', 'AGG' : 'R',

            'GTT' : 'V', 'GCT' : 'A', 'GAT' : 'D', 'GGT' : 'G',
            'GTC' : 'V', 'GCC' : 'A', 'GAC' : 'D', 'GGC' : 'G',
            'GTA' : 'V', 'GCA' : 'A', 'GAA' : 'E', 'GGA' : 'G',
            'GTG' : 'V', 'GCG' : 'A', 'GAG' : 'E', 'GGG' : 'G'
    }

    p = ''
    for i in range(int(len(seq) / 3)):
        codon = seq[3 * i : 3 * i + 3]
        if codon in DNA2Protein:
            p += DNA2Protein[codon]
        else:
            p += '*'
    return p

def threepack(seq):
    p1 = n2p(seq)
    p2 = n2p(seq[1:])
    p3 = n2p(seq[2:])
    P = [p1, p2, p3]
    c1 = p1.count('*')
    c2 = p2.count('*')
    c3 = p3.count('*')
    C = [c1, c2 ,c3]
    m = min(C)
    index = C.index(m)
    return P[index], seq[index:], m

def revseq(seq):
    rev = seq[::-1].upper()
    rev = rev.replace("C","g").replace("G","c").replace("A","t").replace("T","a")
    return rev.upper()

def parse_gtf(gtf_file):
    contig = None
    ID = None
    all_pos = []
    pos = []
    for i in open(gtf_file, "r"):
        if i[0] == "#" and len(pos) > 0:
            all_pos.append({"contig": contig, "strand": strand, "pos": pos, "ID": ID})
            pos = []
        elif i[0] != "#":
            values = i.rstrip().split("\t")
            feat = values[2]
            if feat == "CDS" or feat == "coding_exon" or feat == "cds":
                contig = values[0]
                start = int(values[3])
                end = int(values[4])
                strand = values[6]
                ID = values[8]
                pos.append([start, end])
    return all_pos

def parse_fasta(fasta_file):
    fastas = open(fasta_file, "r").read().split(">")[1:]
    dic = {}
    for f in fastas:
        s = f.split("\n")
        ID = s[0].split(" ")[0]
        seq = "".join(s[1:])
        dic[ID] = seq.upper()
    return dic

def chain_seq(seq, pos):
    out = ""
    for p in pos:
        out += seq[p[0] - 1 : p[1]]
    return out

def _test(seq, pos):
    def _revseq(seq):
        rev = seq.replace("C","g").replace("G","c").replace("A","t").replace("T","a")
        return rev.upper()
    out = ""
    for p in pos:
        out += seq[- p[1]: - p[0]]
    return _revseq(out)

def main(gtf_file, fasta_file, n_output, p_output):
    all_pos = parse_gtf(gtf_file)
    fasta = parse_fasta(fasta_file)
    n_out = open(n_output, "w")
    p_out = open(p_output, "w")
    for i, pos in enumerate(all_pos):
        seq = fasta[pos["contig"]]
        if pos["strand"] == "+":
            seq = chain_seq(seq, pos["pos"])[:-3]
        elif pos["strand"] == "-":
            seq = revseq(chain_seq(seq, pos["pos"]))[:-3]
            #seq = _test(seq, pos["pos"])[:-3]
        else:
            raise ValueError("gff/gtf file may be incorrect")
            
        pro = n2p(seq)
        if "*" in pro or len(seq) % 3 != 0:
            pro, seq, m = threepack(seq)
            if m == 0:
                gene_name = ">" + pos["contig"] + " gene_" + str(i) + "_incomplete\n"
                n_out.write(gene_name + seq + "\n")
                p_out.write(gene_name + pro + "\n")
            else:
                print "gene_%i may be frame-shifted" % i 
                print pos
                print "\n"
        else:
            gene_name = ">" + pos["contig"] + " gene_" + str(i) + "\n"
            n_out.write(gene_name + seq + "\n")
            p_out.write(gene_name + pro + "\n")
    n_out.close()
    p_out.close()

if __name__ == "__main__":
    try:
        gtf_file = sys.argv[1]
        fasta_file = sys.argv[2]
    except:
        print "Usage: python %s GFF FASTA" % sys.argv[0]
        quit()
    n_output = fasta_file + ".cds.fna"
    p_output = fasta_file + ".cds.faa"
    main(gtf_file, fasta_file, n_output, p_output)


