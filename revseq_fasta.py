import sys

pair = {'A':'T', 'a':'t', 
		'T':'A', 't':'a',
		'G':'C', 'g':'c',
		'C':'G', 'c':'g',
		'N':'N', 'n':'n'}
def complement(char):
	return pair[char]

def rev(string):
	string = string.rstrip()
	l = list(string)
	l.reverse()
	m = map(complement,l)
	return "".join(list(m)) + '\n'

counter = 0
output = open(sys.argv[1] + '.rev','w')

fasta = open(sys.argv[1], 'r').read().split('>')[1:]

for f in fasta:
	head = '>' + f.split('\n')[0]
	seq = "".join(f.split('\n')[1:])
	seq = rev(seq)
	output.write(head + "\n" + seq + '\n')

