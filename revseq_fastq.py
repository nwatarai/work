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
for i in open(sys.argv[1],'r'):
	line = counter % 4
	if line == 0:
		output.writelines(i)
	elif line == 1:
		output.writelines(rev(i))
	elif line == 2:
		output.writelines(i)
	elif line == 3:
		output.writelines(i.rstrip()[::-1] + '\n')
	counter += 1
