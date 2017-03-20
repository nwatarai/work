import sys
import numpy as np
import pandas as pd
from random import random, shuffle

argvs = sys.argv
if len(argvs) < 3:
	print "Usage: python %s nmr_of_chm nmr_of_spc nmr_of_reactions output_filename" %argvs[0]
	quit()

nmr_of_chm = int(argvs[1])
nmr_of_spc = int(argvs[2])
nmr_of_reactions = int(argvs[3])
fix = argvs[4]

if nmr_of_reactions >= (nmr_of_chm * (nmr_of_chm -1)):
	sys.exit("nmr_of_reactions must be less than possible")
if nmr_of_reactions < nmr_of_chm - 1:
	sys.exit("nmr_of_reactions must be more than (nmr_of_chm - 1)")

def chm():
	S = 0
	while S <= 0:
		chms = np.random.gamma(1.0,1.0,nmr_of_chm)
		chms[-1] = 0.0
		chms = np.sort(chms)[::-1]
		S = np.sum(chms)
	chms = chms / np.sum(chms) * nmr_of_chm

	index = np.arange(0,nmr_of_chm)
	columns = ["initial"]
	pd.DataFrame(
		chms, index=index, columns=columns).to_csv("param_chm_" + fix + ".csv")

def param_sets(sort=False):
	initial = np.random.gamma(1,0.1,nmr_of_spc) / 10000000
	growthk = np.random.gamma(10,0.1,nmr_of_spc) 
	#growthk = np.random.gamma(2,1,nmr_of_spc)
	transk = growthk
	#sa = np.random.uniform(0.1,0.4,nmr_of_spc)
	sa = np.random.uniform(0.2,0.3,nmr_of_spc)
	#si_rate = np.random.uniform(0.001,0.01,nmr_of_spc)
	si_rate = np.random.uniform(0.01,0.1,nmr_of_spc)

	si = sa * si_rate

	#sorting
	if sort:
		growthk = np.sort(growthk)
		transk = np.sort(transk)
		sa = sa[np.argsort(si)]
		si = np.sort(si) 
	
	index = np.arange(0,nmr_of_spc)
	columns = ["initial_i","growth_k","trans_k","sa","si"]
	data = np.transpose([initial, growthk, transk, sa, si])
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_" + fix + ".csv")

def reaction_sets(random_network=True):
	def dicide_output_index(input_index):
		if random_network:
			output_index = input_index
			while output_index == input_index:
				output_index = np.random.randint(nmr_of_chm)
			return output_index
		else:
			if input_index == nmr_of_chm - 1:
				return input_index
			point = np.random.randint( int((nmr_of_chm - input_index) * (nmr_of_chm - input_index - 1) * 0.5 ) ) 
			#for i in range(input_index + 1, nmr_of_chm):
			for i in range(input_index + 1, nmr_of_chm)[::-1]:
				point -=  ( nmr_of_chm - i )
				if point < 0: 
					return i
			sys.exit("error")

	reactions = [[] for i in range(nmr_of_chm)]
	for i in range(nmr_of_chm - 1):
		input_index = i
		output_index = dicide_output_index(input_index)
		reactions[input_index].append(output_index)
	if random_network:
		output_index = dicide_output_index(nmr_of_chm - 1)
		reactions[nmr_of_chm - 1].append(output_index)
		rest = nmr_of_reactions - nmr_of_chm
	else:
		rest = nmr_of_reactions - nmr_of_chm + 1
	for i in range(rest):
		flag = True
		while flag:
			if random_network:
				input_index = np.random.randint(nmr_of_chm)
			else:
				input_index = np.random.randint(nmr_of_chm - 1)
			output_index = dicide_output_index(input_index)
			if output_index not in reactions[input_index]:
				flag = False
		reactions[input_index].append(output_index)
	data = []
	for i in range(nmr_of_chm):
		for value in reactions[i]:
			data.append([i, value])

	index = np.arange(0,nmr_of_reactions)
	columns = ["input_index","output_index"]
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_reaction_" + fix + ".csv")

def _1_reaction_sets():
	data = []
	if nmr_of_reactions != (nmr_of_chm - 2) *2:
		sys.exit("In this option, nmr_of_reactions must be (nmr_of_chm - 2) * 2")
	for i in range(1, nmr_of_chm - 1):
		input_index = 0
		interval_index = i
		output_index = nmr_of_chm - 1
		data.append([input_index, interval_index])
		data.append([interval_index, output_index])

	index = np.arange(0,nmr_of_reactions)
	columns = ["input_index","output_index"]
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_reaction_" + fix + ".csv")

#all chemical transformed to one available chemical
def _2_reaction_sets():
	data = []
	if nmr_of_reactions != nmr_of_chm - 1:
		sys.exit("In this option, nmr_of_reactions must be (nmr_of_chm - 1)")
	for i in range(nmr_of_chm - 2):
		input_index = i
		output_index = nmr_of_chm - 2
		data.append([input_index, output_index])
	data.append([nmr_of_chm -2, nmr_of_chm -1])

	index = np.arange(0,nmr_of_reactions)
	columns = ["input_index","output_index"]
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_reaction_" + fix + ".csv")

#no chemical interaction
def _3_reaction_sets():
	data = []
	if nmr_of_reactions != nmr_of_chm - 1:
		sys.exit("In this option, nmr_of_reactions must be (nmr_of_chm - 1)")
	for i in range(nmr_of_chm - 1):
		input_index = i
		output_index = nmr_of_chm - 1
		data.append([input_index, output_index])

	index = np.arange(0,nmr_of_reactions)
	columns = ["input_index","output_index"]
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_reaction_" + fix + ".csv")

def _3_reaction_sets():
	data = []
	if nmr_of_reactions != nmr_of_chm - 1:
		sys.exit("In this option, nmr_of_reactions must be (nmr_of_chm - 1)")
	for i in range(nmr_of_chm - 1):
		input_index = i
		output_index = nmr_of_chm - 1
		data.append([input_index, output_index])

	index = np.arange(0,nmr_of_reactions)
	columns = ["input_index","output_index"]
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_reaction_" + fix + ".csv")

#each chemical fransformed to any downstream chemcical by some species
def _4_reaction_sets():
	data = []
	if nmr_of_reactions != (nmr_of_chm - 1) * nmr_of_chm /  2:
		sys.exit("In this option, nmr_of_reactions must be (nmr_of_chm - 1) * nmr_of_chm / 2")
	for i in range(nmr_of_chm - 1):
		input_index = i
		for j in range(i + 1, nmr_of_chm):
			output_index = j
			data.append([input_index, output_index])

	index = np.arange(0,nmr_of_reactions)
	columns = ["input_index","output_index"]
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_reaction_" + fix + ".csv")

def _linear_reaction_sets():
	data = []
	if nmr_of_reactions != nmr_of_chm - 1:
		sys.exit("In this option, nmr_of_reactions must be (nmr_of_chm - 1)")
	for i in range(nmr_of_chm - 1):
		input_index = i
		output_index = input_index + 1
		data.append([input_index, output_index])

	index = np.arange(0,nmr_of_reactions)
	columns = ["input_index","output_index"]
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_reaction_" + fix + ".csv")

def metabolism_sets():
	data = np.random.normal(loc=-2,scale=1,size=(nmr_of_spc, nmr_of_reactions))
	data[data<0] = 0
	for i in range(data.shape[0]):
		S = 0
		if np.sum(data[i]) == 0:
			while S == 0:
				data[i] = np.random.normal(loc=-2,scale=1,size=nmr_of_reactions)
				data[i][data[i]<0] = 0
				S = np.sum(data[i])

	index = np.arange(0,nmr_of_spc)
	columns = np.arange(0,nmr_of_reactions)
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_metabolism_" + fix + ".csv")

"""
def _param_metabolism_sets():
	initial = np.random.gamma(1,0.1,nmr_of_spc) / 10000000
	growthk = np.random.gamma(1,0.5,nmr_of_spc) + np.ones(nmr_of_spc)*0.5
	transk = growthk
	sa = []

	data = np.zeros((nmr_of_spc, nmr_of_reactions))
	for i in range(nmr_of_spc):
		index = np.random.randint(nmr_of_reactions)
		data[i][index] = 1.0
		if index == nmr_of_chm - 2:
			sa.append(0.5 / np.sqrt(float(index + 10)))

	si_rate = np.random.uniform(0.01,0.1,nmr_of_spc)
	si = sa * si_rate	

	index = np.arange(0,nmr_of_spc)
	columns = np.arange(0,nmr_of_reactions)
	pd.DataFrame(
		data, index=index, columns=columns).to_csv("param_metabolism_" + fix + ".csv")

	param_index = np.arange(0,nmr_of_spc)
	param_columns = ["initial_i","growth_k","trans_k","sa","si"]
	param_data = np.transpose([initial, growthk, transk, sa, si])
	pd.DataFrame(
		param_data, index=param_index, columns=param_columns).to_csv("param_" + fix + ".csv")
"""

if __name__ == "__main__":
	chm()
	param_sets(sort=True)
	#reaction_sets(random_network=False)
	#_linear_reaction_sets()
	_3_reaction_sets()
	metabolism_sets()
