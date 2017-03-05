""" Code to generate and record the amplitude at coalescence time.
"""

import pylab
from pycbc.waveform import get_td_waveform
import numpy
import myconstants
import sys
import os
import time
import multiprocessing
from multiprocessing import Pool
from multiprocessing import Semaphore
import copy

""" 
	mass of first body = m1
	mass of second body = m2
	spin1z = spz
	delta_t = del_t
	f_lower = f_low

	(Default) range of m1 = (1, 100), increment 1
	(Default) range of m2 = (1,100), increment 1
	(Default) range of spz = (0.0, 10.0), increment 0.1
	(Default) del_t = 1.0/4096 (const)
	(Default) f_low = const (60)
"""

# Define a function frange() for floating point increments
def frange(start, stop, step):
	x = start
	temp= []
	while x < stop:
		temp.append(x)
		x += step
	return temp

def computeHPHC(argList):
	for apx in ['SEOBNRv2']:
		hp, hc = get_td_waveform(approximant=apx, \
			mass1 = argList[0], \
			mass2 = argList[1], \
			spin1z = argList[2], \
			delta_t = argList[3],\
			f_lower = argList[4])
	return max(hp)

def genDesignMatrix(m1Range, m2Range, spzRange, del_tRange, f_lowRange, strict = False):
	try:
		if len(del_t) not in [1, 3] or len(f_low) not in [1, 3]:
			raise Exception('Syntax for arguments in genDesignMatrix()')
	except Exception:
		print 'del_t and f_low could either have a single value or a range, \
			in range specify min, max and increment'
	EPSILON = myconstants.epsilon
	if len(del_tRange) == 1:
		del_tRange.append(del_tRange[0] + EPSILON)
		del_tRange.append(EPSILON)
	if len(f_lowRange) == 1:
		f_lowRange.append(f_lowRange[0] + EPSILON)
		f_lowRange.append(EPSILON)

	designMatrix = dict()
	designMatrix['amplitude'] = list()
	designMatrix['m1'] = list()
	designMatrix['m2'] = list()
	designMatrix['spz'] = list()
	designMatrix['del_t'] = list()
	designMatrix['f_low'] = list()
	# Size of the matrix is six because of,
	# m1, m2, spz, del_t, f_low and the amplitude (during coalescence)
	del_tRange = frange(del_tRange[0], del_tRange[1], del_tRange[2])
	f_lowRange = frange(f_lowRange[0], f_lowRange[1], f_lowRange[2])

	print '\nComputation START (Generation Stage)'
	argListParallel = list()
	loop = 1
	for m1 in xrange(m1Range[0], m1Range[1], m1Range[2]):
		for m2 in xrange(m2Range[0], m2Range[1], m2Range[2]):
			if strict == True:
				if m2 > m1:
					continue
			for spz in frange(spzRange[0], spzRange[1], spzRange[2]):
				for del_t in del_tRange:
					for f_low in f_lowRange:		
						print 'Entered loop: '+str(loop)
						loop += 1
						temp = list()
						designMatrix['m1'].append(m1)
						designMatrix['m2'].append(m2)
						designMatrix['spz'].append(spz)
						designMatrix['del_t']. append(del_t)
						designMatrix['f_low'].append(f_low)
						temp.append(m1)
						temp.append(m2)
						temp.append(spz)
						temp.append(del_t)
						temp.append(f_low)
						argListParallel.append(temp)
	pool = Pool(multiprocessing.cpu_count())
	# pool = Semaphore(multiprocessing.cpu_count()) # use if above doesn't work
	print '\nComputation START (Evaluation Stage)'
	print 'Number of cores (usage): '+str(multiprocessing.cpu_count())
	response = pool.map(computeHPHC, argListParallel)
	designMatrix['amplitude'] = copy.deepcopy(response)
	return designMatrix

if __name__ == "__main__":
	# Create and assign default range
	# NOTE: del_t and f_low could be provided the range values.
	m1Range = [10, 20, 1]
	m2Range = [10, 20, 1]
	spzRange = [0.1, 0.2, 0.1]
	del_tRange = [1.0/4096]
	f_lowRange = [60]

	# Execute the function and evaluate the time
	# NOTE: If strict == True, then always m2 >= m1
	strictToggle = True
	start = time.time()
	matrix = genDesignMatrix(m1Range, m2Range, spzRange, del_tRange, f_lowRange, strict=strictToggle)
	execTime = time.time() - start

	# Check if path to the file exists or if the file exists
	filePath = "./Data/matData.csv"
	incFile = 1
	while (os.path.exists(filePath)):
		if incFile == 1:
			filePath = filePath[:-4]+str(incFile)+".csv"
		else:
			filePath = filePath[:-5]+str(incFile)+".csv"
		incFile += 1
	output = open(filePath, "w+")
	outputDoc = open(filePath[:-4]+"Docs.md", "w+")

	# Output the information about the data
	outputDoc.write("============\n Data Range\n============\n\n")
	outputDoc.write("m1: ("+str(m1Range[0])+", "+str(m1Range[1])+", "+str(m1Range[2])+")\n")
	outputDoc.write("m2: ("+str(m2Range[0])+", "+str(m2Range[1])+", "+str(m2Range[2])+")\n")
	outputDoc.write("spz: ("+str(spzRange[0])+", "+str(spzRange[1])+", "+str(spzRange[2])+")\n")
	del_tL = len(del_tRange)
	f_lowL = len(f_lowRange)
	outputDoc.write("del_t: (")
	for i in xrange(del_tL):
		if i == del_tL-1:
			outputDoc.write(str(del_tRange[i])+")\n")
		else:
			outputDoc.write(str(del_tRange[i])+", ")
	outputDoc.write("f_low: (")
	for i in xrange(f_lowL):
		if i == f_lowL-1:
			outputDoc.write(str(f_lowRange[i])+")\n")
		else:
			outputDoc.write(str(f_lowRange[i])+", ")
	outputDoc.write("Execution Time: "+str(execTime)+"\n")
	outputDoc.write("Strict Toggle: "+str(strictToggle)+"\n")
	outputDoc.close()

	# Output the generated data to a file
	output.write("#amplitude, "+"#m1, "+"#m2, "+"#spz, "+"#del_t, "+"#f_low,\n")
	for i in xrange(len(matrix['amplitude'])):
		output.write(str(matrix['amplitude'][i])+", "+str(matrix['m1'][i])+", "+str(matrix['m2'][i])+\
			", "+str(matrix['spz'][i])+", "+str(matrix['del_t'][i])+", "+str(matrix['f_low'][i])+"\n")

	output.close()

	print "\nExecuted Successfully"
	print "Time Taken: "+str(execTime)+"\n"
