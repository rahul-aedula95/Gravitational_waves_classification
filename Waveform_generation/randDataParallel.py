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
import random

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

def computeHPHC(argList):
	for apx in ['SEOBNRv2']:
		hp, hc = get_td_waveform(approximant=apx, \
			mass1 = argList[0], \
			mass2 = argList[1], \
			spin1z = argList[2], \
			delta_t = argList[3],\
			f_lower = argList[4])
	return hp, hc

def computeMaxHPHC(argList):
	hp, hc = computeHPHC(argList)
	return max(hp)

def findAmplitude(argListParallel, fast = False):
	if fast == False:
		pool = Pool(multiprocessing.cpu_count())
		# pool = Semaphore(multiprocessing.cpu_count()) # use if above doesn't work

		print '\nComputation START (Evaluation Stage)'
		print 'Number of cores (usage): '+str(multiprocessing.cpu_count())
		return pool.map(computeMaxHPHC, argListParallel)
	else:
		pass		


def randGenValues(size):
	m1 = list()
	spz = list()
	del_t = list()
	f_low = list()
	argListParallel = list()
	temp = list()
	spz = 0.0
	del_t = 1.0/4096
	f_low = 60

	mat = dict()
	mat['m1'] = list()
	mat['m2'] = list()
	mat['spz'] = list()
	mat['del_t'] = list()
	mat['f_low'] = list()

	m2 = random.sample(range(10, 10+size), size)
	for val in m2:
		temp = list()
		while True:
			m = random.randint(val, val+size) 
			if m in m1:
				continue
			m1.append(m)
			break
		mat['m1'].append(val)
		mat['m2'].append(m)
		mat['spz'].append(spz)
		mat['del_t'].append(del_t)
		mat['f_low'].append(f_low)
		temp.append(val)
		temp.append(m)
		temp.append(spz)
		temp.append(del_t)
		temp.append(f_low)
		argListParallel.append(temp)
	return mat, argListParallel

if __name__ == "__main__":
	# Argument to randGenValues - number of samples
	start = time.time()
	fastToggle = True
	matrix, argListParallel = randGenValues(45)
	matrix['amplitude'] = findAmplitude(argListParallel, fastToggle)
	execTime = time.time() - start

	# Check if path to the file exists or if the file exists
	filePath = "./Data/randMatData.csv"
	incFile = 1
	while (os.path.exists(filePath)):
		if incFile == 1:
			filePath = filePath[:-4]+str(incFile)+".csv"
		else:
			filePath = filePath[:-5]+str(incFile)+".csv"
		incFile += 1
	output = open(filePath, "w+")

	# Output the generated data to a file
	output.write("amplitude, "+"m1, "+"m2, "+"spz, "+"del_t, "+"f_low,\n")
	for i in xrange(len(matrix['amplitude'])):
		output.write(str(matrix['amplitude'][i])+", "+str(matrix['m1'][i])+", "+str(matrix['m2'][i])+\
			", "+str(matrix['spz'][i])+", "+str(matrix['del_t'][i])+", "+str(matrix['f_low'][i])+",\n")

	output.close()

	print "\nExecuted Successfully"
	print "Time Taken: "+str(execTime)+"\n"