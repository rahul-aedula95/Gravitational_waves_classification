from pycbc.waveform import get_td_waveform
import matplotlib.pyplot as plt
import pylab
import time
from math import log
from math import exp
import numpy as np
from scipy.optimize import curve_fit
import myconstants
from scipy.interpolate import spline
from scipy.interpolate import interp1d
from scipy.interpolate import BarycentricInterpolator
from scipy.interpolate import KroghInterpolator
from scipy import interpolate
import copy
import myconstants
import multiprocessing
from multiprocessing import Pool
from multiprocessing import Semaphore

def computeHPHC(argList):
	for apx in ['SEOBNRv2']:
		hp, hc = get_td_waveform(approximant=apx, \
			mass1 = argList[0], \
			mass2 = argList[1], \
			spin1z = argList[2], \
			delta_t = argList[3],\
			f_lower = argList[4])
	return hp, hc 

def writeData1(data, output):
	keys = data.keys()
	keys = ['y', 'x']
	for i, key in enumerate(keys):																																																																																																																																																																																																																																																																												
		if i == len(keys)-1:
			output.write(str(key)+"\n")
		else:
			output.write(str(key)+", ")
		# data[key] = data[key][::-1]

	for i in xrange(len(data[keys[0]])):
		for j, key in enumerate(keys):
			if j == len(keys)-1:
				output.write(str(data[key][i])+"\n")
			else:
				output.write(str(data[key][i])+"	")
	return

def writeData(data, output):

	keys = data.keys()
	# keys = ['posTime', 'posAmp']
	for i, key in enumerate(keys):																																																																																																																																																																																																																																																																												
		if i == len(keys)-1:
			output.write(str(key)+"\n")
		else:
			output.write(str(key)+", ")

	for i in xrange(len(data[keys[0]])):
		for j, key in enumerate(keys):
			if j == len(keys)-1:
				output.write(str(data[key][i])+"\n")
			else:
				output.write(str(data[key][i])+", ")
	return

def peaks(amp, delT):
	ampPeaks = dict()
	ampPeaks['amp'] = list()
	ampPeaks['time'] = list()
	ampPeaks['posAmp'] = list()
	ampPeaks['negAmp'] = list()
	ampPeaks['posTime'] = list()
	ampPeaks['negTime'] = list()

	# delT = mat['del_t'][0]
	time = delT
	amp = list(amp)
	initValue = max(amp)
	coalesceIndex = amp.index(initValue)

	# True = amplitude (ascending)
	# False = amplitude (descending)
	Toggle = True

	for value in amp[coalesceIndex::-1]:
		time -= delT
		if Toggle == True:
			if value > initValue:
				initValue = value
				continue
		else:
			if value < initValue:
				initValue = value
				continue
		if Toggle == True:
			ampPeaks['posAmp'].append(initValue)
			ampPeaks['posTime'].append(time)
		else:
			ampPeaks['negAmp'].append(initValue)
			ampPeaks['negTime'].append(time)

		Toggle = not Toggle
		ampPeaks['amp'].append(initValue)
		ampPeaks['time'].append(time)
	return ampPeaks

def regFunc(x, a, b, c):
	return a*np.exp(b*x)+c

def predict(x, a, b, c):
	return [(1.0/b)*log((t-c)/a) for t in x]

def regFuncCoalesce(x, a, b):
	return a*(x**b)

def predictCoalesce(x, a, b):
	return [(t/a)**(1.0/b) for t in x]

#==================================================================================
# WAVE ANALYSIS FUNCTION
#==================================================================================
def waveAnalysis(argList):
	hp, hc = computeHPHC([argList[0], argList[1], argList[2], argList[3], argList[4]])
	h = [t*10**19 for t in hp]
	ampPeaks = peaks(h, argList[3])

	y = np.array(ampPeaks['posAmp'])
	x = np.array(ampPeaks['posTime'])	

	try:
		popt, pcov = curve_fit(regFunc, y, x, [-1, -1, -1])
	except Exception:
		return None
	a0 = popt[0]
	b0 = popt[1]
	c0 = popt[2]
	cleanTime = [t for t in ampPeaks['posTime'] if (t-c0)/a0 > 0.0]

	predAmp = predict(cleanTime, a0, b0, c0)

	# CHECK THIS
	y = np.array([-t for t in ampPeaks['negAmp']])
	x = np.array(ampPeaks['negTime'])

	try:
		popt, pcov = curve_fit(regFunc, y, x, [-1, -1, -1])
	except Exception:
		return None
	a = popt[0]
	b = popt[1]
	c = popt[2]
	cleanTimeN = [t for t in ampPeaks['negTime'] if (t-c)/a > 0.0]

	predAmpN = [-t for t in predict(cleanTimeN, a, b, c)]

	# CHECK THIS
	# Separating the non-coalescence Period.
	numEleCleaned = len(ampPeaks['posAmp'])-len(predAmp)
	init = 0
	while ampPeaks['posAmp'][numEleCleaned+init] < predAmp[init]:
		init+=1
	# Marks the start of the nonCoalesce phase for poAmp
	nonCoalescePeriod = numEleCleaned+init
	predAmp = predAmp[init::][::-1]
	cleanTime = cleanTime[init::][::-1]

	numEleCleanedN = len(ampPeaks['negAmp'])-len(predAmpN)
	init = 0
	while ampPeaks['negAmp'][numEleCleanedN+init] < predAmpN[init]:
		init+=1
	# Marks the start of the nonCoalesce phase for poAmp
	nonCoalescePeriodN = numEleCleanedN+init
	predAmpN = predAmpN[init::][::-1]
	cleanTimeN = cleanTimeN[init::][::-1]

	if len(predAmp) > len(predAmpN):
		predAmp = predAmp[:len(predAmpN)]
		cleanTime = cleanTime[:len(predAmpN)]
		nonCoalesceLen = len(predAmpN)
	elif len(predAmp) < len(predAmpN):
		predAmpN = predAmpN[:len(predAmp)]
		cleanTimeN = cleanTimeN[:len(predAmpN)]
		nonCoalesceLen = len(predAmp)
	else:
		nonCoalesceLen = len(predAmp)
	nonCoalescePeriod = len(ampPeaks['posTime'])-nonCoalesceLen
	nonCoalescePeriodN = len(ampPeaks['negTime'])-nonCoalesceLen

	time = list()
	amp = list()
	length = max([len(cleanTime), len(cleanTimeN)])
	if abs(cleanTime[0]) > abs(cleanTimeN[0]):
		for t in xrange(length):
			if t < len(cleanTime):
				time.append(cleanTime[t])
				amp.append(predAmp[t])
			if t < len(cleanTimeN):
				time.append(cleanTimeN[t])
				amp.append(predAmpN[t])
	else:
		for t in xrange(length):
			if t < len(cleanTimeN):
				time.append(cleanTimeN[t])
				amp.append(predAmpN[t])
			if t < len(cleanTime):
				time.append(cleanTime[t])
				amp.append(predAmp[t])

	#======================================
	# During Coalescence
	#======================================
	y = ampPeaks['posAmp'][::-1]
	x = ampPeaks['posTime'][::-1]
	slopeChangeIndex = 0
	for i in xrange(len(y)):
		if (y[i+1]-y[i])/(x[i+1]-x[i]) >= 1.0:
			slopeChangeIndex = i
			break
	slopeChangeIndex = len(y)-slopeChangeIndex

	slopeChangeIndex = min(slopeChangeIndex, nonCoalescePeriod)
	slopeChangeIndex = max(slopeChangeIndex, 4)

	# CHANGE INTRODUCED HERE
	temp = len(ampPeaks['posAmp']) - len(cleanTime)
	# y1 = ampPeaks['posAmp'][:slopeChangeIndex][::-1][:-1]
	# x1 = ampPeaks['posTime'][:slopeChangeIndex][::-1][:-1]
	y1 = ampPeaks['posAmp'][:temp][::-1]
	x1 = ampPeaks['posTime'][:temp][::-1]

	tempDict = dict()
	tempDict['y'] = y1
	tempDict['x'] = x1

	y = np.array(y1)
	x = np.array(x1)	

	try:
		popt, pcov = curve_fit(regFuncCoalesce, y, x, [-1, -1])
	except Exception:
		return None
	a1 = popt[0]
	b1 = popt[1]

	x1 = x1[:-1]
	y1 = y1[:-1]

	predAmpCoalesce = predictCoalesce(x1, a1, b1)
	predAmpCoalesce.append(max(ampPeaks['posAmp']))
	x1.append(max(ampPeaks['posTime']))

	# CHANGE INTRODUCED HERE
	y2 = [-t for t in ampPeaks['negAmp'][:temp][::-1]]
	x2 = ampPeaks['negTime'][:temp][::-1]
	y = np.array(y2)
	x = np.array(x2)	

	try:
		popt, pcov = curve_fit(regFuncCoalesce, y, x, [-1, -1])
	except Exception:
		return None
	a = popt[0]
	b = popt[1]

	predAmpNCoalesce = predictCoalesce(x2, a, b)
	predAmpNCoalesce = [-t for t in predAmpNCoalesce]

	timeCoalesce = list()
	ampCoalesce = list()
	length = max([len(x1), len(x2)])	
	for t in xrange(length):
		if t < len(x2):
			timeCoalesce.append(x2[t])
			ampCoalesce.append(predAmpNCoalesce[t])
		if t < len(x1):
			timeCoalesce.append(x1[t])
			ampCoalesce.append(predAmpCoalesce[t])	
	fH_Coal = 1.0/(x1[len(x1)-1] - x1[len(x1)-2])
	fL_Coal = 1.0/(x1[1] - x1[0])
	fH_nCoal = 1.0/(cleanTime[len(cleanTime)-1] - cleanTime[len(cleanTime)-2])
	fL_nCoal = 1.0/(cleanTime[1]-cleanTime[0])
	print 'computed:', [argList[0], argList[1], argList[4], a0, b0, c0, a1, b1]
	return [argList[0], argList[1], argList[4], a0, b0, c0, a1, b1, len(cleanTime), len(x1), \
	fH_Coal, fL_Coal, fH_nCoal, fL_nCoal]



#===================================================================================

if __name__ == "__main__":
	regAnalysis = dict()
	regAnalysis['eA'] = list()
	regAnalysis['eB'] = list()
	regAnalysis['eC'] = list()
	regAnalysis['xA'] = list()
	regAnalysis['xB'] = list()
	regAnalysis['f'] = list()
	regAnalysis['m1'] = list()
	regAnalysis['m2'] = list()
	regAnalysis['cLen'] = list()
	regAnalysis['ncLen'] = list()
	regAnalysis['fH_Coal'] = list()
	regAnalysis['fL_Coal'] = list()
	regAnalysis['fH_nCoal'] = list()
	regAnalysis['fL_nCoal'] = list()

	argListParallel = list()
	delT = 1.0/4096
	# for f in xrange(20, 41):
	# 	for m1 in xrange(10, 31):
	# 		for m2 in xrange(10, 31):
	for f in xrange(20, 22):
		for m1 in xrange(10, 12):
			for m2 in xrange(10, 12):
				if m2 > m1:
					pass
				temp = list()
				temp.append(m1)
				temp.append(m2)
				temp.append(0.0)
				temp.append(delT)
				temp.append(f)
				argListParallel.append(temp)
	pool = Pool(multiprocessing.cpu_count())
	# pool = Semaphore(multiprocessing.cpu_count()) # use if above doesn't work
	print '\nComputation START (Evaluation Stage)'
	print 'Number of cores (usage): '+str(multiprocessing.cpu_count())
	print ''
	response = pool.map(waveAnalysis, argListParallel)
	for res in response:
		if None in res:
			continue
		regAnalysis['eA'].append(res[3])
		regAnalysis['eB'].append(res[4])
		regAnalysis['eC'].append(res[5])
		regAnalysis['xA'].append(res[6])
		regAnalysis['xB'].append(res[7])
		regAnalysis['f'].append(res[2])
		regAnalysis['m1'].append(res[0])
		regAnalysis['m2'].append(res[1])
		regAnalysis['cLen'].append(res[8])
		regAnalysis['ncLen'].append(res[9])
		regAnalysis['fH_Coal'].append(res[10])
		regAnalysis['fL_Coal'].append(res[11])
		regAnalysis['fH_nCoal'].append(res[12])
		regAnalysis['fL_nCoal'].append(res[13])

	outFile = open("./Data/WaveCharacteristics.csv","w+")
	writeData(regAnalysis, outFile)
	outFile.close()

	print ''

	# #================================================================================
	# # REAL ANALYSIS
	# #================================================================================
	# print 'length of Coalescence Period:', len(x1)
	# print 'length of non Coalescence period:', len(cleanTime)
	# #================================================================================

	# plt.figure(figsize=(10,8))
	# plt.scatter(ampPeaks['posTime'], ampPeaks['posAmp'], color='g')
	# plt.scatter(ampPeaks['negTime'], ampPeaks['negAmp'], color='g')
	# plt.plot(hp.sample_times, h, color='b')
	# # plt.scatter(ampPeaks['posTime'], predicted, color = 'r')
	# plt.scatter(cleanTime, predAmp, color = 'r')
	# plt.scatter(cleanTimeN, predAmpN, color ='r')
	# plt.scatter(x1, predAmpCoalesce, color = 'black')
	# plt.scatter(x2, predAmpNCoalesce, color='black')
	# # plt.plot(xnew, splineFit, color = 'orange')
	# plt.plot(time, amp, color = 'orange')
	# plt.plot(timeCoalesce, ampCoalesce, color='orange')
	# plt.xlabel("Time")
	# plt.ylabel("Amplitude")
	# plt.savefig('peaks.png')
	# plt.show()


	# print '\nf_lower (approx.):', 1.0/(ampPeaks['posTime'][-2]-ampPeaks['posTime'][-1])
	# print 'f_higher (approx.):', 1.0/(ampPeaks['posTime'][1] - ampPeaks['posTime'][2])
	# print ''




