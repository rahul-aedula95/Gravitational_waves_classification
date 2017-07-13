"""
	* Find the separation between coalescence and non-coalescence period using x = A.e^(B.y) and 
		reversing the results to obtain y from x.
	* Train the model x = A.y^B and obtain y from x for all points. Then use predicted values before 
		coalescence to plot before the separation.
	* Train the model x = A.y^B and obtain y from x for points after coalescence. Use these values
		to plot after coalescence.
"""

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
	keys = ['posTime', 'posAmp']
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
				output.write(str(data[key][i])+"	")
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

if __name__ == "__main__":
	delT = 1.0/4096
	hp, hc = computeHPHC([10, 2, 0.0, delT, 60])
	# hp, hc = computeHPHC([10, 10, 0.0, delT, 30])
	# hp, hc = computeHPHC([53, 44, 0.0, delT, 10])
	# hp, hc = computeHPHC([30, 30, 0.0, delT, 20])
	h = [t*10**19 for t in hp]
	ampPeaks = peaks(h, delT)

	outFile = open("./Data/ampPeaks.csv","w+")
	writeData(ampPeaks, outFile)
	outFile.close()

	y = np.array(ampPeaks['posAmp'])
	x = np.array(ampPeaks['posTime'])	

	popt, pcov = curve_fit(regFunc, y, x, [-1, -1, -1])
	a = popt[0]
	b = popt[1]
	c = popt[2]
	cleanTime = [t for t in ampPeaks['posTime'] if (t-c)/a > 0.0]
	cleanTimeLen = len(cleanTime)
	cleanTimeCopy = cleanTime

	y = np.array(ampPeaks['posAmp'][:len(cleanTime)])
	x = np.array(cleanTime)

	# NEW MODIFICATION
	popt, pcov = curve_fit(regFuncCoalesce, ampPeaks['posAmp'], ampPeaks['posTime'], [-1, -1])
	a = popt[0]
	b = popt[1]
	# cleanTime = [t for t in ampPeaks['posTime'] if (t-c)/a > 0.0]

	predAmp = predictCoalesce(cleanTime, a, b)

	# REPEATING THE PROCEDURE FOR NEGATIVE PEAKS
	y = np.array([-t for t in ampPeaks['negAmp']])
	x = np.array(ampPeaks['negTime'])

	popt, pcov = curve_fit(regFunc, y, x, [-1, -1, -1])
	a = popt[0]
	b = popt[1]
	c = popt[2]
	cleanTimeN = [t for t in ampPeaks['negTime'] if (t-c)/a > 0.0]
	cleanTimeNLen = len(cleanTimeN)

	y = np.array(ampPeaks['negAmp'][:len(cleanTimeN)])
	x = np.array(cleanTimeN)

	popt, pcov = curve_fit(regFuncCoalesce, np.array([-t for t in ampPeaks['negAmp']]), np.array(ampPeaks['negTime']), [-1, -1])
	a = popt[0]
	b = popt[1]

	predAmpN = [-t for t in predictCoalesce(cleanTimeN, a, b)]


	# Make them the same length
	if len(predAmp) > len(predAmpN):
		predAmp = predAmp[::-1][:len(predAmpN)]
		cleanTime = cleanTime[::-1][:len(predAmpN)]
		predAmpN = predAmpN[::-1]
		cleanTimeN = cleanTimeN[::-1]
	elif len(predAmp) < len(predAmpN):
		predAmpN = predAmpN[:len(predAmp)]
		cleanTimeN = cleanTimeN[:len(predAmpN)]
		predAmp = predAmp[::-1]
		cleanTime = cleanTime[::-1]

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

	# CHANGE INTRODUCED HERE
	temp = len(ampPeaks['posAmp']) - cleanTimeLen + 1
	tempN = len(ampPeaks['posAmp']) - cleanTimeLen

	y1 = ampPeaks['posAmp'][:temp][::-1][:-1]
	x1 = ampPeaks['posTime'][:temp][::-1][:-1]	

	tempDict = dict()
	tempDict['y'] = y1
	tempDict['x'] = x1

	outFile = open("./Data/ampPeaks1.csv","w+")
	writeData1(tempDict, outFile)
	outFile.close()

	y = np.array(y1)
	x = np.array(x1)	

	popt, pcov = curve_fit(regFuncCoalesce, y, x, [-1, -1])
	a = popt[0]
	b = popt[1]

	predAmpCoalesce = predictCoalesce(x1, a, b)
	predAmpCoalesce.append(max(ampPeaks['posAmp']))
	x1.append(max(ampPeaks['posTime']))

	# CHANGE INTRODUCED HERE (Negative Peaks)
	y2 = [-t for t in ampPeaks['negAmp'][:tempN][::-1]]
	x2 = ampPeaks['negTime'][:tempN][::-1]
	y = np.array(y2)
	x = np.array(x2)	

	popt, pcov = curve_fit(regFuncCoalesce, y, x, [-1, -1])
	a = popt[0]
	b = popt[1]

	predAmpNCoalesce = predictCoalesce(x2, a, b)
	predAmpNCoalesce = [-t for t in predAmpNCoalesce]

	# MODIFIED CODE
	timeCoalesce = list()
	ampCoalesce = list()
	length = max([len(x1), len(x2)])
	if abs(x1[0]) > abs(x2[0]):
		for t in xrange(length):
			if t < len(x1):
				timeCoalesce.append(x1[t])
				ampCoalesce.append(predAmpCoalesce[t])
			if t < len(x2):
				timeCoalesce.append(x2[t])
				ampCoalesce.append(predAmpNCoalesce[t])
	else:
		for t in xrange(length):
			if t < len(x2):
				timeCoalesce.append(x2[t])
				ampCoalesce.append(predAmpNCoalesce[t])
			if t < len(x1):
				timeCoalesce.append(x1[t])
				ampCoalesce.append(predAmpCoalesce[t])

	#=======================================================
	# MERGE the separation
	#=======================================================
	time = time + timeCoalesce
	amp = amp + ampCoalesce
	
	time = np.array(time)
	amp = np.array(amp)

	plt.figure(figsize=(10,8))
	plt.scatter(ampPeaks['posTime'], ampPeaks['posAmp'], color='g')
	plt.plot(hp.sample_times, h, color='b')
	# plt.scatter(ampPeaks['posTime'], predicted, color = 'r')
	# plt.scatter(cleanTime, predAmp, color = 'r')
	# plt.scatter(cleanTimeN, predAmpN, color ='r')
	# plt.scatter(x1, predAmpCoalesce, color = 'black')
	# plt.scatter(x2, predAmpNCoalesce, color='black')
	# plt.plot(xnew, splineFit, color = 'orange')
	plt.plot(time, amp, color = 'orange')
	# plt.plot(timeCoalesce, ampCoalesce, color='orange')
	plt.xlabel("Time")
	plt.ylabel("Amplitude")
	plt.savefig('peaks.png')
	plt.show()

	print '\nf_lower (approx.):', 1.0/(ampPeaks['posTime'][-2]-ampPeaks['posTime'][-1])
	print 'f_higher (approx.):', 1.0/(ampPeaks['posTime'][1] - ampPeaks['posTime'][2])



