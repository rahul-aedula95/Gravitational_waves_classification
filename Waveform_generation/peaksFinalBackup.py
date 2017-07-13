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
	print 'delT:', delT
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

def regressionPeaks(amp, time):
	y = [log(t) for t in amp]
	x = time
	# computing coefficient b
	n = len(y)
	sigmaXY = sum([y[t]*x[t] for t in xrange(n)])
	sigmaX = sum(x)
	sigmaY = sum(y)
	sigmaX2 = sum([x[t]**2 for t in xrange(n)])
	b = (sigmaXY - (1.0/n)*(sigmaX)*(sigmaY))/(sigmaX2-(1.0/n)*(sigmaX**2))
	# computing intercept a
	a = (1.0/n)*(sigmaY - b*sigmaX)
	# Computing A and B in exponential function.
	A = exp(a)
	B = b
	return A, B

def regFunc(x, a, b, c):
	return a*np.exp(b*x)+c
	# return (1.0/b)*np.log((x-c)/a)

def predict(x, a, b, c):
	return [(1.0/b)*log((t-c)/a) for t in x]

# def regFuncCoalesce(x, a, b, c):
# 	return a*np.exp(b*x)+c

# def predictCoalesce(x, a, b, c):
# 	return [a*np.exp(b*t)+c for t in x]

# def regFuncCoalesce(x, a, b, c, d):
# 	# return a*np.exp(b*x)+c
# 	return 1/(a*x**2 + b*x +c) + d

# def predictCoalesce(x, a, b, c, d):
# 	# return [a*np.exp(b*t)+c for t in x]
# 	return [(1/(a*t**2 + b*t +c) + d) for t in x]

def regFuncCoalesce(x, a, b):
	# return a*np.exp(b*x)+c
	# return a*np.exp(b*x)
	return a*(x**b)

def predictCoalesce(x, a, b):
	# return [a*np.exp(b*t)+c for t in x]
	# return [(1.0/b)*np.log(t/a) for t in x]
	return [(t/a)**(1.0/b) for t in x]

if __name__ == "__main__":
	delT = 1.0/4096
	hp, hc = computeHPHC([10, 2, 0.0, delT, 60])
	# hp, hc = computeHPHC([10, 2, 0.0, delT, 40])
	# hp, hc = computeHPHC([10, 10, 0.0, delT, 30])
	# hp, hc = computeHPHC([53, 44, 0.0, delT, 10])
	h = [t*10**19 for t in hp]
	ampPeaks = peaks(h, delT)
	print 'iN:', h.index(min(h))
	print 'i:', h.index(max(h))

	A, B = regressionPeaks(ampPeaks['posAmp'], ampPeaks['posTime'])
	predicted = [A*exp(B*t) for t in ampPeaks['posTime']]

	# plt.figure(figsize=(10,8))
	# plt.scatter(ampPeaks['posTime'], ampPeaks['posAmp'], color='g')
	# plt.plot(hp.sample_times, h, color='b')
	# plt.scatter(ampPeaks['posTime'], predicted, color = 'r')
	# plt.xlabel("Time")
	# plt.ylabel("Amplitude")
	# plt.savefig('peaks.png')
	# plt.show()

	outFile = open("./Data/ampPeaks.csv","w+")
	writeData(ampPeaks, outFile)
	outFile.close()

	# Symmetrical Sigmoidal
	# ampPeaks['posTime'] = [t-min(ampPeaks['posTime']) for t in ampPeaks['posTime']]
	y = np.array(ampPeaks['posAmp'])
	x = np.array(ampPeaks['posTime'])	

	popt, pcov = curve_fit(regFunc, y, x, [-1, -1, -1])
	print popt
	a = popt[0]
	b = popt[1]
	c = popt[2]
	cleanTime = [t for t in ampPeaks['posTime'] if (t-c)/a > 0.0]

	predAmp = predict(cleanTime, a, b, c)
	print predAmp

	# CHECK THIS
	# if max(ampPeaks['posTime']) not in cleanTime:
	# 	cleanTime.insert(0, max(ampPeaks['posTime']))
	# 	predAmp.insert(0, max(ampPeaks['posAmp']))

	# Predicting negative amplitude valuese()
	y = np.array([-t for t in ampPeaks['negAmp']])
	x = np.array(ampPeaks['negTime'])

	popt, pcov = curve_fit(regFunc, y, x, [-1, -1, -1])
	print popt
	a = popt[0]
	b = popt[1]
	c = popt[2]
	cleanTimeN = [t for t in ampPeaks['negTime'] if (t-c)/a > 0.0]

	predAmpN = [-t for t in predict(cleanTimeN, a, b, c)]
	print predAmpN
	# CHECK THIS
	# if max(ampPeaks['negTime']) not in cleanTime:
	# 	cleanTimeN.insert(0, max(ampPeaks['negTime']))
	# 	predAmpN.insert(0, min(ampPeaks['negAmp']))


	# Separating the non-coalescence Period.
	numEleCleaned = len(ampPeaks['posAmp'])-len(predAmp)
	init = 0
	while ampPeaks['posAmp'][numEleCleaned+init] < predAmp[init]:
		init+=1
	# Marks the start of the nonCoalesce phase for poAmp
	nonCoalescePeriod = numEleCleaned+init
	predAmp = predAmp[init::][::-1]
	cleanTime = cleanTime[init::][::-1]
	print predAmp, cleanTime

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
	print nonCoalescePeriod
	print 'len:', nonCoalescePeriod

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
		print predAmp
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
	print 's:', slopeChangeIndex
	print 'l:', len(predAmp)
	slopeChangeIndex = len(y)-slopeChangeIndex

	slopeChangeIndex = min(slopeChangeIndex, nonCoalescePeriod)
	slopeChangeIndex = max(slopeChangeIndex, 4)

	# print 'comp:', slopeChangeIndex, nonCoalescePeriod
	# slopeChangeIndex = nonCoalescePeriod

	# CHANGE INTRODUCED HERE
	temp = len(ampPeaks['posAmp']) - len(cleanTime)
	# y1 = ampPeaks['posAmp'][:slopeChangeIndex][::-1][:-1]
	# x1 = ampPeaks['posTime'][:slopeChangeIndex][::-1][:-1]
	y1 = ampPeaks['posAmp'][:temp][::-1]
	x1 = ampPeaks['posTime'][:temp][::-1]
	print 'y1:', y1
	print 'x1:', x1

	tempDict = dict()
	tempDict['y'] = y1
	tempDict['x'] = x1
	# tempDict['x'] = [t+abs(min(x1)) for t in x1]

	outFile = open("./Data/ampPeaks1.csv","w+")
	writeData1(tempDict, outFile)
	outFile.close()

	y = np.array(y1)
	x = np.array(x1)	

	popt, pcov = curve_fit(regFuncCoalesce, y, x, [-1, -1])
	print popt
	a = popt[0]
	b = popt[1]
	# c = popt[2]
	# d = popt[3]

	x1 = x1[:-1]
	y1 = y1[:-1]

	predAmpCoalesce = predictCoalesce(x1, a, b)
	predAmpCoalesce.append(max(ampPeaks['posAmp']))
	x1.append(max(ampPeaks['posTime']))
	print predAmpCoalesce

	# CHANGE INTRODUCED HERE
	# y2 = [-t for t in ampPeaks['negAmp'][:slopeChangeIndex][::-1]]
	# x2 = ampPeaks['negTime'][:slopeChangeIndex][::-1]
	# ampPeaks['negAmp'].insert(0, min(hp))
	# ampPeaks['negTime'].insert(0, 0)
	y2 = [-t for t in ampPeaks['negAmp'][:temp][::-1]]
	x2 = ampPeaks['negTime'][:temp][::-1]
	y = np.array(y2)
	x = np.array(x2)	

	popt, pcov = curve_fit(regFuncCoalesce, y, x, [-1, -1])
	print popt
	a = popt[0]
	b = popt[1]
	# c = popt[2]
	# d = popt[3]

	# x2 = x2[:-1]
	# y2 = y2[:-1]
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

	# ================================================

	#===================================================
	# Merging NonCoalescence with Coalescence
	#===================================================

	# time = time + timeCoalesce
	# amp = amp + ampCoalesce

	#===================================================


	# ampPeaks1 = copy.deepcopy(ampPeaks)
	# ampPeaks1['posAmp'] = ampPeaks1['posAmp'][0:nonCoalescePeriod]
	# ampPeaks1['posTime'] = ampPeaks1['posTime'][0:nonCoalescePeriod]
	# outFile = open("./Data/ampPeaks.csv","w+")
	# writeData(ampPeaks1, outFile)
	

	# Cubic Spline Fitting
	time = np.array(time)
	amp = np.array(amp)
	# amp1 = [y for (x, y) in sorted(zip(time, amp))][::-1]
	# time1 = sorted(time, reverse=True)
	# time1 = np.array(time1)
	# amp1 = np.array(amp1)

	print len(predAmp), len(predAmpN)
	# print 'time:', time
	# print 'amp:', amp
	xnew = np.linspace(time.min(), time.max(), 1000)																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																																							
	# splineFit = spline(time, amp, xnew)

	fit = interp1d(time, amp, kind='cubic')
	splineFit = fit(xnew)
	# fit2 = BarycentricInterpolator(time, amp)
	# fit1 = KroghInterpolator(time, amp)
	# splineFit1 = fit1(xnew)
	# splineFit2 = fit2(xnew)
	# fit3 = interpolate.InterpolatedUnivariateSpline(time, amp)
	# splineFit3 = fit3(xnew)
	# fit4 = interpolate.splrep(time[::-1], amp[::-1], s=10000)
	# splineFit4 = interpolate.splev(xnew[::-1], fit4, der=3)

	plt.figure(figsize=(10,8))
	plt.scatter(ampPeaks['posTime'], ampPeaks['posAmp'], color='g')
	plt.scatter(ampPeaks['negTime'], ampPeaks['negAmp'], color='g')
	plt.plot(hp.sample_times, h, color='b')
	# plt.scatter(ampPeaks['posTime'], predicted, color = 'r')
	plt.scatter(cleanTime, predAmp, color = 'r')
	plt.scatter(cleanTimeN, predAmpN, color ='r')
	plt.scatter(x1, predAmpCoalesce, color = 'black')
	plt.scatter(x2, predAmpNCoalesce, color='black')
	# plt.plot(xnew, splineFit, color = 'orange')
	plt.plot(time, amp, color = 'orange')
	plt.plot(timeCoalesce, ampCoalesce, color='orange')
	plt.xlabel("Time")
	plt.ylabel("Amplitude")
	plt.savefig('peaks.png')
	plt.show()


	print '\nf_lower (approx.):', 1.0/(ampPeaks['posTime'][-2]-ampPeaks['posTime'][-1])
	print 'f_higher (approx.):', 1.0/(ampPeaks['posTime'][1] - ampPeaks['posTime'][2])




