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

def computeHPHC(argList):
	for apx in ['SEOBNRv2']:
		hp, hc = get_td_waveform(approximant=apx, \
			mass1 = argList[0], \
			mass2 = argList[1], \
			spin1z = argList[2], \
			delta_t = argList[3],\
			f_lower = argList[4])
	return hp, hc 

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

if __name__ == "__main__":
	delT = 1.0/4096
	hp, hc = computeHPHC([10, 10, 0.0, delT, 60])
	h = [t*10**19 for t in hp]
	ampPeaks = peaks(h, delT)

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


	time = list()
	amp = list()
	length = max([len(cleanTime), len(cleanTimeN)])
	for t in xrange(length):
		if t < len(cleanTime):
			time.append(cleanTime[t])
			amp.append(predAmp[t])
		if t < len(cleanTimeN):
			time.append(cleanTimeN[t])
			amp.append(predAmpN[t])

	ampPeaks1 = copy.deepcopy(ampPeaks)
	ampPeaks1['posAmp'] = ampPeaks1['posAmp'][0:nonCoalescePeriod]
	ampPeaks1['posTime'] = ampPeaks1['posTime'][0:nonCoalescePeriod]
	outFile = open("./Data/ampPeaks.csv","w+")
	writeData(ampPeaks1, outFile)
	

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
	plt.plot(hp.sample_times, h, color='b')
	# plt.scatter(ampPeaks['posTime'], predicted, color = 'r')
	plt.scatter(cleanTime, predAmp, color = 'r')
	plt.scatter(cleanTimeN, predAmpN, color ='black')
	# plt.plot(xnew, splineFit, color = 'orange')
	plt.plot(time, amp, color = 'orange')
	plt.xlabel("Time")
	plt.ylabel("Amplitude")
	plt.savefig('peaks.png')
	plt.show()


	print '\nf_lower (approx.):', 1.0/(ampPeaks['posTime'][-2]-ampPeaks['posTime'][-1])
	print 'f_higher (approx.):', 1.0/(ampPeaks['posTime'][1] - ampPeaks['posTime'][2])



