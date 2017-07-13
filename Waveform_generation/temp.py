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

def regFuncCoalesce(x, a, b, c, d):
	# return a*np.exp(b*x)+c
	return 1/(a*x**2 + b*x +c) + d

def predictCoalesce(x, a, b, c, d):
	# return [a*np.exp(b*t)+c for t in x]
	return [(1/(a*t**2 + b*t +c) + d) for t in x]

if __name__ == "__main__":
	delT = 1.0/4096
	hp, hc = computeHPHC([10, 10, 0.0, delT, 60])
	h = [t*10**19 for t in hp]
	ampPeaks = peaks(h, delT)

	y = ampPeaks['posAmp'][::-1]
	x = ampPeaks['posTime'][::-1]
	slopeChangeIndex = 0
	for i in xrange(len(y)):
		if (y[i+1]-y[i])/(x[i+1]-x[i]) >= 1.0:
			slopeChangeIndex = i
			break
	slopeChangeIndex = len(y)-slopeChangeIndex


	y1 = ampPeaks['posAmp'][:slopeChangeIndex][::-1]
	x1 = ampPeaks['posTime'][:slopeChangeIndex][::-1]
	y = np.array(y1)
	x = np.array(x1)	

	popt, pcov = curve_fit(regFuncCoalesce, x1, y1, [1, 1, 1, 1])
	print popt
	a = popt[0]
	b = popt[1]
	c = popt[2]
	d = popt[3]

	predAmpCoalesce = predictCoalesce(x1, a, b, c, d)

	y2 = [-t for t in ampPeaks['negAmp'][:slopeChangeIndex][::-1]]
	x2 = ampPeaks['negTime'][:slopeChangeIndex][::-1]
	y = np.array(y2)
	x = np.array(x2)	

	popt, pcov = curve_fit(regFuncCoalesce, x, y, [1, 1, 1, 1])
	print popt
	a = popt[0]
	b = popt[1]
	c = popt[2]
	d = popt[3]

	predAmpNCoalesce = predictCoalesce(x2, a, b, c, d)


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

	# timeCoalesce = timeCoalesce[::-1]
	# ampCoalesce = ampCoalesce[::-1]


	plt.figure(figsize=(10,8))
	plt.scatter(ampPeaks['posTime'], ampPeaks['posAmp'], color='g')
	plt.plot(hp.sample_times, h, color='b')
	# plt.scatter(ampPeaks['posTime'], predicted, color = 'r')
	plt.scatter(x1, predAmpCoalesce, color = 'r')
	plt.scatter(x2, predAmpNCoalesce, color='r')
	# plt.scatter(cleanTimeN, predAmpN, color ='black')
	# plt.plot(xnew, splineFit, color = 'orange')
	plt.plot(timeCoalesce, ampCoalesce, color = 'orange')
	plt.xlabel("Time")
	plt.ylabel("Amplitude")
	# plt.savefig('peaks.png')
	plt.show()