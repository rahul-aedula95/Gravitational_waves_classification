""" Program to find the correlation between the strain data,
	and to design a GLM to predict data.
"""

from sklearn import linear_model as lm
import scipy
import numpy as np
from numpy import reshape
import matplotlib.pyplot as plt
import myconstants
import pylab
from pycbc.waveform import get_td_waveform
from math import sqrt
from math import log

def computeHPHC(argList):
	for apx in ['SEOBNRv2']:
		hp, hc = get_td_waveform(approximant=apx, \
			mass1 = argList[0], \
			mass2 = argList[1], \
			spin1z = argList[2], \
			delta_t = argList[3],\
			f_lower = argList[4])
	return hp, hc 

def readData(fileName, strict = False):
	""" 
		Assumes fileName as the input,
		returns a dictionary consisting of lists 
		whose key element is the name of the column,
		followed by list of values belonging to that column
	"""
	mat = dict()
	matHeader = fileName.readline()[:-2].split(",")
	for header in matHeader:
		mat[header.strip()] = list()
	for line in fileName:
		# Remove ",\n" from the line
		line = line[:-2]
		temp = line.split(",")
		for index, value in enumerate(temp):
			if strict == False:
				mat[matHeader[index].strip()].append(float(value.strip()))
			else:
				mat[matHeader[index].strip()].append(value.strip())
	return mat

def writeData(data, output):
	keys = data.keys()
	for i, key in enumerate(keys):
		if i == len(keys)-1:
			output.write(str(key)+",\n")
		else:
			output.write(str(key)+", ")

	for i in xrange(len(data[keys[0]])):
		for j, key in enumerate(keys):
			if j == len(keys)-1:
				output.write(str(data[key][i])+",\n")
			else:
				output.write(str(data[key][i])+", ")
	return

def exoToPyCBC(data):
	data['pMassSU'] = [float(t)*myconstants.EU_To_SU for t in data['pMassEU']]
	data['pFreqSeconds'] = [1/float(t)*myconstants.daysToSeconds for t in data['pPeriodDays']]
	data['GWFreq'] = [2*t for t in data['pFreqSeconds']]
	# Converting to float for the prediction
	data['sMassSU'] = [float(t) for t in data['sMassSU']]
	return data


def regression(amplitude, *args):
	"""
		Assumes samples of amplitude are provided along with arrays
		of factors. Returns a regression fit and the score.
	"""
	reg = lm.LinearRegression()
	y = np.array(amplitude).reshape(len(amplitude), 1)
	X = np.matrix(args).transpose()
	reg.fit(X, y)

	return reg, reg.score(X, y)

def predict(reg, *args):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	X = np.matrix(args).transpose()
	return list(reg.predict(X)[:, 0])

def anova(reg, y, *args):
	"""
		Assumes samples of amplitude are provided along with the
		regression model and the factors. Returns the anova table.
	"""
	anovaResults = dict()
	yMean = sum(y)/len(y)
	yHat = list()
	# Generate list of predicted values
	yHat = predict(reg, *args)
	# Compute Sum of Squared Residuals
	anovaResults['SST'] = sum([(t-yMean)**2 for t in y])
	# Compute Sum of Squares Total
	anovaResults['SSR'] = sum([(t-yMean)**2 for r in yHat])
	# Compute Sum of Sqaured Errors
	anovaResults['SSE'] = sum([(y[t]-yHat[t]) for t in xrange(len(y))])
	anovaResults['SSRegression'] = anovaResults['SSR']
	anovaResults['SSResidual'] = anovaResults['SSE']
	anovaResults['SSTotal'] = anovaResults['SST']
	anovaResults['RegressionFreedom'] = len(reg.coef_[0])
	anovaResults['ResidualFreedom'] = len(y) - anovaResults['RegressionFreedom'] 
	anovaResults['MSR'] = anovaResults['SSR']/anovaResults['RegressionFreedom']
	anovaResults['MSE'] = anovaResults['SSE']/anovaResults['ResidualFreedom']
	anovaResults['F'] = anovaResults['MSR']/anovaResults['MSE']
	anovaResults['R2'] = anovaResults['SSR']/anovaResults['SST']
	anovaResults['R'] = sqrt(anovaResults['R2'])
	print 'Hello:', anovaResults['SSR']+anovaResults['SSE']
	return anovaResults

def slope(y, x):
	""" 
		Assumes y and x are lists of values for two variables.
		Returns the slope (dy/dx) for the corresponding values.
	"""
	return [(y[t]-y[t-1])/(x[t]-x[t-1]) for t in xrange(1, len(y))]

def genSlopes(y, x, order):
	""" 
		Assumes y and x are lists of values for two variables.
		Returns the slopes - first order, second order and so on,
		as per the input order variable.
	"""
	slopes = dict()
	# m1 = dy/dx (first order)
	# m2 = second order differential and so on.
	m = slope(y, x)
	slopes['m1'] = m
	for i in xrange(2, order+1):
		xStart = len(x)-len(m)
		newSlope = slope(m, x[xStart:])
		index = 'm'+str(i)
		m = newSlope
		slopes[index] = newSlope
	return slopes

if __name__ == "__main__":
	# matData2.csv has m1 = (10, 20)
	inputData = open("./Data/testData.csv", "r")
	# strict=False let's you convert values to float (useful for computations)
	mat = readData(inputData, strict=False)
	mat['amplitude'] = [log(t*10**19) for t in mat['amplitude']]
	# print mat
	# regression coulld be called with variable number of factors
	reg, score = regression(mat['amplitude'], mat['m1'], mat['m2'])
	print 'Coefficients: '+str(reg.coef_)
	print 'Intercept: '+str(reg.intercept_)
	print 'Score: '+str(score)
	anovaResults = anova(reg, mat['amplitude'], mat['m1'], mat['m2'])
	# predict(reg, mat['m1'], mat['m2'])

	print '\nanovaResults:'
	print anovaResults

