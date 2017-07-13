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
from math import log
from math import sqrt
from sklearn.datasets import load_iris
from sklearn.preprocessing import PolynomialFeatures

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
	data['pFreqHertz'] = [1/(float(t)*myconstants.daysToSeconds) for t in data['pPeriodDays']]
	data['GWFreq'] = [2*t for t in data['pFreqHertz']]
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
	# y = list()
	# for i in xrange(len(args)):
	# 	y.append(amplitude)
	# print y
	# y = np.matrix(y).transpose()
	return reg, reg.score(X, y)

def regressionL(amplitude, *args):
	"""
		Assumes samples of amplitude are provided along with arrays
		of factors. Returns a regression (logistic) fit and the score.
	"""
	y = np.array(amplitude)
	reg = lm.LogisticRegression(C=1., solver='lbfgs')
	# y = np.array(amplitude).reshape(len(amplitude), 1)
	# X = np.matrix(args).transpose()
	iris = load_iris()
	X = iris.data[:len(amplitude), :2]
	y = iris.target
	reg.fit(X, y)
	print X
	print np.matrix(args).transpose()
	# y = list()
	# for i in xrange(len(args)):
	# 	y.append(amplitude)
	# print y
	# y = np.matrix(y).transpose()
	return X, y, reg, reg.score(X, y)

# def regressionTest(amplitude, *args):
# 	"""
# 		Assumes samples of amplitude are provided along with arrays
# 		of factors. Returns a regression fit and the score.
# 	"""
# 	reg = lm.LinearRegression()
# 	print args
# 	y = list()
# 	y.append(amplitude)
# 	# y = np.array(amplitude).reshape(len(amplitude), 1)
# 	y = np.matrix(y).transpose()
# 	X = np.matrix(args[0]).transpose()
# 	reg.fit(X, y)

# 	return reg, reg.score(X, y)

def predict(reg, *args):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	X = np.matrix(args).transpose()
	return list(reg.predict(X)[:, 0])

def predictL(reg, X):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	# X = np.matrix(args).transpose()
	return list(reg.predict(X))

def predictE(reg, *args):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	X = np.matrix(args).transpose()
	pred = list(reg.predict(X)[:, 0])
	return [myconstants.e**(-t) for t in pred]

def predictP(reg, *args):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	x = np.zeros(shape=(len(args[0]), len(args)))
	for i in xrange(len(args[0])):
		for j in xrange(len(args)):
			x[i][j] = args[j][i]
	x1 = poly.fit_transform(x)		

	x2 = np.zeros(shape=(len(x1[0]), len(x1)))
	for i in xrange(len(x1[0])):
		for j in xrange(len(x1)):
			x2[i][j] = x1[j][i]
	X = np.matrix(x2).transpose()
	return list(reg.predict(X)[:, 0])

def anova(reg, y, *args):
	"""
		Assumes samples of amplitude are provided along with the
		regression model and the factors. Returns the annova table.
	"""
	x = args[0]
	xMean = sum(x)/len(x)
	sX = [t-xMean for t in x]
	sXX = sum([t**2 for t in sX])
	yMean = sum(y)/len(y)
	sY = [t-yMean for t in y]
	yHat = list()
	anovaResults = dict()
	# Generate list of predicted values
	yHat = predict(reg, *args)
	# print 'yHat:', yHat
	# print 'y:', y
	anovaResults['SSR'] = (sum([(sX[t]*sY[t]) for t in xrange(len(sX))])**2)/sXX
	anovaResults['SST'] = sum([t**2 for t in y])
	anovaResults['SSE'] = anovaResults['SST'] - anovaResults['SSR']
	anovaResults['R2'] = anovaResults['SSR']/anovaResults['SST']
	anovaResults['R'] = sqrt(anovaResults['R2'])
	# Compute Sum of Squared Residuals
	# anovaResults['SST'] = sum([(t-yMean)**2 for t in y])
	# anovaResults['SST'] = len(y)*(yMean**2)
	# # Compute Sum of Squares Total
	# anovaResults['SSR'] = sum([(t-yMean)**2 for r in yHat])
	# Compute Sum of Sqaured Errors
	# anovaResults['SSE'] = sum([(y[t]-yHat[t]) for t in xrange(len(y))])
	# anovaResults['SSRegression'] = anovaResults['SSR']
	# anovaResults['SSResidual'] = anovaResults['SSE']
	# anovaResults['SSTotal'] = anovaResults['SST']
	# anovaResults['RegressionFreedom'] = len(reg.coef_[0])
	# anovaResults['ResidualFreedom'] = len(y) - anovaResults['RegressionFreedom'] 
	# anovaResults['MSR'] = anovaResults['SSR']/anovaResults['RegressionFreedom']
	# anovaResults['MSE'] = anovaResults['SSE']/anovaResults['ResidualFreedom']
	# anovaResults['F'] = anovaResults['MSR']/anovaResults['MSE']
	# anovaResults['R2'] = anovaResults['SSR']/anovaResults['SST']
	# anovaResults['R'] = sqrt(anovaResults['R2'])
	return anovaResults

def slope(y, x):
	""" 
		Assumes y and x are lists of values for two variables.
		Returns the slope (dy/dx) for the corresponding values.
	"""
	return [(y[t]-y[t-1])/(x[t]-x[t-1]) for t in xrange(1, len(y))]
	# return [math.abs((y[t]-y[t-1]))/(x[t]-x[t-1]) for t in xrange(1, len(y))]

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
	# # matData2.csv has m1 = (10, 20)
	# inputData = open("./Data/matData2.csv", "r")
	# # strict=False let's you convert values to float (useful for computations)
	# mat = readData(inputData, strict=False)
	# print mat
	# # regression coulld be called with variable number of factors
	# reg, score = regression(mat['amplitude'], mat['m1'], mat['m2'])
	# print 'Coefficients: '+str(reg.coef_)
	# print 'Intercept: '+str(reg.intercept_)
	# print 'Score: '+str(score)
	# anovaResults = anova(reg, mat['amplitude'], mat['m1'], mat['m2'])
	# # predict(reg, mat['m1'], mat['m2'])

	# inputExoData = open("./Data/cleanExoData.csv", "r")
	# """ 
	# 	Convert the exo planets data compatible to the
	# 	PyCBC Data. Convert Earth-units to Sun-units. Introduce GW Frequency.
	# """
	# exoData = exoToPyCBC(readData(inputExoData, strict=True))

	# # Predict the GW amplitude for the exo planets using the linear model (trained)
	# exoData['GWCoalesceAmp'] = predict(reg, exoData['sMassSU'], exoData['pMassSU'])
	# print "The predicted GW Amplitude for exoplanets is: "
	# # print exoData['GWCoalesceAmp']
	# print exoData['GWCoalesceAmp']

	# # Write the predicted amplitude to a .csv file
	# outputExoPredicted = open("./Data/ExoDataPredicted.csv", "w+")
	# writeData(exoData, outputExoPredicted)
	# outputExoPredicted.close()

	#===============================================================
	# matDataNeo Analysis and Scatter plot
	#===============================================================
	# matData2.csv has m1 = (10, 20)
	inputData = open("./Data/testData.csv", "r")
	# strict=False let's you convert values to float (useful for computations)
	mat = readData(inputData, strict=False)
	# print mat
	# regression coulld be called with variable number of factors
	mat['amp'] = list()
	mat['m11'] = list()
	mat['amplitude'] = [(t*10**19) for t in mat['amplitude']]
	for i, num in enumerate(mat['m1']):
		if mat['m2'][i] == 10:
			mat['m11'].append(num)
			mat['amp'].append(mat['amplitude'][i])

	mat['m11'] = [(mat['m11'][t]) for t in xrange(81)]
	mat['amp'] = [(mat['amp'][t]) for t in xrange(81)]
	# mat['m2'] = [(mat['m2'][t]) for t in xrange(len(mat['m2']))]
	reg, score = regression(mat['amp'], mat['m11'])	
	print '==========================================================='
	print '                Regression Analysis:'
	print '==========================================================='
	print '\n             Simple Linear Model'
	print '             ----------------------'
	print 'Coefficients:', reg.coef_
	print 'Intercept:', reg.intercept_
	print 'score (coefficient of determination):', score
	print 'correlation:', sqrt(score)
	# test = int(raw_input('\npress 1 to show plot: '))
	# if test == 1:
		# plt.figure(figsize=(10,8))
		# plt.scatter(mat['m11'], mat['amp'], color='g')
		# print 'm1:', mat['m1']
		# predAmplitude = predict(reg, mat['m11'])
		# anovaResults = anova(reg, mat['amplitude'], mat['m11'])
		# print 'anova:', anovaResults


	#===============================================================
	# Analysis of slopes
	#===============================================================
	
	delT = mat['del_t'][0]
	time = +delT
	ampPeaks = dict()
	ampPeaks['amp'] = list()
	ampPeaks['time'] = list()
	ampPeaks['posAmp'] = list()
	ampPeaks['negAmp'] = list()
	ampPeaks['posTime'] = list()
	ampPeaks['negTime'] = list()

	hp, hc = computeHPHC([10, 10, 0.0, 1.0/4096, 60])
	hp = [t*10**19 for t in hp]
	# print hp
	hp = list(hp)
	# print hp
	initValue = max(hp)
	# print 'Hello:', initValue
	coalesceIndex = hp.index(initValue)
	# print initValue
	print 'coalesceIndex:', coalesceIndex


	# True = amplitude (ascending)
	# False = amplitude (descending)
	Toggle = True

	for value in hp[coalesceIndex::-1]:
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

	# Generate Slopes
	m = genSlopes(ampPeaks['posAmp'], ampPeaks['posTime'], 1)
	# print 'm:', m

	ampPeaksExpo = [log(1.0/t) for t in ampPeaks['posAmp']]
	regSlope, score = regression(ampPeaksExpo, ampPeaks['posTime'])
	print 'coef:', regSlope.coef_
	print 'intercept:', regSlope.intercept_
	print ampPeaks['posAmp']


	test = int(raw_input('\npress 1 to show plot: '))
	if test == 1:
		plt.figure(figsize=(10,8))
		plt.scatter(ampPeaks['posTime'], ampPeaks['posAmp'], color='g')
		# print 'm1:', mat['m1']
		predAmplitude = predictE(regSlope, ampPeaks['posTime'])
		# anovaResults = anova(reg, mat['amplitude'], mat['m11'])
		# print 'anova:', anovaResults

		plt.plot(ampPeaks['posTime'], predAmplitude)
		plt.xlabel("posTime")
		plt.ylabel("Amplitude")
		plt.savefig('slopeAnalysis.png')
		plt.show()

	filePath = "./Data/waveData.csv"
	output = open(filePath, "w+")
	for i in xrange(len(ampPeaks['posTime'])):
		output.write(str(ampPeaks['posTime'][i])+", "+str(ampPeaks['posAmp'][i])+"\n")
	output.close()

	# Code 
	row = len(ampPeaks['posTime'])
	print row, len(ampPeaks['posAmp'])
	col = 1
	x = np.zeros(shape=(row, col))

	poly = PolynomialFeatures(degree=2)

	for i in range(row):
		x[i][0] = ampPeaks['posTime'][i]

	# print my_data
	print x        

	y = np.zeros(shape=(1,row))


	x1 = poly.fit_transform(x)
	for k in range(row):
	    y[0][k] = ampPeaks['posAmp'][k]
	    
	# print x1
	# print y
	# y = map(lambda x: x * 10**19, y)
	y = np.array(y)

	'''
	for k in range(0,row):
	    x[k] = my_data[k][2]
	'''

	# print np.shape(x1)
	# print np.shape(y)
	y = y.transpose()
	reg = lm.LinearRegression()
	reg.fit(x1, y)

	print '\nreg.coef_:', reg.coef_

	print '\nreg.intercept_:', reg.intercept_

	print '\nreg.score_:', reg.score(x1, y)










	
