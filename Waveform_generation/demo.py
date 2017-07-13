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
import time

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
	matHeader = fileName.readline()[:-1].split(",")
	print matHeader
	for header in matHeader:
		mat[header.strip()] = list()
	for line in fileName:
		# Remove ",\n" from the line
		line = line[:-1]
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

def writeDataSelected(data, output):
	keys = ['GWCoalesceAmp', 'sMassSU', 'pMassSU', 'GWFreq', 'sName', 'pName', 'pMassClass']
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

def exoToPyCBC(data):
	data['pMassSU'] = [float(t)*myconstants.EU_To_SU for t in data['pMassEU']]
	data['pFreqSeconds'] = [1.0/(float(t)*myconstants.daysToSeconds) for t in data['pPeriodDays']]
	data['GWFreq'] = [2*t for t in data['pFreqSeconds']]
	# Converting to float for the prediction
	data['sMassSU'] = [float(t) for t in data['sMassSU']]
	return data


def regression(amplitude, *args):
	"""
		Assumes samples of amplitude are provided along with arrays
		of factors. Returns a regression fit and the score.
	"""
	reg = lm.LinearRegression(fit_intercept = False)
	y = np.array(amplitude).reshape(len(amplitude), 1)
	X = np.matrix(args).transpose()
	reg.fit(X, y)
	return reg, reg.score(X, y)

# def regressionL(amplitude, *args):
# 	"""
# 		Assumes samples of amplitude are provided along with arrays
# 		of factors. Returns a regression (logistic) fit and the score.
# 	"""
# 	y = np.array(amplitude)
# 	reg = lm.LogisticRegression(C=1., solver='lbfgs')
# 	iris = load_iris()
# 	X = iris.data[:len(amplitude), :2]
# 	y = iris.target
# 	reg.fit(X, y)
# 	return X, y, reg, reg.score(X, y)

def predict(reg, *args):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	X = np.matrix(args).transpose()
	return list(reg.predict(X)[:, 0])

def predictL(reg, *args):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	X = np.matrix(args).transpose()
	# return list(reg.predict(X))
	return list(reg.predict(X)[:, 0])

def predictE(reg, *args):
	"""
		Assumes regression model is provided as input along with
		the predictors. Returns the predicted values as a list.
	"""
	X = np.matrix(args).transpose()
	pred = list(reg.predict(X)[:, 0])
	return [myconstants.e**t for t in pred]

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

def chripMass(m1, m2):
	prodMass = [(m1[t]*m2[t])**(3.0/5) for t in range(len(m1))]
	sumMass = [(m1[t]+m2[t])**(1.0/5) for t in range(len(m1))]
	return [(prodMass[t]/sumMass[t]) for t in range(len(prodMass))]

if __name__ == "__main__":
	"""
		Exo Planet computation disabled
	"""
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
	# testData Analysis and Scatter plot
	#===============================================================
	# testData.csv has m1 = (10, 99)
	# inputData = open("./Data/testData.csv", "r")
	# inputData = open("./Data/randMatData.csv", "r")
	inputData = open("./Data/datatest.csv", "r")
	# inputData = open("./Data/matData4.csv", "r")
	mat = readData(inputData, strict=False)
	# mat['f_low'] = [30 for t in mat['f_low']]

	mat['amplitude'] = [(t*10**19) for t in mat['amplitude']]
	mat['chirpMass'] = chripMass(mat['m1'], mat['m2'])
	reg, score = regression(mat['amplitude'], mat['chirpMass'])	
	print '==========================================================='
	print '                Regression Analysis:'
	print '==========================================================='
	print '\n             Simple Linear Model'
	print '            ---------------------'
	print 'Coefficients:', reg.coef_
	print 'Intercept:', reg.intercept_
	print 'score (coefficient of determination):', score
	test = int(raw_input('\npress 1 to show plot: '))
	if test == 1:
		plt.figure(figsize=(10,8))
		plt.scatter(mat['chirpMass'], mat['amplitude'], color='g')
		# print 'm1:', mat['m1']
		predAmplitude = predict(reg, mat['chirpMass'])
		# print 'predAmplitude:', predAmplitude
		# anovaResults = anova(reg, mat['amplitude'], mat['m11'])
		# print 'anova:', anovaResults

		plt.plot(mat['chirpMass'], predAmplitude)
		plt.xlabel("Mass chirpMass")
		plt.ylabel("Amplitude")
		plt.savefig('linearReg.png')
		plt.show()

		test = int(raw_input('\npress 1 to continue: '))
		if test == 1:
			mat['m12'] = [log(mat['m1'][t]) for t in xrange(len(mat['m1']))]
			mat['m22'] = [log(mat['m2'][t]) for t in xrange(len(mat['m2']))]
			mat['f_low2'] = [log(mat['f_low'][t]) for t in xrange(len(mat['f_low']))]
			mat['chirpMass2'] = [log(mat['chirpMass'][t]) for t in xrange(len(mat['chirpMass']))]
			regL, score = regression(mat['amplitude'], mat['chirpMass2'])

			print '\n             Logarithmic Regression'
			print '            ------------------------'
			print 'Coefficients:', regL.coef_
			print 'Intercept:', regL.intercept_
			print 'score (coefficient of determination):', score
			test = int(raw_input('\npress 1 to show plot: '))
			if test == 1:
				plt.figure(figsize=(10,8))
				plt.scatter(mat['m22'], mat['amplitude'], color='g')
				# print 'm1:', mat['m1']
				predAmplitude = predictL(regL, mat['chirpMass2'])
				# anovaResults = anova(reg, mat['amp'], mat['m12'])
				# print 'anova:', anovaResults

				plt.plot(mat['m22'], predAmplitude)
				plt.xlabel("Mass m2")
				plt.ylabel("Amplitude")
				plt.savefig('logReg.png')
				plt.show()

				test = int(raw_input('\npress 1 to continue: '))
				if test == 1:	
					mat['amp'] = [log(mat['amplitude'][t]) for t in xrange(len(mat['amplitude']))]
					regE, score = regression(mat['amp'], mat['m1'], mat['m2'], mat['f_low'])

					print '\n             Exponential Regression'
					print '            ------------------------'
					print 'Coefficients:', regE.coef_
					print 'Intercept:', regE.intercept_
					print 'score (coefficient of determination):', score
					test = int(raw_input('\npress 1 to show plot: '))
					if test == 1:
						plt.figure(figsize=(10,8))
						plt.scatter(mat['m2'], mat['amplitude'], color='g')
						# print 'm1:', mat['m1']
						predAmplitude = predictE(regE, mat['m1'], mat['m2'], mat['f_low'])
						# anovaResults = anova(reg, mat['amp'], mat['m12'])
						# print 'anova:', anovaResults

						plt.plot(mat['m2'], predAmplitude)
						plt.xlabel("Mass m2")
						plt.ylabel("Amplitude")
						plt.savefig('expoReg.png')
						plt.show()



	print 'Computation Start (Predict Exo Planet Amplitude)\n'
	inputExoData = open("./Data/cleanExoData1.csv", "r")
	""" 
		Convert the exo planets data compatible to the
		PyCBC Data. Convert Earth-units to Sun-units. Introduce GW Frequency.
	"""
	exoData = exoToPyCBC(readData(inputExoData, strict=True))

	start = time.time()
	# Predict the GW amplitude for the exo planets using the linear model (trained)
	# exoData['sMassSU1'] = [log(t) for t in exoData['sMassSU']]
	# exoData['pMassSU1'] = [log(t) for t in exoData['pMassSU']]
	# exoData['GWFreq1'] = [log(t) for t in exoData['GWFreq']]
	exoData['chirpMass'] = chripMass(exoData['sMassSU'], exoData['pMassSU'])
	exoData['GWCoalesceAmp'] = predict(reg, exoData['chirpMass'])
	exoData['GWCoalesceAmp'] = [t*(10**-19) for t in exoData['GWCoalesceAmp']]
	print regE.coef_
	print regE.intercept_
	# print exoData['GWCoalesceAmp']
	end = time.time() - start
	# Write the predicted amplitude to a .csv file
	outputExoPredicted = open("./Data/ExoDataPredicted.csv", "w+")
	writeDataSelected(exoData, outputExoPredicted)
	outputExoPredicted.close()
	print 'Computed Successfully'
	print 'Computation time:', end
	print '\n\n'

	ch = chripMass([87	], [1.99])
	test = [ch]
	print predict(reg, ch)










	
