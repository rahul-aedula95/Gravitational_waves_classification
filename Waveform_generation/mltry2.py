import numpy as np 
import csv
from sklearn import linear_model as lm
filename = './Data/datatest.csv'
from numpy import genfromtxt
from sklearn.preprocessing import PolynomialFeatures
from numpy import reshape
import matplotlib.pyplot as plt
import myconstants
import time

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
	keys = ['GWCoalesceAmp', 'sMassSU', 'pMassSU', 'GWFreq']
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

my_data = genfromtxt(filename, delimiter=',',comments='#')

row,column = np.shape(my_data)
#b = pow(10,19)
col = 2
x = np.zeros(shape=(row-1, col))

poly = PolynomialFeatures(degree=2)

for i in range(1, row):
	# Considering only m1 and m2 while testing
	for j in range(0, col):
		x[i-1][j]=my_data[i][j+1]

# print my_data
print x        

y = np.zeros(shape=(1,row-1))


x1 = poly.fit_transform(x)
for k in range(1,row):
    y[0][k-1] = my_data[k][0]
    
# print x1
# print y
y = map(lambda x: x * 10**19, y)
y = np.array(y)

'''
for k in range(0,row):
    x[k] = my_data[k][2]
'''

# print np.shape(x1)
# print np.shape(y)
y = y.transpose()
reg = lm.LinearRegression(fit_intercept = False)
reg.fit(x1, y)


# print reg.coef_

# y = np.array(y).reshape(len(y), 1)
#pred = reg.predict([7,8,0.0,0.00006103515625,20])
# score1= reg.score(x,y)

#print pred   

# print score1

print '\nreg.coef_:', reg.coef_

print '\nreg.intercept_:', reg.intercept_

print '\nreg.score_:', reg.score(x1, y)

print ''

inputData = open("./Data/datatest.csv", "r")
mat = readData(inputData, strict=False)
mat['amplitude'] = [t*10**19 for t in mat['amplitude']]

plt.figure(figsize=(10,8))
plt.scatter(mat['m2'], mat['amplitude'], color='g')
		# print 'm1:', mat['m1']
predAmplitude = predictP(reg, mat['m1'], mat['m2'])
		# print 'predAmplitude:', predAmplitude
		# anovaResults = anova(reg, mat['amplitude'], mat['m11']
		# print 'anova:', anovaResults

plt.plot(mat['m2'], predAmplitude)
plt.xlabel("Mass m2")
plt.ylabel("Amplitude")
plt.savefig('linearReg.png')
plt.show()

print 'Computation Start (Predict Exo Planet Amplitude)\n'
inputExoData = open("./Data/cleanExoData.csv", "r")
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
# exoData['chirpMass'] = chripMass(exoData['sMassSU'], exoData['pMassSU'])
exoData['GWCoalesceAmp'] = predictP(reg, exoData['sMassSU'], exoData['pMassSU'])
exoData['GWCoalesceAmp'] = [t*(10**-19) for t in exoData['GWCoalesceAmp']]

# print exoData['GWCoalesceAmp']
end = time.time() - start
# Write the predicted amplitude to a .csv file
outputExoPredicted = open("./Data/ExoDataPredicted.csv", "w+")
writeDataSelected(exoData, outputExoPredicted)
outputExoPredicted.close()
print 'Computed Successfully'
print 'Computation time:', end
print '\n\n'