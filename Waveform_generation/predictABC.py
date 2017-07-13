"""
	Code to predict the coefficients of the Fitted Amplitudes for Non-Coalescence phase
	and the Coalescence phase.
"""

from sklearn import linear_model as lm
import scipy
import numpy as np
from numpy import reshape
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import myconstants
import pylab
from pycbc.waveform import get_td_waveform
from math import log
from math import sqrt
from sklearn.datasets import load_iris
from sklearn.preprocessing import PolynomialFeatures

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
	keys = ["xB", "chirpMass", "f"]
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

def chirpMass(m1, m2):
	prodMass = [(m1[t]*m2[t])**(3.0/5) for t in range(len(m1))]
	sumMass = [(m1[t]+m2[t])**(1.0/5) for t in range(len(m1))]
	return [(prodMass[t]/sumMass[t]) for t in range(len(prodMass))]

if __name__ == "__main__":
	inFile = open("./Data/WaveCharacteristics.csv", "r")
	outFile = open("./Data/xB.csv", "w+")
	data = readData(inFile)
	data["chirpMass"] = chirpMass(data["m1"], data["m2"])
	reg, score = regression(data["xB"], data["chirpMass"], data["f"])
	writeData(data, outFile)
	outFile.close()
	inFile.close()

	# plt.figure(figsize=(10,8))

	fig = plt.figure(figsize=(10,8))
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(data["chirpMass"], data["f"], data["xB"])
	ax.set_xlabel("chirpMass")
	ax.set_ylabel("f")
	ax.set_zlabel("xA")
	plt.show()

	print score
