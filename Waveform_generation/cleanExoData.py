""" 
	Program to find clean the exo-planet data.
	The code removes insufficient information from the data file.
"""

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

def exoDataClean(data, output):
	"""
		Assumes a dictionary (exo-planet) as the input,
		returns clean data without insufficient information
		for the GW Prediction.
		Also, the output filehandle is input and data is written to
		it as a .csv
	"""
	keys = data.keys()
	cleanData = dict()
	for i, key in enumerate(keys):
		if i == len(keys)-1:
			output.write(str(key)+",\n")
		else:
			output.write(str(key)+", ")
		cleanData[key] = list()

	for i in xrange(len(data[keys[0]])):
		flag = True
		for key in keys:
			if data[key][i].strip() == "":
				flag = False
				break
		if flag == True:
			for j, key in enumerate(keys):
				if j == len(keys)-1:
					output.write(str(data[key][i])+",\n")
				else:
					output.write(str(data[key][i])+", ")
				cleanData[key].append(data[key][i])
	return cleanData

if __name__ == "__main__":
	inputExoData = open("./Data/exoData1.csv", "r")
	outputExoData = open("./Data/cleanExoData1.csv", "w+")
	exoData = exoDataClean(readData(inputExoData, strict=True), outputExoData)
	outputExoData.close()
	inputExoData.close()