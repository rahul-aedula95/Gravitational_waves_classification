import pandas_ml as pdl
import pandas as pd
import numpy as np
from sklearn.metrics import accuracy_score
from imblearn.over_sampling import SMOTE
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
from sklearn.neighbors import KernelDensity
from sklearn.neighbors import KNeighborsRegressor
from sklearn.preprocessing import normalize 

def KernelDensityEstimation(disp_arr,sub_matrix,max_count):

	#print len(disp_arr)
	#print sub_matrix
	sub_matrix = np.delete(sub_matrix,[i for i in disp_arr],axis=0)
	sub_matrix = np.delete(sub_matrix,0,axis=1)
	#x = np.delete(sub_matrix,0,axis=1)
	#y = np.delete(sub_matrix,[1,2,3,4],axis=1)

	print max_count
	newlen = max_count - len(sub_matrix)
	print newlen
	kde = KernelDensity(kernel='gaussian')
	#print sub_matrix
	kde.fit(sub_matrix)

	sam = kde.sample(newlen)

	#print sam
		
	return sam


def SampleGenerator(clean_matrix):

	#z = np.delete(clean_matrix,0,axis=1)
	#print clean_matrix
	
	z = clean_matrix

	

	max_count = 0
	for i in xrange(0,len(clean_matrix)):
		if(clean_matrix[i][0] ==1):
			max_count = max_count + 1
	#clean_matrix = np.delete(clean_matrix,0,axis=1)
	arr = []
	#print len(clean_matrix)
	for j in xrange(2,7):
		for i in xrange(0,len(clean_matrix)):
			if(clean_matrix[i][0] != j):
				arr.append(i)
		#print j
		#print len(arr)
		kd = KernelDensityEstimation(arr,clean_matrix,max_count)
		kd = AddLabel(kd,j)
		z = np.concatenate((z,kd),axis=0)
		
		arr = []
	return z	


def AddLabel(kd_sample,label):

	row,col = np.shape(kd_sample)

	new_kd = np.zeros(shape=(row,col+1))

	for i in xrange(0,row):
		new_kd[i][0] = label
	for i in xrange(0,row):
		for j in xrange(0,col):
			new_kd[i][j+1] = kd_sample[i][j]

	return new_kd		


def smote(clean_data):

	new_arr = []
	new_arr2 = []
	new_arr3 = []
	for i in xrange(0,len(clean_data)):
		if(clean_data[i][0] != 1):
			new_arr.append(i)

		if(clean_data[i][0] != 6):
			new_arr2.append(i)

		if(clean_data[i][0] == 6):
			new_arr3.append(i)
		
	ulimit = len(new_arr)

	z = np.delete(clean_data,[i for i in new_arr],axis=0)

	temp = np.delete(clean_data,[i for i in new_arr2],axis=0)

	z = np.append(z,temp,axis=0)

	
	smini = SMOTE(kind='regular',k_neighbors=4,random_state=42)

	x_ini = np.delete(z,[0],axis=1)
	y_ini = np.delete(z,[1,2,3,4],axis=1)

	xres_ini , yres_ini = smini.fit_sample(x_ini,y_ini)
	yres_ini = np.reshape(yres_ini,(len(yres_ini),1))

	
	z = np.append(yres_ini,xres_ini,axis=1)

	z = np.delete(z,[i for i in xrange(0,992)],axis=0)
	
	smfinal = SMOTE(kind='regular',k_neighbors=47,random_state=42)

	z_final = np.delete(clean_data,[i for i in new_arr3],axis=0)
	
	x_final = np.delete(z_final,[0],axis=1)
	y_final = np.delete(z_final,[1,2,3,4],axis=1)

	xres_final , yres_final = smfinal.fit_sample(x_final,y_final)

	yres_final = np.reshape(yres_final,(len(yres_final),1))
	z_final = np.append(yres_final,xres_final,axis=1)

	z = np.append(z_final,z,axis=0)

	return z




'''
	const_matrix = np.delete(clean_data,[i for i in new_arr],axis=0)

	for j in xrange(2,7):
		for i in xrange(0,len(clean_matrix)):
			if(clean_matrix[i][0] != j):
				arr.append(i)
		
		ad = KernelDensityEstimation(arr,clean_matrix)
		ad = AddLabel(kd,j)
		z = np.concatenate((z,kd),axis=0)
		
		arr = []
	return z	

'''




def ClassSeparator(disp_arr,sub_matrix):
	
	sub_matrix = np.delete(sub_matrix,[i for i in disp_arr],axis=0)
	
	


def Normalizer(initclean_data):

	x = np.delete(initclean_data,0,axis=1)
	y = np.delete(initclean_data,[1,2,3,4],axis=1)
	x = normalize(x,norm='l1',axis=0)
	y = np.append(y,x,axis=1)
	return y


def DataClean(filename):
	
	df1 = pd.read_csv(filename)

	lis = [' Jovia',' Terra',' Superterra',' Neptunia',' Subterra',' Mercuria']


	for i in range(0,len(df1)):
		val = df1[' pMassClass'][i]
 	 	if(val == lis[0]):
			df1[' pMassClass'][i] = 1

 	 	if(val == lis[1]):
 	 		df1[' pMassClass'][i] = 2
 	 	if(val == lis[2]):
 	 		df1[' pMassClass'][i] = 3
			
 	 	if(val == lis[3]):
 	 		df1[' pMassClass'][i] = 4

 	 	if(val == lis[4]): 	
 	 		df1[' pMassClass'][i] = 5

 	 	if(val == lis[5]):
 	 		df1[' pMassClass'][i] = 6 	 	

	df2 = df1[[' pMassClass',' sMassSU',' pMassSU','GWCoalesceAmp',' GWFreq']].copy()

	x1= df2.values

	return x1


def classify(xo,yo,xn):
	
	clf = RandomForestClassifier(n_estimators=10, max_depth=None,min_samples_split=2,random_state=0)

	clf.fit(xo,yo.astype(int).ravel())

	print "The mean accuracy score for the training data set is "
	print clf.score(xo,yo.astype(int))

	y_pred =  clf.predict(xn)
	print y_pred
	return y_pred




def DivideData(clean_matrix):
	'''
	clean_matrix = np.delete(clean_matrix,[i for i in arr],axis=0)

	x = np.delete(clean_matrix,0,axis=1)

	y = np.delete(clean_matrix,[1,2,3,4],axis=1)
    
    return (x,y)
    
    '''
	np.set_printoptions(threshold=np.nan)
	#print clean_matrix
	np.random.shuffle(clean_matrix)
	

	
    
	x = np.delete(clean_matrix,0,axis=1)
	y = np.delete(clean_matrix,[1,2,3,4],axis=1)
	xo = np.delete(x,[i for i in range(5039,len(clean_matrix))],axis=0)
	yo = np.delete(y,[i for i in range(5039,len(clean_matrix))],axis=0)

	xn = np.delete(x,[i for i in range(0,5039)],axis=0)
	yn = np.delete(y,[i for i in range(0,5039)],axis=0)

	return (xo,yo,xn,yn)



def ConfusionMatrixResults(y_pred,yn):
	cm = pdl.ConfusionMatrix(yn.astype(int).ravel(),y_pred.astype(int))
	cm.print_stats()

if __name__ == "__main__":

	filename = 'ExoDataPredicted.csv'
	
	cleandata = DataClean(filename)
	#normalize_data = Normalizer(cleandata)

	newdata = smote(cleandata)
	

	#newdata = SampleGenerator(normalize_data)

	

	x_old,y_old,x_new,y_new = DivideData(newdata)

	y_predicted = classify(x_old,y_old,x_new)

	ConfusionMatrixResults(y_predicted,y_new)
	
	#print len(cleandata)
	#print len(newdata)



#print newdata

#print len(newdata)
