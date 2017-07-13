import pandas_ml as pdl
from imblearn.over_sampling import     SMOTE 
import pandas as pd
import numpy as np
from sklearn.naive_bayes import GaussianNB
from sklearn.decomposition import PCA
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import AdaBoostClassifier
import subprocess
from imblearn.over_sampling import     RandomOverSampler 
from sklearn.neighbors import KernelDensity
filename = 'ExoDataPredicted.csv'

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






#x1 = np.delete(x1,[i for i in arr],axis=0)

x = np.delete(x1,0,axis=1)

y = np.delete(x1,[1,2,3,4],axis=1)


'''
print len(x)
print len(y)

print x
print y


kde = KernelDensity(kernel='gaussian')
kde.fit(x,y)

sam = kde.sample(max_count)

print sam 

print max_count
'''

#sm = SMOTE(kind='regular')
#ros = RandomOverSampler(random_state=0)
#x,y = ros.fit_sample(x, y.ravel())



#print len(x)

#print len(y)



xo = np.delete(x,[i for i in range(3000,len(df1))],axis=0)
yo = np.delete(y,[i for i in range(3000,len(df1))],axis=0)

xn = np.delete(x,[i for i in range(0,2999)],axis=0)
yn = np.delete(y,[i for i in range(0,2999)],axis=0)


clf = RandomForestClassifier(n_estimators=10, max_depth=None,min_samples_split=2,random_state=0)

clf.fit(xo,yo.astype(int).ravel())

print "The mean accuracy score for the training data set is "
print clf.score(xo,yo.astype(int))


y_pred =  clf.predict(xn)
''

cm = pdl.ConfusionMatrix(yn.astype(int).ravel(),y_pred)

cm.print_stats()

#print df2
