import numpy as np 
from sklearn.decomposition import PCA
filename = './Data/xA.csv'

my_data = np.genfromtxt(filename,delimiter=',',comments='#')
row, col = np.shape(my_data)

test =  np.zeros(shape=(row,col-1))



for i in range(0,row):
    for j in range(0,col-1):

         test[i][j]=my_data[i][j+1]


pca = PCA(n_components=2,svd_solver='arpack')
pca.fit(test)



x = pca.transform(test)

print x 

with open('xA_Analysis.csv','w') as f:
     np.savetxt(f, x,delimiter=',')




