import numpy as np
import matplotlib.pyplot as plt
from pycbc.waveform import get_td_waveform
import pylab
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model
from scipy.interpolate import interp1d
def computeHPHC(argList):
	for apx in ['SEOBNRv2']:
		hp, hc = get_td_waveform(approximant=apx, \
			mass1 = argList[0], \
			mass2 = argList[1], \
			spin1z = argList[2], \
			delta_t = argList[3],\
			f_lower = argList[4])
		plt.plot(hp.sample_times, hp, label=apx)
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
	x = np.array(ampPeaks['posTime'])
	x = np.append(x,np.array(ampPeaks['negTime']))
	y = np.array(ampPeaks['posAmp'])
	y = np.append(y,np.array(ampPeaks['negAmp']))
	z = interp1d(x,y,kind='cubic')
	x_new = np.linspace(x[0], x[-1], 10000)
	f2 = z(x_new)
	#y_new = z(x_new)
	plt.plot(x,y,'o', x_new, f2)
	#plt.xlim([x[0]-1, x[-1] + 1 ])
	plt.scatter(ampPeaks['posTime'],ampPeaks['posAmp'])	
	return ampPeaks

if __name__ == "__main__":
	delT = 1.0/16384
	hp, hc = computeHPHC([10, 2, 0.0, delT, 30])
	# hp = [t*10**19 for t in hp]
	#print peaks(hp, delT)['posAmp']
	count = 0
	for i in peaks(hp,delT)['posAmp']:
         count=count + 1
         

print count
     
plt.show()
#plt.ylabel('Strain')
#plt.xlabel('Time (s)')
#plt.legend()
#pylab.show()