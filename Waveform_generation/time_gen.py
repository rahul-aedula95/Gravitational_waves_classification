import numpy as np
import pylab
from pycbc.waveform import get_td_waveform


filename = "sampletime_"
m1=10
m2=10
m1_step=5    #step in change of mass1
m2_step=5    #step in change of mass2
#extrapolate the step depending on what is required
count=1     # file name count for creation 

#Add more loops to add more features to generate more samples 
for apx in ['SEOBNRv2']:
    for i in xrange(0,2):
        for j in xrange(0,2):
           hp, hc = get_td_waveform(approximant=apx,
                                 mass1=m1,
                                 mass2=m2,
                                 spin1z=0.6,
                                 delta_t=1.0/16384,
                                 f_lower=45,f_ref=200)

           
           with open(filename+str(count),'w') as f:
                 np.savetxt(f, hp.sample_times)
            m1=m1+m1_step
            count=count+1
        m2=m2+m2_step    




print done



