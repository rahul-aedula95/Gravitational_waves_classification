import pylab
from pycbc.waveform import get_td_waveform

for apx in ['SEOBNRv2']:
    hp, hc = get_td_waveform(approximant=apx,
                                 mass1=53,
                                 mass2=43	,
                                 spin1z=0.0,
                                 delta_t=1.0/4096,
                                 f_lower=10)

    pylab.plot(hp.sample_times, hp, label=apx)

print 'max(hp):', max(hp)
print 'min(hp):', min(hp)
print hp

pylab.ylabel('Strain')
pylab.xlabel('Time (s)')
pylab.legend()
pylab.show()
