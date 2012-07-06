#!/usr/bin/env python
import numpy as np
import pylab as plt
#
# Some plotting of the synthetic and retrieved data

# generate the clean data
doys  = np.arange ( 1, 365+1, 1)[1:]
ndvi_clean =  np.clip(np.sin((doys-1)/72.), 0,1)

# read the state file
state = (open('output/Identity/NDVI_Identity2.params').readlines())[1:]
sdata = np.array([np.array(i.split()).astype(float) for i in state[1:]])[1:]
sNdvi = sdata[:,2][1:]
sdNdvi = sdata[:,4][1:]
doyss = sdata[:,0][1:]
# noisy sample data
noisy = (open('output/Identity/NDVI_fwd2.params_orig').readlines())
ndata = np.array([np.array(i.split()).astype(float) for i in noisy[1:]])
nNdvi = ndata[:,2]
sdnNdvi = ndata[:,3]
ndoys = ndata[:,0]

# retrieved
plt.fill_between(doyss,y1=sNdvi-1.96*sdNdvi,y2=sNdvi+1.96*sdNdvi,facecolor='0.8')
p1 = plt.plot(doyss,sNdvi,'r')
# original
p2 = plt.plot(doys,ndvi_clean,'g')
# samples used
p3 = plt.errorbar(ndoys,nNdvi,yerr=sdnNdvi*1.96,fmt='bo')

plt.legend([p2,p3,p1],['Original state','Sampled noisy state','Retrieved state'])

#plt.show()
plt.savefig('images/example1plot2.png')
