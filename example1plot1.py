#!/usr/bin/env python
import numpy as np
import pylab as plt
#
# Some plotting of the synthetic and retrieved data

# generate the clean data
doys  = np.arange ( 1, 365+1, 1)
ndvi_clean =  np.clip(np.sin((doys-1)/72.), 0,1)

# read the state file
state = open('output/Identity/NDVI_Identity1.params').readlines()
sdata = np.array([np.array(i.split()).astype(float) for i in state[1:]])
sNdvi = sdata[:,2]
sdNdvi = sdata[:,4]

# noisy sample data
noisy = open('output/Identity/NDVI_fwd1.params_orig').readlines()
ndata = np.array([np.array(i.split()).astype(float) for i in noisy[1:]])
nNdvi = ndata[:,2]
sdnNdvi = ndata[:,3]
ndoys = ndata[:,0]

# retrieved
plt.fill_between(doys,y1=sNdvi-1.96*sdNdvi,y2=sNdvi+1.96*sdNdvi,facecolor='0.8')
p1 = plt.plot(doys,sNdvi,'r')
# original
p2 = plt.plot(doys,ndvi_clean,'g')
# samples used
p3 = plt.errorbar(ndoys,nNdvi,yerr=sdnNdvi*1.96,fmt='bo')

plt.legend([p2,p3,p1],['Original state','Sampled noisy state','Retrieved state'])

#plt.show()
plt.savefig('images/example1plot1.png')
