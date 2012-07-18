#!/usr/bin/env python
"""
A script to process MODIS surface BRDF data using EOLDAS. note that
this method uses a different strategy for xvalidation to the 
Lewis et al 2012 paper: We first invert each independent observation,
and smoothly interpolate to get a continuous time series of parameters
This is then used to feed into EOLDAS as a starting first guess. We
hope that this guess has already put the minimisation of the cost
function close to the minimum, so it's relatively fast to optimise
the full cost function, and then calculate the Hessian, etc.
"""
import numpy as np
import matplotlib.pyplot as plt
from sentinel import Sentinel


def main (ifile = 'data/brdf_WW_1_A_1.kernelFiltered.dat', \
        confFile='config_files/semid_default.conf', \
        solve = ['xlai','xkab','scen','xkw','xkm','xleafn','xs1']) :
  '''
  Show that we can use this same setup to solve 
  for the field data (MODIS) by using a different config file
  '''
  
  
  ifile = 'data/brdf_WW_1_A_1.kernelFiltered.dat'
  ofileSingle = 'output/brdf_WW_1_A_1.kernelFilteredSingle.dat'
  ofile = 'output/brdf_WW_1_A_1.kernelFiltered.dat'

  s = Sentinel(solve=solve,confFile=confFile)
  # as above, solve for initial estimate
  s.solveSingle(ifile,ofileSingle)
  s.paramPlot(s.loadData(ofileSingle),s.loadData(ofileSingle),\
                 filename='plots/testFieldDataSingle.png')
  s.smooth(ofileSingle,ofile=ofileSingle+'_smooth')
  s.paramPlot(s.loadData('input/truth.dat'),\
              s.loadData(ofileSingle+'_smooth'),\
              filename='plots/%s.png'%(ofileSingle+'_smooth'))

  gamma = np.median(np.sqrt(s.gammaSolve.values()) * \
        np.array(s.wScale.values()))
  gamma = int((np.sqrt(s.gammaSolve[solve[0]]) * \
        np.array(s.wScale[solve[0]]))+0.5)
  s.solveRegular(ifile,ofile,modelOrder=2,gamma=gamma, \
        initial=ofileSingle+'_smooth')
  s.paramPlot(s.loadData(ofileSingle+'_smooth'),\
              s.loadData(ofile + '_result.dat'),\
              filename='plots/%s_Gamma%08d.png'%(ofile,gamma))

if __name__ == "__main__":
    main()