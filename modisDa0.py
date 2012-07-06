import numpy as np
import eoldas
import sys,tempfile
import os

gamma = 100/np.sqrt(28.)
model_order = 1
BRF = 'data/brdf_WW_1_A_1.kernelFiltered.dat'

names = 'gamma_time,xlai, xhc,  rpl,  xkab, scen, xkw, xkm,   xleafn, xs1,xs2,xs3,xs4,lad'
default = '%f,0.99,5,0.01,0.99,0.001,0.99,0.35,1.5,1,0.001,0,0,5'%gamma

#default = '%f,0.99,5,0.01,0.5,0.001,0.5,0.35,1.5,1,0.001,0,0,5'%gamma

cmd = 'eoldas ' +  \
                '--conf=config_files/eoldas_config.conf  ' + \
                '--conf=config_files/semid_default.conf ' + \
                '--parameter.x.default=%s '%default + \
                '--parameter.solve=0,1,0,0,1,1,1,1,1,1,1,0,0,0 ' +\
                '--logfile=logs/file.log ' + \
                '--parameter.limits=[[1,365,28]] ' +\
                '--parameter.result.filename=output/gamma%f/state.dat '%gamma + \
                '--operator.obs.y.result.filename=output/gamma%f/obs.dat '%gamma + \
                '--operator.obs.y.state=%s '%BRF + \
                '--operator.modelt.rt_model.model_order=%d '%model_order + \
                '--calc_posterior_unc'

# initialise eoldas
m = eoldas.eoldas(cmd)
# solve DA
m.solve(write=True)

