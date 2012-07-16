#!/usr/bin/env python
import matplotlib
matplotlib.use("agg")
import numpy as np
import eoldas
import sys,tempfile
import os

gamma = 50
model_order = 1
dayStep  = 2
BRF = 'data/brdf_WW_1_A_1.kernelFiltered.dat'

names = 'gamma_time,xlai, xhc,  rpl,  xkab, scen, xkw, xkm,   xleafn, xs1,xs2,xs3,xs4,lad'
defaultp = '0.99,0.01,0.01,0.99,0.001,0.99,0.35,1.0,1,0.001,0.001,0.001,5'
default = '%f,%s'%(gamma,defaultp)
priorsd = ''.join(['0.33'] + [',0.33']*(len(names.split(','))-2))

#default = '%f,0.99,5,0.01,0.5,0.001,0.5,0.35,1.5,1,0.001,0,0,5'%gamma

cmd = 'eoldas ' +  \
                '--conf=config_files/eoldas_config.conf  ' + \
                '--conf=config_files/semid_default.conf ' + \
                '--parameter.x.default=%s '%default + \
                '--parameter.solve=0,1,0,0,1,1,1,1,1,1,1,0,0,0 ' +\
                '--logfile=logs/file.log ' + \
                '--parameter.limits=[[1,365,%d]] '%dayStep +\
                '--parameter.result.filename=output/gamma%f/state.dat '%gamma + \
                '--operator.obs.y.result.filename=output/gamma%f/obs.dat '%gamma + \
                '--operator.obs.y.state=%s '%BRF + \
                '--operator.prior.y.sd=%s '%priorsd + \
                '--operator.prior.y.state=%s '%defaultp + \
                ' --operator.modelt.rt_model.model_order=%d '%model_order + \
                '--calc_posterior_unc'

# '--parameter.x.state=output/gamma%f/state.dat '%gamma + \
# '--conf=confs/semid_default_state.conf' + \

# initialise eoldas
self = eoldas.eoldas(cmd)
# solve DA
self.solve(write=True)

