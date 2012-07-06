#!/usr/bin/env python
import pdb
import numpy as np

if __name__ == "__main__":
    import sys,tempfile
    this = sys.argv[0]
    import  eoldas
    from solve_eoldas_identity import *
    # import the setup methods from solve_eoldas_identity
    import pylab as plt

    # SD of noise
    noise=0.15
    # nominal sampling period
    n_per=7
    # proportion of missing samples
    obs_off=0.33
    # order of differential model (integer)
    model_order=1
    # sgwindow is larger to create larger data gaps
    sgwindow=10

    # set up data for this experiment
    file, ideal_gamma,doys,ndvi_clean,ndvi,qa_flag  = \
		prepare(noise=noise,n_per=n_per,\
			obs_off=obs_off,model_order=model_order,\
			sgwindow=sgwindow)

    # set gamma to thge theoretical value
    gamma = ideal_gamma

    # set op file names
    xfile = 'output/Identity/NDVI_Identity_a.params'
    yfile = 'output/Identity/NDVI_fwd_a.params'

    # initialise options for DA overriding any in config files
    # make sure we use some different output file names to othe scripts
    cmd = 'eoldas ' +  \
                ' --conf=config_files/eoldas_config.conf --conf=config_files/Identity.conf ' + \
                ' --logfile=mylogs/Identity.log ' + \
                ' --calc_posterior_unc ' + \
                ' --parameter.solve=0,1 ' + \
		' --parameter.result.filename=%s '%xfile +\
                ' --parameter.x.default=%f,0.0 '%(gamma) + \
		' --operator.obs.y.result.filename=%s'%yfile +\
                ' --operator.obs.y.state=%s'%file+\
                ' --operator.modelt.rt_model.model_order=%d '%model_order 

    # initialise eoldas
    self = eoldas.eoldas(cmd)
    # solve DA
    self.solve(write=True)
    
    # now pull some data out of the eoldas

    # the 'root' of the DA data structure is in self.solver.root

    root = self.solver.root

    # The state vector data are stored in root.x
    # with ancillary information in root.x_meta
    # so the state vector is e.g. in root.x.state
    # and the names are in root.x_meta.state

    state_names = root.x_meta.state
    state = root.x.state

    # The sd representation of the posterior is in root.x.sd
    # This is all set up in eoldas_Solver.py
    # All storage is of type ParamStorage, an extended
    # dictionary structure. You can explore it interactively
    # with e.g. root.dict().keys() or self.keys() since
    # self here is a straight dictionary
    sd = root.x.sd

    # The full inverse Hessian is in self.solver.Ihessian
    Ihessian = self.solver.Ihessian

    # A mask to reduce this to only the state variables
    # being targeted (solve == 1) is through a call to:
    NDVI_Ihessian = self.solver.unloader(None,Ihessian,M=True)

    # This is of shape (365,365) here.
    # so now lets produce an image of it
    # to visualise the structure
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cax = ax.imshow(NDVI_Ihessian,interpolation='nearest')
    ax.set_title('Posterior Uncertainty matrix for NDVI')
    # Add colorbar, make sure to specify 
    # tick locations to match desired ticklabels
    cbar = fig.colorbar(cax) #), ticks=[-1, 0, 1])
    #cbar.ax.set_yticklabels(['< -1', '0', '> 1'])# vertically oriented colorbar
    # see http://matplotlib.sourceforge.net/plot_directive/mpl_examples/pylab_examples/colorbar_tick_labelling_demo.py
    # save it   
    plt.savefig('output/IHessianNDVI_expt1.png')

 
   
 

    
