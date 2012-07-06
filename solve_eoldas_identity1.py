#!/usr/bin/env python
import pdb
import numpy as np

def create_data ( n_per=4, noise=0.15, obs_off=0.33, \
					window_size=0, order=4):
    """
    Create synthetic "NDVI-like" data for a fictitious time series. We return
    the original data, noisy data (using IID Gaussian noise), the QA flag as well
    as the time axis.
    
    Missing observations are simulated by drawing a random number between 0 and 1
    and checking against obs_off.
    
    Parameters
    ----------
    n_per : integer
    Observation periodicity. By default, assumes every 8 days

    noise : float
    The noise standard deviation. By default, 0.15

    obs_off : float
    The threshold to decide on missing observations in the time series.

    window_size : integer, odd
    window size for savitzky_golay filtering. A large window size will lead
    to larger data gaps by correlating the noise. Set to zero by default
    which applies no smoothing.

    order : integer
    order of the savitzky_golay filter. By default 4.
    
    """
    from savitzky_golay import savitzky_golay
    import numpy as np

    doys  = np.arange ( 1, 365+1, n_per)
    ndvi_clean =  np.clip(np.sin((doys-1)/72.), 0,1) 
    ndvi =  np.clip(np.sin(doys/72.), 0,1) 
    # add Gaussian noise of sd noise
    ndvi = np.random.normal(ndvi,noise,ndvi.shape[0])
     
    # set the qa flags for each sample to 1 (good data)
    qa_flag = np.ones_like ( ndvi).astype( np.int32 )
    passer = np.random.rand( ndvi.shape[0])
    if window_size >0:
        # force odd
        window_size = 2*(window_size/2)+1
        passer = savitzky_golay(passer, window_size=window_size, order=order)
    # assign a proportion of the qa to 0 from an ordering of the smoothed 
    # random numbers
    qa_flag[np.argsort(passer)[:passer.size * obs_off]]  = 0

    return ( doys, ndvi_clean, ndvi, qa_flag )

def prepare(noise=0.15,n_per=7,obs_off=0.33,model_order=1,\
					sgwindow=0,sgorder=4):
    import os,tempfile

    # Start by creating some univarite time series data
    ( doys, ndvi_clean, ndvi, qa_flag ) = \
                create_data ( n_per=n_per, noise=noise, \
			obs_off=obs_off,\
			window_size=sgwindow, order=sgorder)

    # calculate the gradient of the ndvi
    dydx = np.gradient(ndvi_clean)/np.gradient(doys)
    # to the model_order-th order
    for i in np.arange(1,model_order):
        dydx  = np.gradient(dydx)
    # calculate the 'true' gamma
    model_unc = np.sqrt(np.mean(dydx**2))
    gamma = 1./model_unc

    temp, fp_brf = tempfile.mkstemp()    
    os.write (temp, "BRDF %d 1 NDVI %s\n" % (len(ndvi),str(noise)) )
    [ os.write (temp,"%f %d 0 0 0 0 %f\n" % ( doys[i], qa_flag[i], ndvi[i])) \
        for i in xrange (len(ndvi))]
    os.close(temp)

    return fp_brf,gamma,doys,ndvi_clean,ndvi,qa_flag

if __name__ == "__main__":
    import sys,tempfile
    this = sys.argv[0]
    import eoldas

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

    # set gamma to less than the the theoretical value
    gamma = ideal_gamma*0.45

    # set op file names
    xfile = 'output/Identity/NDVI_Identity1.params'
    yfile = 'output/Identity/NDVI_fwd1.params'

    # initialise optuions for DA overriding any in config files
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
    
