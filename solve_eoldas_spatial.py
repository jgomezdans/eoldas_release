#!/usr/bin/env python
import pdb
import numpy as np

def gauss_kern(size, sizey=None):
    """ 
    Returns a normalized 2D gauss kernel array for convolutions 
    From http://www.scipy.org/Cookbook/SignalSmooth
    """
    import numpy as np
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def create_data ( nr,nc,n_per=1, noise=0.15, obs_off=0.33, \
					window_size=3):
    """
    Create synthetic "NDVI-like" data for a fictitious space series. We return
    the original data, noisy data (using IID Gaussian noise), the QA flag as well
    as the time axis.
    
    Missing observations are simulated by drawing a random number between 0 and 1
    and checking against obs_off.
    
    Parameters
    ----------
    nc : number of columns
    nr : number of rows

    n_per : integer
    Observation periodicity. By default, assumes every sample

    noise : float
    The noise standard deviation. By default, 0.15

    obs_off : float
    The threshold to decide on missing observations in the series.

    
    """
    import numpy as np
    from scipy import signal

    r,c = np.mgrid[1:nr+1:n_per,1:nc+1:n_per]
    ndvi_clean =  np.clip(np.sin(np.pi*2*(r-1)/nr)*np.sin(np.pi*2*(c-1)/nc), 0,1) 
    ndvi = ndvi_clean.copy()
    # add Gaussian noise of sd noise
    ndvi = np.random.normal(ndvi,noise,ndvi.shape)
     
    # set the qa flags for each sample to 1 (good data)
    qa_flag = np.ones_like ( ndvi).astype( np.int32 )
    passer = np.random.rand(ndvi.shape[0],ndvi.shape[1])
    if window_size >0:
        # force odd
        window_size = 2*(window_size/2)+1
        # smooth passer
        g = gauss_kern(window_size)
        passer = signal.convolve(passer,g,mode='same')
    # assign a proportion of the qa to 0 from an ordering of the smoothed 
    # random numbers
    qf = qa_flag.flatten()
    qf[np.argsort(passer,axis=None)[:passer.size * obs_off]]  = 0
    qa_flag = qf.reshape(passer.shape)
    return ( r,c , ndvi_clean, ndvi, qa_flag )

def prepare(nr=50,nc=100,noise=0.15,n_per=1,obs_off=0.33,model_order=1,\
					filtwindow=3):
    import os,tempfile

    # Start by creating some univarite time series data
    ( r, c, ndvi_clean, ndvi, qa_flag ) = \
                create_data ( nr,nc,n_per=n_per, noise=noise, \
			obs_off=obs_off,\
			window_size=filtwindow)

    # calculate the gradient of the ndvi
    dydx = np.gradient(ndvi_clean.flatten())/np.gradient(c.flatten())
    # to the model_order-th order
    for i in np.arange(1,model_order):
        dydx  = np.gradient(dydx)
    # calculate the 'true' gamma
    model_unc = np.sqrt(np.mean(dydx**2))
    gamma = 1./model_unc

    temp, fp_brf = tempfile.mkstemp()    
    os.write (temp, \
        "#PARAMETERS mask row col NDVI sd-NDVI\n")
    [ os.write (temp,"%d %f %f %f %f\n" % \
        ( int(qa_flag.flatten()[i]), r.flatten()[i], c.flatten()[i], \
          ndvi.flatten()[i], noise)) \
        for i in xrange (len(ndvi.flatten()))]
    os.close(temp)

    return fp_brf,gamma,r,c,ndvi_clean,ndvi,qa_flag

if __name__ == "__main__":
    import sys,tempfile
    this = sys.argv[0]
    import eoldas

    # SD of noise
    noise=0.15
    # nominal sampling period
    n_per=1
    # proportion of missing samples
    obs_off=0.333
    # order of differential model (integer)
    model_order=1
    # sgwindow is larger to create larger data gaps
    filtwindow=13

    # set up data for this experiment
    file, ideal_gamma,r,c,ndvi_clean,ndvi,qa_flag  = \
		prepare(noise=noise,n_per=n_per,\
			obs_off=obs_off,model_order=model_order,\
			filtwindow=filtwindow)

    # set gamma to less than the the theoretical value
    gamma = ideal_gamma*0.2
    # set op file names
    xfile = 'output/Identity/NDVI_Identity_spatial.params'
    yfile = 'output/Identity/NDVI_fwd_spatial.params'

    # initialise optuions for DA overriding any in config files
    cmd = 'eoldas_run.py ' +  \
                ' --conf=config_files/eoldas_config.conf --conf=config_files/IdentitySpatial.conf ' + \
                ' --logfile=mylogs/IdentitySpatial.log ' + \
                ' --calc_posterior_unc ' + \
                ' --parameter.solve=0,0,1 ' + \
		' --parameter.result.filename=%s '%xfile +\
                ' --parameter.x.default=%f,%f,0.0 '%(gamma,gamma) + \
		' --operator.obs.y.result.filename=%s'%yfile +\
                ' --operator.obs.y.state=%s'%file+\
                ' --parameter.x.state=%s'%file+\
                ' --operator.modelr.rt_model.model_order=%d '%model_order +\
                ' --operator.modelc.rt_model.model_order=%d '%model_order

    # initialise eoldas
    self = eoldas.eoldas(cmd)
    # solve DA
    self.solve(write=True)
    # ok, now read the results back in
    fwd = (np.array([i.split() for i in open(xfile).readlines()[1:]]).astype(float))
    fwdndvi = fwd[:,5].reshape(ndvi_clean.shape)
    fwdndvisd = fwd[:,8].reshape(ndvi_clean.shape)

    import pylab as plt

    plt.clf()
    plt.xlabel('original ndvi')
    plt.ylabel('modelled ndvi')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.plot([0,1],[0.,1.],'k')
    #plt.plot(ndvi[qa_flag].flatten(),fwdndvi[qa_flag].flatten(),'r.',label='fwd at obs')
    plt.errorbar(ndvi_clean.flatten(),fwdndvi.flatten(),yerr=fwdndvisd.flatten()*1.96,fmt='g.') 
    plt.plot(ndvi_clean.flatten(),fwdndvi.flatten(),'r.')
    #plt.plot(ndvi_clean.flatten(),fwdndvi.flatten(),'g.',label='fwd')
    plt.savefig('%s.png'%xfile)

    plt.clf()
    plt.xlabel('original ndvi')
    plt.ylabel('modelled ndvi')
    plt.plot([0,1],[0.,1.],'k')
    plt.xlim(0,1)
    plt.ylim(0,1)
    #plt.plot(ndvi[qa_flag].flatten(),fwdndvi[qa_flag].flatten(),'r.',label='fwd at obs')
    plt.errorbar(ndvi_clean.flatten(),fwdndvi.flatten(),yerr=fwdndvisd.flatten()*1.96,fmt='g.') 
    plt.plot(ndvi_clean.flatten(),fwdndvi.flatten(),'r.')
    plt.errorbar(ndvi[qa_flag].flatten(),fwdndvi[qa_flag].flatten(),\
                 yerr=(fwdndvisd[qa_flag].flatten()*1.96),xerr=(fwdndvisd[qa_flag].flatten()*0+noise*1.96),fmt='b.')
    plt.plot(ndvi[qa_flag].flatten(),fwdndvi[qa_flag].flatten(),'g.')
    plt.xlim(0,1.0)
    plt.ylim(0,1.0)
    plt.savefig('%s_withdata.png'%xfile)
    plt.clf()
    plt.imshow(ndvi_clean,interpolation='nearest',vmin=-1,vmax=1)
    plt.colorbar()
    plt.savefig('%s_clean.png'%xfile)
    plt.clf()
    plt.imshow(fwdndvisd,interpolation='nearest',vmin=0,vmax=0.15)
    plt.colorbar()
    plt.savefig('%s_fwdsd.png'%xfile)

    plt.clf()
    plt.imshow(fwdndvi,interpolation='nearest',vmin=-1,vmax=1)
    plt.colorbar()
    plt.savefig('%s_fwd.png'%xfile)
    plt.clf()
    xx = ndvi.copy()
    xx[(1-qa_flag).astype(bool)] = -1
    plt.imshow(xx,interpolation='nearest',vmin=-1,vmax=1)
    plt.colorbar()
    plt.savefig('%s_data.png'%xfile)
    plt.clf()
    x = np.arange(len(fwdndvi[13])) + 1
    plt.plot(x,ndvi_clean[13],label='clean')
    plt.plot(x,xx[13],'.',label='data')
    plt.plot(fwdndvi[13],label='fwd')
    plt.ylim(-0.1,1)
    plt.errorbar(x,fwdndvi[13],yerr=(fwdndvisd[13]))
    plt.xlabel('column')
    plt.ylabel('ndvi')
    plt.legend()
    plt.savefig('%s_transect.png'%xfile)

