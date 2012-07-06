***************************************************************************
A synthetic experiment: simulating the performance of Sentinel-2
***************************************************************************

A useful application of the EOLDAS tools is to simulate how different 
sensors or combinations of sensors might be expected to work together.
Simply, we could generate a set of artificial observations, acquired
by a sensor (or sensors) from which the spectral, noise and orbital 
acquisition characteristics are known (or we can try different set ups
to see their trade offs). We can then use EOLDAS to invert these
synthetic observations and compare to the reality that went into them.
This type of experiments can be very useful to test, for example, the
usefulness of constellation concepts. The current example is a 
reduced version of the experiments that are presented in 
`Lewis et al (2010) <http://dx.doi.org/10.1016/j.rse.2011.12.027>`_, and
show the use of EOLDAS in simulating ESA's upcoming Sentinel-2 performance
for the inversion of biophysical parameters. In this document, we try
to reproduce the middle column in Figure 2 in `Lewis et al (2010) <http://dx.doi.org/10.1016/j.rse.2011.12.027>`_

The command for this is ::
    
    eoldas_run.py --conf=config_files/sentinel.conf --no_calc_posterior_unc
    

        





