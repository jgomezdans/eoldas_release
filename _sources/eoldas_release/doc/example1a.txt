Application to MODIS reflectance data
------------------------------------------
We can now apply the concepts demonstrated above to real EO data from the MODIS sensors.

This is actually quite trivial to achieve as it involves the same basic
configuration file as previously. This time, we will access it from ``eoldas_run.py``.


An example data file then is ``data/modis_botswana.dat`` which is MODIS reflectance data for a site in Botswana. The format of this file is
'BRDF', which looks like:

.. literalinclude:: ../data/modis_botswana.dat
  :lines: 1-23
  :language: python

The first line is a header, as with the previous 'PARAMETERS' format we saw.
The first word on the header must be ``BRDF`` (a # can be included for compatitbilty with
some other formats).

The second item on the header line is the number of samples in the file (92 here). The third item is the 
number of wavebands represented (7 here) then this is followed by wavelengths, wavebands, or waveband identifiers.

If a single float number is given, it is assumed to represent a narrow waveband centred around
that wavelength. Wavelength is assumed to be in nm. If it is two float numbers connected by a dash (e.g. 45-550)
then it is taken to be a tophat function representation of a waveband, where the two number represent 
the minimum and maximum wavelength respectively.
If some other code is found (e.g. ``B2``), it is assumed to be a tag associated with a particular waveband.
In this case, any spectral libraries loaded in the configuration files are searched to see if they contain
a definition of this waveband. If nothing is found and it cannot be interpreted as a number, it is assigned an arbitrary wavelength (index starting from 0) and should not be interpreted as a wavelength identifier.

The 7 columns following the wavebands specify the standard deviation in the reflectance data.
All subsequent lines contain reflectance, as a function of 'locational', 'control' and spectral information. The first column
here defines location (time), columns 2 to 6 are 'control' variables: things the observation waries as a function of that are not
spectral or locational. In this case, they are (all angles are in degrees): 

# a data mask (1 for good),
# view zenith angle, 
# view azimuth angle, 
# solar zenith angle,  
# solar azimuth angle. 

The final 7 columns give the reflectance in the 7 wavebands for each location defined.

The data span days of year 181 to 272 inclusive. There is a lot of day to day variation in the data, but a clear
underlying trend in the Near infrared. We will set ourselves the task of trying to use eoldas to filter the dataset
so that the trend becomes more apparent. We will also try to extrapolate from the sample days to the whole year (which is clearly
quite a challenge, but instructive for understanding uncertainty).

We will apply the same data assimilation components as previously, i.e. an Identity operator for the NIR reflectance
(in other words, we take the reflectance as a prior estimate of what we wish to estimate) and a regularisation filter
implemented as a constraint on first order derivatives at lag one day. We will suppose that we know the uncertainty in this
model (the temporal constraint) to be expressed as a standard deviation of .002 (1/500), i.e. the root mean squared deviation from
one day to the next is expected to be around .002. An examination of the dataset will show that the NIR reflectance roughly
increases from around 0.1 to 0.2 over about 100 days, so we could use a gamma value of 1000. However, we don't want to impose that constraint too strongly, so we will instead use a gamma of 500. We are not trying to fit a trendline
here, just to smooth the dataset so that we can detect the underlying behaviour. The result should however not be very sensitive to this guess and can always be explored using approaches such as cross validation.
  
We have a configuration file set up in ``config_files/Identity2.conf``. This provides a full
description of the problem we wish to solve, which is to reduce noise in the NIR observation series.

.. literalinclude:: ../config_files/Identity2.conf

In this, we can see the required location field defined in the ``[parameters]`` section (i.e. ``parameters.location=[[1,365,1]]``).
The control variables are associated with  the observation operator, here given as ``[mask,vza,vaa,sza,saa]``.

To run eoldas with this configuration then all we need do is type::

    eoldas_run.py --conf=config_files/eoldas_config.conf \
      --conf=config_files/Identity2.conf --calc_posterior_unc


This writes the files:

``output/Identity/MODIS_botswana.params``: state estimation (smoothed reflectance)

.. literalinclude:: ../output/Identity/MODIS_botswana.params
    :lines: 182-192

``output/Identity/Botswana_fwd.params``: forward modelling of y

.. literalinclude:: ../output/Identity/Botswana_fwd.params
    :lines: 1-10

as well as the orignal data in ``output/Identity/Botswana_fwd.params_orig``

.. literalinclude:: ../output/Identity/Botswana_fwd.params_orig
    :lines: 1-10


with appropriate graphics for the state:

.. image:: ../output/Identity/MODIS_botswana.params.plot.x.png

and the y data:

.. image:: ../output/Identity/Botswana_fwd.params.plot.y2.png

.. image:: ../output/Identity/Botswana_fwd.params.plot.x.png

The resultant state data are quite instructive: where we have observations, the uncertainty
is reduced from 0.015 to around 0.004 (the actual degree of noise reduction
depends on the value of gamma used). Where there are no data, the uncertainty 
grows to around 0.017. It is slightly reduced at the year start/end because of the wrparound condition used here.

In this case, we have used a first order differential constraint with a periodic boundary condition. These
are quite important in this case: we only have observations in a limited time window, so using a periodic
boundary condition is one way to place some form of constraint at what happens when we have no data. The
first order differential model will in essence perform a linear interpolation where there are no data, which is
probably appropriate for this case. You can try changing the model order to see what happens. Run e.g.::

  eoldas_run.py --conf=config_files/eoldas_config.conf \
    --conf=config_files/Identity2.conf --calc_posterior_unc \
    --operator.modelt.rt_model.model_order=2 \
    --parameter.x.default=5000,0.1 \
    --operator.obs.y.result.filename=output/Identity/Botswana_fwd.params2 \
    --parameter.result.filename=output/Identity/MODIS_botswana.params2
  
and have a look at ``output/Identity/MODIS_botswana.params2.plot.x.png``:

.. image:: ../output/Identity/MODIS_botswana.params2.plot.x.png

The assumed behaviour is quite different to the first order differential constraint outside of the observations. With a high value of gamma, the result is essentially a straight line where there are observations.
The influence of the wraparound condition is also clear here.


If we used a lower gamma, we would see some features of using a second order model::

    eoldas_run.py --conf=config_files/eoldas_config.conf \
      --conf=config_files/Identity2.conf --calc_posterior_unc \
      --operator.modelt.rt_model.model_order=2 \
      --parameter.x.default=200,0.1  \
      --operator.obs.y.result.filename=output/Identity/Botswana_fwd.params3  \
      --parameter.result.filename=output/Identity/MODIS_botswana.params3
    


.. image:: ../output/Identity/MODIS_botswana.params3.plot.x.png


Where there are observations, the second order difference constraint would have high frequency oscillations. This can be a positive feature or an annoyance: it all depends on how you expect the function to behave. It is
not an arbitrary choice: the user has imposed a particular expectation of behaviour here through the model. There clearly needs to be some evidence for choosing one form of model over another, although that is not always straightforward.

One other interesting feature of this result is that the mean estimate is bounded (0.0,1.0) which is reasonable for reflectance data (theoretically, BRF can go above 1, but this is rarely observed). This means that the mean value (the red line) is forced to stay above 0 (days 100 to 150) even though the apparent trend from the observations might otherwise make it go below zero. This condition is imposed in the configuration by the line:

.. literalinclude:: ../config_files/Identity2.conf
    :lines: 25-25

One further thing to note about this result is that when the uncertainty is so large (as is the case here when we are extrapolating) that the confidence interval spans more than the entire data range (0,1) the mean and standard deviation described in the Gaussian statistics are a little meaningless or at least should be interpreted with caution. It might in this case be better to describe the data (in the extrapolated region) as simply being somewhere between the physical limits.

A final comment is that actually, a large amount of the departure of the signal from our smooth trajectory, the high frequency variation that we have treated as noise here, can most likely be described by considering the physics os scattering by the surface (there are, on the whole, BRDF effects). We will return to this later.

As well as raising a few issues with regard to model selection, we have demonstrated further use of the eoldas command line for quickly changing some terms in an experiment.
However, it is all very well in showing that eoldas can ingest satellite data and
apply constraints to provide an estimate of state variables, but the
state variables here (reflectance) do not directly help us monitor the
properties of the land surface or link to process models.

To do that using EO data, we generally need more complex obeseravtion operators. These are dealt with in the next section.

