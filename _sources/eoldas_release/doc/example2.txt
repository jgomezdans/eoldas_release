***********************************************************
Radiative transfer modelling for optical remote sensing
***********************************************************

Radiative transfer modelling for optical remote sensing
=======================================================

Introduction
------------

Our main focus will be on using the EOLDAS data assimilation techniques to
understand the land surface using optical remote sensing data. Optical sensors
capture reflected solar radiation in the 400 to 2500nm range, over a number
of wavebands. The signal on these bands is a combination of processes that
influence the scattering of photons, such as interactions of photons within
vegetation canopies, interactions with aerosols in the atmosphere, gaseous
absorption processes, etc. Since all these processes contribute to the signal,
it is important to have a physical understanding to separate their contributions
and produce accurate estimates of the combined state of the atmosphere and
land surface.

Typically, RT models of the land surface model the interactions of photons and
plant components (leaves, trunks, etc.) above a surface layer (usually soil, but
also snow). Atmospheric radiative transfer models are needed to account for
water vapour and ozone absorption in the atmosphere, as well as for the effect
of aerosols. One can think of this arrangement as a stacking of models, going
from bottom to top as: soil model, vegetation model and atmospheric model.

Our main focus is the study of the land surface using DA techniques. 
The proposed RT models predict, for a given wavelength and acquisition geometry, directional
reflectance of the land surface. The parameters that govern these models are
typically to do with soil structure (soil brightness and roughness terms) and
the vegetation architecture (leaf area index (LAI), leaf angle distribution,
...),
as well as leaf biochemical parameters (chlorphyll a+b concentration, dry
matter,
leaf structural parameters...). These paraeters are the ones that form our land
surface state vector, and they will vary with time and location, but provide a
full description of the surface.

We note that these models are highly non-linear, and combined with only modest
observations available in terms of angular and spectral sampling, as well as the
contribution of thermal noise in the acquisitions, result in an*ill-posed
problem*: typically, the available observations do not provide enough constrains
to the models, and this results in state vector elements being poorly
determined, equivalently having large uncertainties. We shall see that these
situation can be improved by the use of data assimilation techniques.

We will first consider a simple RT model for continuous canopies under a soil
layer, with a leaf optical model embedded. The soil layer is Lambertian, and
it is assumed that its brightness at a specific waveband is a linear combination
of a number of spectral basis functions. By default, these are the basis
functions defined in Price 1991 , but they can be easily changed by users.
Note that the model will oversimplify the directional effects of typical rough
land surfaces, and that the Price spectral basis functions really only apply to
soils at field capacity. Also, the effect of e.g. snow is not considered.

The leaf model is the widely used
PROSPECT model, which treats leaves as a stack of thin plates. Other parameters
control the absorption of energy at different wavelengths. The parameters of
the PROSPECT model are

#. :math:`N`, the number of leaf layers,
#. :math:`C_{ab}`, the concentration of chlorophyll a+b :math:`[\mu g cm^{-2}]`,
#. :math:`C_{car}`, the concentration of carotenoids :math:`[\mu g cm^{-2}]`,
#. :math:`C_{w}`, equivalent leaf water :math:`[cm]`,
#. :math:`C_{dm}`, dry matter :math:`[\mu g cm^{-2}]`,
#. :math:`C_{sen}`, the proportion of senescent material (fractional, between 0 and 1).

The canopy model chosen is Semidiscrete by `Gobron et al. 1997 <http://www.google.co.uk/url?sa=t&source=web&cd=2&ved=0CCMQFjAB&url=http%3A%2F%2Ffapar.jrc.ec.europa.eu%2FWWW%2FData%2FPages%2FFAPAR_Software%2FFAPAR_SOFTWARE%2FNADIM%2FSemidiscret_96JD04013.pdf&rct=j&q=semidiscrete%20gobron&ei=YxRcTq7HK8bG8QP9vJGZDA&usg=AFQjCNElypSrP7JeUMxbkGMK01LaDgTqnA&sig2=p5OSBuf_Pl0Nxl21MWVgTw&cad=rja>`_, for which an adjoint has been develoepd here to allow more rapid state estimation. 
This model assumes a continuous canopy
and is governed by LAI, leaf area distribution (discretised into five different
classes), two terms that control the so-called "hotspot" effect, and three 
spectral terms (one for the soil and two for leaves), which are fed in from 
PROSPECT and the spectral soil model. In total, the state vector has thirteen or
fourteen components (four related to the soil, six or seven to the leaf 
properties and three to the canopy). This allows the simulation of the 
directional reflectance for narrow spectral wavebands and for a specified 
illumination/acquisition geometry.


Solving for a single observation
================================
We first apply the eoldas to the problem of attempting a state estimate from
a single satellite observation. These data have multiple wavebands, but we consider here only a single location and date.

An example data file for this is some MERIS data over a field site in Germany.
The file is ``MERIS_WW_1_A_1.brf``, the start of which looks like:

.. literalinclude:: ../data/meris/MERIS_WW_1_A_1.brf
  :lines: 1-11

We build a suitable configuration file ``config_files/meris_single.conf``

.. literalinclude:: ../config_files/meris_single.conf

You can note several things about the semidiscrete observation oeprator from this:

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 14-14

The (apparent) state variables, other than gamma_time that is (potentially) used for the differential operator are:

#. xlai : :math:`LAI`, the single sided leaf area per unit ground area
#. xhc    : The canopy height
#. rpl    : The leaf radius/dimension (same units as canopy height)
#. xkab   : :math:`C_{ab}`, the concentration of chlorophyll a+b :math:`[\mu g cm^{-2}]`,
#. scen   : :math:`C_{sen}`, the proportion of senescent material (fractional, between 0 and 1).
#. xkw    : :math:`C_{w}`, equivalent leaf water :math:`[cm]`,
#. xkm    : :math:`C_{dm}`, dry matter :math:`[\mu g cm^{-2}]`,
#. xleafn : :math:`N`, the number of leaf layers,
#. xs1    : Soil PC1 (soil brightness)
#. xs2    : Soil wetness
#. xs3    : not used
#. xs4    : not used
#. lad    : leaf angle distribution (coded 0-5). See Semidiscrete code for details.

However, the section ``[parameter.x.assoc_transform]`` contains:

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 40-41

which sets default state transformations and inverse transforms to be identity (e.g. xlai = xlai).

Then the section:

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 43-47

overrides this for the specifc states mentioned to
tell the eoldas that there is a transformation (in the RT model code) between true state
(e.g. :math:`LAI`) and the state we solve for. The reason for using a transformation is
in an attempt to approximately linearise the sensitivity of :math:`H(x)` to the state variable.
This is a good idea generally, for avoiding getting trapped in local minima when optimising
and also because Gaussian statistics are more appropriate for the linearised state.
The inverse transforms are covered in the sections ``[parameter.x.assoc_invtransform]``

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 49-53

The state vectors will be saved to file in their transformed form (as this is appropriate for the
stored uncertainty information), but plots are generated for the inverse transformed states.
The user has access to these inverse transforms in the code as x.transform and x.invtransform.

Second, we note that some of the state variables are 'switched off' for optimisation in the section ``parameter.solve``. Initially, all solve states are set to 1:

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 15-15

Then we switch selected states to 0:

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 23-29

In the ``[general]`` section, we have set the flags:

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 93-95

which tell eoldas to generate plots, and to keep updating them every 'plotmod' calls to the operator.
This allows the user to visualise how the data assimilation is proceeding, by viewing the appropriate png files.
The file names are formed from any results.filename fields found in a state.

Finally, we note the section ``[operator.obs.rt_model]``. This tells the eoldas parameters specific to
the operator being used. In the case of an observation operator, this must include ``operator.obs.rt_model.model``. Here, we have:

.. literalinclude:: ../config_files/meris_single.conf
  :lines: 104-111

where we set the model to ``semidiscrete1``, which will have been installed with the eoldas package (it is a dependency)

We can run this configuration with::

    ~/.local/bin/eoldas_run.py --conf=config_files/eoldas_config.conf \
      --conf=config_files/meris_single.conf \
      --parameter.limits='[[232,232,1]]' --calc_posterior_unc


which tells eoldas to solve for 8 state varibles: ``xlai, xkab, scen, xkw, xkm, xleafn, xs1, xs2`` from
any meris data found on day 232.

Here, since we are using an observation on a single data, we do not include the dynamic model, and
merely attempt to solve with the observation operator.
This process, with a single operator (the radiative transfer model), can also be called 'model inversion' where we
attempt to estimate the state vector elements known (and expressed in the radiative transfer model) to control
the spectral bidirectional reflectance of vegetation canopies). In truth we cannot claim that this is the
optimal solution based only on the model and observations, as we also have the implicit constraint of bounds to
the state estimate. This is a form of prior information that we use to help constrain the problem. Using an Identity Operator
as we have done in previous exercises is, in many ways just another way of expressing this same idea, just with a different
statistical distribution. With 'hard' bounds to the state vector elements we are assuming a 'flat' probability of
occurrence between those bounds. If we use an Identity Operator and describe :math:`y` as a Gaussian distribution, we have 
a 'soft' bound to the problem. In many cases there are physical limits to states, for example, xlai must lie in the bound [0,1]
or more practically away from the edges of those bounds, (0,1). Similarly, plant height cannot physically be negative.

Results are written to ``output/meris/MERIS_WW_1_A_1.params``:

.. literalinclude:: ../output/meris/MERIS_WW_1_A_1.params

with the foward modelled results in ``output/meris/MERIS_WW_1_A_1.fwd``:

.. literalinclude:: ../output/meris/MERIS_WW_1_A_1.fwd

and the original data in ``output/meris/MERIS_WW_1_A_1.fwd_orig``:

.. literalinclude:: ../output/meris/MERIS_WW_1_A_1.fwd_orig

with graphics:

.. image:: ../output/meris/MERIS_WW_1_A_1.fwd.plot.y.png


which is a plot of MERIS reflectance as a function of wavelength. The 
red dots show the orginal data points with assumed uncertainty shown 
as 95% confidence intervals (0.01 * 1.96 here). The green line is the 
modelled version, with the error bars indicating the 95% confidence 
intervals for each waveband (slightly, but not greatly lower than the 
prior uncertainties).

This shows that the eoldas is able to replicate the observed 
reflectance very well. If we can do this, then why do we need data 
assimilation? The answer to that comes from looking back at the 
solved state information.

For ``xlai`` (i.e. :math:`\exp(-LAI/2)`) we retrieve a value of 
0.721891, i.e. LAI of 0.652. This might well be plausible for
the particular location, but we see that we have an uncertainty of 
0.654474 on this. Since xlai is bounded [0,1] this means that this is 
a very unreliable estimate. 

For some of the model states that we requested to be solved
we do even worse, ``xkw`` (also bounded [0,1]) has an apparent 
uncertainty of 46.326810, though the reason for this 
is simply that the MERIS instrument does not sample at wavelengths 
at which reflectance is very sensitive to leaf water content. 
Interestingly, one of the lower uncertainties is 
0.166161 for ``xkab``. This is because MERIS *does* sample at several 
wavebands that are sensitive to chlorophyll content. The retrieved 
value here is 0.780574 which corresponds to a concentration of 
24.8 :math:`[\mu g cm^{-2}]`.

So, not only is our knowledge of these state variables we require 
actually quite poor from a single observation, but also, if we were to 
look at the full inverse Hessian, we would see high degrees of 
correlation between the state estimates.

This is not surprising: we are trying here to retrieve  8 states 
from measurements in 15 wavebands. We have assumed here an uncertainty
of 0.01 in the surface reflectance, which could be a little high for 
some bands, but note that we have assumed the error independent for 
each wavelength, which it most certainly would not be (in that case,
and in the absence of any further information on the error correlation
structure, it may be appropriate to inflate the standard deviation).

On the positive side, the optimisation did at least converge to a 
viable solution, and gives a final value of :math:`J` of around 3.83.
Remember that :math:`J` is essentially half the sum of the squared 
differences between modelled and measured reflectance, over 15 
wavebands here, gives an MSE equivalent of 0.51 and a RMSE of 0.71,
which we can loosely interpret as meaning that the 'average' departure 
between what we measure and what we simulate is around 0.71 standard 
deviations.

Simulating surface reflectance with the EOLDAS model operator
==============================================================

Even though the state estimate is a little poor, we can take it and simulate what we would see with
other sensors (i.e. simulate different sensors). This demonstrates how the state vector, phrased
as biophysical variables, allows us to translate information from one sensor to another (and thus
allows us to combine information from multiple sensors as we shall see).

We build a configuration file for the MODIS instrument as ``config_files/modis_single.conf``.

This is mainly the same as the meris configuration file other that a few items refering to
the input state (which is the output of the meris run here: ``output/meris/MERIS_WW_1_A_1.params``),
the band names and uncertainty, and the input observations file, ``data/brdf_WW_1_A_1.kernelFiltered.dat``.

.. literalinclude:: ../config_files/modis_single.conf

We can, for instance then run this with::

    
    eoldas_run.py --conf=config_files/eoldas_config.conf \
	--conf=config_files/meris_single.conf \
	--parameter.limits='[[232,232,1]]' \
	--passer --conf=config_files/modis_single.conf


We need to load the MERIS configuration file before the MODIS one so 
that the latter overrides and options in the former.  There are some 
dangers to using a whole string of configuration files as it can be 
difficult to keep track of information, but the software allows you 
to define as many as you want (up to any system limits to command 
lines etc.). 

You can confirm that the configuration file loading is as you would 
expect by looking at the log file ``logs/modis_single.log``.  You 
really should look at that the first time you generate a new 
configuration file or you might end up mis-interpreting results. 

Alternative to writing a new configuration file, we could have 
overridden the meris configuration settings with command line options.
That has some of the same issues as having too many configuration 
files though as it could become difficult to keep track of information.

The flag ``--passer`` here tells eoldas *not* to perform parameter 
estimation, but just to do a 'forward' simulation from the state that 
we load.

The result we are interested in is written to ``output/modis/MODIS_WW_1_A_1.fwd`` and ``output/modis/MODIS_WW_1_A_1.fwd_orig``:

.. literalinclude:: ../output/modis/MODIS_WW_1_A_1.fwd

.. literalinclude:: ../output/modis/MODIS_WW_1_A_1.fwd_orig

with the plot in 

.. image:: ../output/modis/MODIS_WW_1_A_1.fwd.plot.y.png

The result is rather poor for both of the MODIS observations on day 
232, especially at longer wavelengths. This could be due to factors 
such as poor atmospheric modelling of either the MODIS or MERIS 
datasets, but it could equally be just due to the large uncertainty in
the state estimates obtained from the MERIS data alone.

So, the information content of the MERIS data alone for a single 
observation only poorly constrain the estimate of the state variables.
If we had some prior expactation of those states, we could use that to
help constrain the estimate and we will proceed onto that in the next 
section.

In the meantime, it is instructive to see how we can use *both* the 
MERIS and MODIS observations to constrain the state estimate for this 
single day dataset.

Before that, to get an estimate from MODIS alone we create a modified 
form of ``config_files/modis_single.conf``, ``config_files/modis_single_a.conf``:

.. literalinclude:: ../config_files/modis_single_a.conf


We have changed the output filenames so as not to overwrite the previous results and remove the state definition for 
parameter.x (i.e. we don't use the information from MERIS in this example).

Similarly to before, we run::

  eoldas_run.py --conf=config_files/eoldas_config.conf \
    --conf=config_files/meris_single.conf --passer \
    --parameter.limits='[[232,232,1]]'  \
    --conf=config_files/modis_single_a.conf

   


The state estimate result is in ``output/modis/MODIS_WW_1_A_1.params_a`` 

.. literalinclude:: ../output/modis/MODIS_WW_1_A_1.params_a

with the graphic in

.. image:: ../output/modis/MODIS_WW_1_A_1.fwd_a.plot.y.png

This result is not particularly good, and if we examine the  s of the optimisation we will see a warning: ABNORMAL_TERMINATION_IN_LNSRCH. This is simly not a viable result, even though we have wider spectral sampling than for MERIS. That may be partly down to
particular characteristoics of these samples, or it might be for example that there just isn't enough information to solve the
problem (8 parameters from 7 wavebands, remember). 

If we suspected it was just that the solution was trapped in some local minimum because of an inappropriate 
starting point for the optimisation, we could include e.g. the lines:

.. literalinclude:: ../config_files/modis_single_a.conf
  :lines: 5-6

in the configuration file, as the optimisation would then start from the state found from the MERIS data, but in this case
that is of little value.

So, we can 'fit' to the MERIS data for a single date very well, but the resultant state estimates have very high uncertainty,
for some state elements simply because of poor spectral sampling. And we have some MODIS data with wider spectral sampling
(though fewer wavebands) but we cannot make direct use of it, even with two observations on a single day.

The obvious thing to try then is to combine the observational datasets, and this is something we can do very easily
in eoldas. Even though the spectral sampling of the data is very different, we can just define different observation operators
(using the same radiative transfer model) to deal with this. Of course, there are other issues to consider, such as
spatial sampling/location, but the same principle applies there in a general DA system.

We can base the new configuration on that we used for MERIS, and simply declare the MODIS operator with a different name:
(confs/modis_single_b.conf)

.. literalinclude:: ../config_files/modis_single_b.conf

We can then run with both meris and modis data::

    eoldas_run.py --conf=config_files/eoldas_config.conf \
      --conf=config_files/meris_single.conf \
      --parameter.limits='[[232,232,1]]' \
      --conf=config_files/modis_single_b.conf



In any case,  the result still does not have a good convergence, so the combination has not been a complete success.

The new (combined) state estimate is in ``output/modis/MODIS_MODIS_WW_1_A_1.params_b``:

.. literalinclude:: ../output/modis/MODIS_MODIS_WW_1_A_1.params_b

The new ``xlai`` value is 0.435037 with a standard deviation of 0.000730. This corresponds to a LAI of 1.66 (confidence interval 1.58 to 
1.66). The ``xkab`` value though is 0.958366 
which is very low (4.25) when inverse transformed. The apparent uncertainty 
for this (sd) is 0.018995 on ``xlai``, so the ``xlai`` confidence interval is 
0.0 to 6.25, since it is bounded at 0. These posterior estimates are 
probably unreliable: we can suggest that because the optimisation 
routine reported that  it did not find a very well converged solution. 

A Coding note and caveat
------------------------

.. warning::

   This needs updating!
   
Note that in the configuration file we have specified the operator.obs_modis.rt_model.model field as rtmodel_ad_trans1, rather than rtmodel_ad_trans2. 

In reality, these are the same models, but because they are FORTRAN shared object libraries containing the semidiscrete model
and because they make use of global variables, one cannot run more than one observation operator that call the same shared object.
To get around this problem, we have essentially made multiple copies of the library, so that different observation operators
(e.g. using data from different sensors) can each have their own version. New operators in the future should be written with more
careful consideration of such issues. Similar coding issues also make it next to impossible to run
multi-core paralellisation optimisations to speed up code execution. The lesson is ... avoid globals like the plague, but there 
is normally some way we can get around the problem. In this case, the implication is that only 8 different observation operators
calling the semidiscrete code can be run in one data assimilation exercise, as we have only generate 8 copies of the library.

Other forms of constraint
-------------------------
We have seen that we can apply eoldas to state vector estimation using the semidiscrete canopy reflectance model. This is 
very much what we want to be able to do from EO data: infer characteristics of the land surface, so that is positive.
However, observations from a single sensor very often have insufficient information content to reliably allow
an estimate of the range of state variables we know to affect the signal. The only constraint we have applied here beyond
the observations (and the model) is that of physical bounds to some of the states (N.B. in eoldas, if you want to set only
an upper or lower bound, you can use None to indicate that no bound should be set). 

Even though we managed a 'stable' retrieval from MERIS data, the resultant uncertainties were very high, which makes the
value of the retrieval rather low for any practical purposes. We saw that in this case, we could not find a viable solution
for the 2 MODIS samples we used. Also, we demonstrated how to use the observation operator can be used to map from
one EO measurement set (the 15 MERIS wavebands) to anoth (MODIS). The comparison was poor however because of the poorly
constrained estimate of the states from the MERIS data.

How then can we improve the situation? The only real answers are: (i) get sensor data with higher information content
(there are practical, physical and financial limitations to that however); (ii) add other sources of information.
In generating global satellite products, for example, it is common practice to assume that some of the state variables
(e.g. those affecting leaf and soil reflectance) are *known*. This is rather a dangerous assumption and it would be better
to at least phrase this supposed insight into some distribution (applying hard physical bounds which we can implement
in eoldas  by changing the limits of validity of the states; or assuming some statistical distribution`)

To demonstrate this, we build a new configuration file, ``config_files/modis_single_c.conf``:

.. literalinclude:: ../config_files/modis_single_c.conf

The only new sections are:

.. literalinclude:: ../config_files/modis_single_c.conf
  :lines: 9-27

which state that the operator (called ``operator.prior``) is to have 
:math:`H(x) = I(x)` (we call the operator Operator, the default)
and we, simplistically here, set the state and standard deviation to 
0.5. This then is a quite loose constraint as most states are bounded [0,1].

We run now with::

    
    eoldas_run.py --conf=config_files/eoldas_config.conf \
      --conf=config_files/meris_single.conf \
      --parameter.limits='[[232,232,1]]' \
      --conf=config_files/modis_single_c.conf

Now, we get convergence, with the results in ``output/modis/MODIS_WW_1_A_1.params_c``:

.. literalinclude:: ../output/modis/MODIS_MODIS_WW_1_A_1.params_c

One particular thing to note here is that *all* of the uncertainties are now less that 0.5. This happens (simply?) because
 we set the prior uncertainty to 0.5. 
 
.. note::

  Think about why that happens, statistically... 
  essentially because :math:`C_{post}^{-1} = C_{a}^{-1} + C_{prior}^{-1}`.

That is fine, and, as we have noted, expected (indeed we can think of that *being* the main purpose of applying the prior constraint), but
it does mean that users of eoldas must bear in mind the consequences of setting too tight constraints as they will result in
*apparently* low posterior uncertainties.

The original and forward modelled MODIS data resulting from this are:

.. image:: ../output/modis/MODIS_WW_1_A_1.fwd_c.plot.y.png

and the MERIS:

.. image:: ../output/meris/MERIS_WW_1_A_1.fwd_c.plot.y.png 

The reproduction of the observations is clearly not wonderful, indicating some other issues
with the data or model, but at least we have demonstrated how a prior can be used in eoldas to help
condition the solution. As indicated above, one has to be careful applying prior constraints, as the results
can be very dependent on the values chosen, but really such care must be taken about all choices here (which model, which data etc.).

Summary
-------

In this section, we have moved from applying simple operators where the domain of the state vector and observations was the same
to a more typical and useful case (for EO interpretation and use), where they are not. We have introduced a physically-based radiative transfer
tool as the operator function (one that is provided with an adjoint for (relatively) rapid optimisation of the combined cost functions :math:`J`.

We have shown how eoldas configuration files can be progressively modified (and used in combination) to set up particular problems.
We have obttained an estimate of biophysical parameters from MODIS data (failed to converge for observations on a single date), from
MERIS data, and from both MODIS and MERIS combined.

This is already quite a powerful tool, but we will see in the next section that other constraints can be interesting as well.



