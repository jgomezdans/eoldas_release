***************************************************************************
Inverting a time series of MODIS observations over agricultural fields
***************************************************************************

Introduction
^^^^^^^^^^^^^^^^

The current example builds on the previous synthetic example by using
data from the MODIS sensor on the TERRA and AQUA platforms to invert
the state of the land surface over an agricultural area in Thuringia
(Germany). We shall use a script similar to `sentinel.py <https://github.com/jgomezdans/eoldas_release/blob/master/sentinel.py>`_
but will use a slightly different inversion strategy. The data are in
the file `brdf_WW_1_A_1.kernelFiltered.dat <https://github.com/jgomezdans/eoldas_release/blob/master/data/brdf_WW_1_A_1.kernelFiltered.dat>`_

The solution strategy tries to overcome one of the weaknesses of the
system so far: the need for very costly cross-validation in order
to estimate the hyperparameters
