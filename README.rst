====================
EOLDAS examples
====================

.. contents:: :local:
    
    
Installation
=============

This repository holds the code and ancillary files that are required to run the
examples in the `EOLDAS user's guide <http://jgomezdans.github.com/eoldas_relase/`_. 
Please visit that website to find out more about EOLDAS!

However, we first need to install EOLDAS and any other modules that are required
(in this case, RT codes that are distributed as independent packages). This is
straightforward, as all these packages are python packages that can be installed
using either ``easy_install`` or ``pip``. The installation will download and 
install both the eoldas package, and the required RT packages. These will be 
compiled and installed (there will be lots of compiler warning errors. These are
normal and you should ignore them). To install using ``easy_install`` do ::
    
easy_install -U eoldas --user

    
This will install things in ``~/.local/lib`` and ``~/.local/bin``. 
The latter might be added to your ``$PATH`` to gain access to the main 
executable file, ``eoldas_run.py``.

Installation using pip is also possible ::

pip install eoldas

Note that a number of auxiliary packages are required. These are 
usually available from a complete python installation and are:

* numpy
* scipy
* matplotlib

The optimisation package `OpenOpt <http://openopt.org/>`_ is not 
required, but allows more flexibility in the choice of optimisation
methods.
    
Installing the example files using git
****************************************

The examples from the user's guide have been revamped. All the files (config and ancillary data files) are in *this* github repository. Once the above are all installed, you can either download them or *clone the repository using git*. The latter method is recommended, so you can track changes etc.

Installing from the Zip file
--------------------------------

Click on the icon that says "ZIP" towards the top left of the page), unzip the resulting file somewhere, and run the examples from this new directory. 

Cloning the repository
--------------------------

Make sure you have a recent version of git installed, then do ::
    
git clone https://github.com/jgomezdans/eoldas_examples.git

If you plan to contribute to the codebase, it's recommended that you fork the repository (see `here for more information <https://help.github.com/articles/fork-a-repo>`_ 



Experiment 1 ( Savitzky-Golay smoothing)
==============================================

The first example (in  `here <http://www2.geog.ucl.ac.uk/~plewis/eoldas/eoldas_guide.html>'_) runs a Savitzky-Golay filter over some MODIS NDVI data. The script is ``savitzy_golay.py``, and it requires the file ``data/FuentesAndalucia_MOD09A1.txt`` that holds the data. You run the script by issuing the command ``python savitzky_golay.py``. This will pop a window, that when closed will save the plot to the images directory (this directory will be created if it doesn't exist). This plot should correspond to `this one <http://www2.geog.ucl.ac.uk/~plewis/eoldas/_images/golay.png>`_

Experiment 2 (EOLDAS with an identity operator)
====================================================

This experiment is described in depth in `here <http://www2.geog.ucl.ac.uk/~plewis/eoldas/example1.html>`_. It requires two configuration files, stored in config_files: ``Identity.conf`` and ``eoldas_config.conf``. It also requires a datafile in ``data/Identity/random_ndvi1.dat``. The script that uses these files is ``solve_eoldas_identity.py``. The script is executed with the command ``python solve_eoldas_identity.py``. After running this, the script will create a number of directories to store logs (``mylogs/``), solutions and further diagnostics (these will be under ``output/Identity/``). The results should correspond to figures in the user's guide.

An extra script is provided in `here <http://www2.geog.ucl.ac.uk/~plewis/eoldas/example1.html#example-plotting-data-from-the-output-files>`_. This shows how to make plots of the output using Python and Matplotlib. The required script is ``example1plot.py``. This script will produce a plot stored under the ``images`` directory. It should be compared to the one on the user's manual. Similarly for ``solve_eoldas_identity1.py``, you can use ``example1plot1.py`` to plot figures simular to those in the user's guide. Similarly, you can try ``solve_eoldas_identity2.py`` and ``example1plot2.py``.

To plot the Hessian, as in `here <http://www2.geog.ucl.ac.uk/~plewis/eoldas/example1.html#interfacing-a-little-more-deeply-with-the-eoldas-code>`_, you can use the script ``solve_eoldas_identity_a.py``. It will save the plot in ``output/IHessianNDVI_expt1.png``.

Experiment 3 (smoothing of MODIS observations)
=====================================================

.. note::
    
   The command line is slightly changed from that in the user's guide  (assuming you have  ``~/.local/bin/`` in your ``$PATH``. If that's not the case, prepend it to the command line)
    
    
The `third experiment <http://www2.geog.ucl.ac.uk/~plewis/eoldas/example1a.html>`_ requires some real MODIS observations. These are given in ``data/modis_botswana.dat``. You will also require the configuration file ``config_files/Identity2.conf``. To produce all the plots in this section, you will need to have the ``eoldas_run.py`` executable in your path (this will usually be ``~/.local/bin/``). ::

    eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/Identity2.conf --calc_posterior_unc

The output of this experiment will come out in ``output/Identity`` (with names such as ``Botswana_fwd.params2.plot.y.png`` or ``NDVI_Identity1.params.plot.x.png``)
To run the other examples, do ::

    eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/Identity2.conf --calc_posterior_unc --operator.modelt.rt_model.model_order=2 --parameter.x.default=5000,0.1 --operator.obs.y.result.filename=output/Identity/Botswana_fwd.params2 --parameter.result.filename=output/Identity/MODIS_botswana.params2
    eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/Identity2.conf --calc_posterior_unc --operator.modelt.rt_model.model_order=2 --parameter.x.default=200,0.1 --operator.obs.y.result.filename=output/Identity/Botswana
    
Check the actual command line to see where output results are going to come out.

Experiment 4 (RT observation operators )
================================================

Radiative transfer modelling for optical remote sensing. In `this experiment <http://www2.geog.ucl.ac.uk/~plewis/eoldas/example2.html>`_, we will use the semidiscrete model to invert and forward model real observations from spaceborne sensors. The first experiment gets a single observation from MERIS (15 bands in the visible/near-infrared range), and inverts this observation. The command to run it is: ::

~/.local/bin/eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/meris_single.conf --parameter.limits='[[232,232,1]]' --calc_posterior_unc

The solution will appear in ``output/meris/``, where you can find both the text files and plots that are in the users' guide.

A second example uses the results from the first, and uses the estimated state of the land surface to provide a prediction of the reflectance that would be seen by the MODIS sensor on that same day. This is then compared to the actual observations. The command is ::

~/.local/bin/eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/meris_single.conf --parameter.limits='[[232,232,1]]' --passer --conf=config_files/modis_single.conf 

Other experiments in that section are: ::

~/.local/bin/eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/meris_single.conf --parameter.limits='[[232,232,1]]' --passer --conf=config_files/modis_single_a.conf 

(the output for this will be in e.g. ``output/modis/MODIS_WW_1_A_1.fwd_a.plot.y.png``). The following experiment will be ::
    
    eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/meris_single.conf --parameter.limits='[[232,232,1]]' --conf=config_files/modis_single_b.conf 
    
Output for MERIS will be in eg ``output/meris/MERIS_WW_1_A_1.fwd_b.plot.y.png`` whereas for MODIS it will be in ``output/modis/MODIS_WW_1_A_1.fwd_b.plot.y.png``

The experiment that demonstrates changing the prior definition is run like ::

    eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/meris_single.conf --parameter.limits='[[232,232,1]]' --conf=config_files/modis_single_c.conf 
    

Experiment 5 (time series assimilation using observation operator )
===========================================================================

Experiment 5.1
******************

This experiment is the same as in the webpage, but with different initialisation methods, prior ranges etc. We hope to expand the user's guide with more comment. The experiment is contained in a single python script, ``modisDa0.py``, and can be executed by calling it as ``python modisDa0.py``.

Experiment 5.2
******************

The current experiment expands on the previous experiments by demonstrating how the smoothing constraint is applied to a problem where a complex observation operator is used to model the observations. The details of this experiment are `in this page <http://www2.geog.ucl.ac.uk/~plewis/eoldas/example3.html>`_ . The first example can be performed by the following command (it takes around 30-45 minutes)::
    
~/.local/bin/eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/semid_default.conf --logfile=logs/file.log --parameter.result.filename=output/gamma2/state.dat --operator.obs.y.result.filename=output/gamma2/obs.dat --parameter.x.default=2,0.99,5,0.01,0.99,0.001,0.99,0.35,1.5,1,0.001,0,0,5 --operator.obs.y.state=data/brdf_WW_1_A_1.kernelFiltered.dat

You will see results appearing in ``output/gamma2``. In particular, the state that we are solving for appears in ``output/gamma2/state.dat.plot.x.png`` (for the plot) and ``output/gamma2/state.dat`` for the raw data. Also note that in this example, the posterior covariance matrix and Hessian have not been calculated.

.. note::
   
   Here's where the hack to do several different gammas starting from a previous gamma needs to go.


Experiment 6
=============

This forms part of the Sentinel experiment from  Lewis et al. (2012).

To run: ::

    eoldas_run.py --conf=config_files/eoldas_config.conf --conf=config_files/sentinel.conf --logfile=logs/file1.log

In ``config_files/sentinel.conf``, we initialise state using the file ``data/rse1_init.dat`` (to spped things up here). The synthetic Sentinel-2 MSI observations are in the file `data/rse1_test.100.dat`.

You can view the results as they process in e.g. ``output/rse1/rse1_test.dat.plot.x.png``.

