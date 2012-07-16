***************************************************************************
Introduction
***************************************************************************

This webiste contains a **live** version of the EOLDAS documentation. 
EOLDAS stands for Earth Observation Land Data Assimilation System, and 
in here we refer to an implementation of a Data Assimilation (DA) 
prototype to be used to monitor the land surface in the optical domain.
The code is fairly generic, and can be used to implement a number of
different variational DA schemes (although we mostly have used weak constraints)

This software was developed

.. note:: 

   Partners, project number, thanks to etc.


Installation
^^^^^^^^^^^^^

EOLDAS has mostly been written in Python, with some radiative transfer
codes written in Fortran and made available to the main Python library.
The requirements for the package are

* Python (versions > 2.5 and < 3.0)
* Scipy
* Numpy
* Matplotlib
* gfortran
* OpenOpt (optional)
* git (version control system. Optional)

The first few packages are usually available in modern Linux distributions,
as well as in MacOsX. In some cases, the use of the `Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_
to install a whole Python ecosystem might be a worthy option. Similar 
efforts for Windows might be found in `Python xy <http://code.google.com/p/pythonxy/>`_.

OpenOpt is an optimisation suite, and is optional, although its use is 
highly recommended. Using ``pip`` or ``easy_install``, the Python
package installers, it can be installed as ::

easy_install -U openopt 

or ::

pip install openopt

Note that if you do not have root/superuser access, you can install it
for your user as ::

easy_install --user -U openopt

or ::

pip install --user openopt

The eoldas distribution can be installed using the same approach as
OpenOpt: ::

easy_install --user -U eoldas

or ::

pip install --user eoldas

This will install as an user, ignore the ``--user`` option to do a 
system-wide install. Note that executables will be installed too, and 
the installation executables path will be installation dependent. Note
that this command will also install the semidiscrete radiative transfer
code and compile it (many warning will show up, but they are safe to
ignore).

If you plan to do development on EOLDAS, we strongly encourage you to
`fork <https://help.github.com/articles/fork-a-repo>`_ the `EOLDAS repository <https://github.com/jgomezdans/eoldas>`_
on github, and feel free to submit patches, suggestions, etc.


Installing the user's guide
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The user's guide is available on  `github.com <http://github.com/jgomezdans/eoldas_release/>`_. 
Although you can just download an archive with the examples, data and
source of this guide, we encourage you to `fork <https://help.github.com/articles/fork-a-repo>`_ 
the `eoldas_release git repository <http://github.com/jgomezdans/eoldas_release/>`_,
 and to work on it and submit bug reports, fixes etc.

The EOLDAS community
^^^^^^^^^^^^^^^^^^^^^^^

EOLDAS is a complex and varied tool that can help solve many problems 
where Gaussian statistics are a good assumption in variational problems.
We would like users to use github as a way to collaborate in fixing bugs,
spreading what you can do with EOLDAS, and to put together tutorials.