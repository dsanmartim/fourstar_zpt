.. FourStar Cookbook documentation master file, created by
   sphinx-quickstart on Mon May 24 08:56:11 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to FourStar Cookbook's documentation!
=============================================

This is a guide on how to observe and reduce spectrophotometric data obtained with
LDSS-3 at Clay Telescope. In its first part (:ref:`photo`) it provides an overview
on how to observe photometric data and how to obtain zero-points from raw images.
In its second part (:ref:`spec`) it provides an overview on how to observe a spectroscopic
standard star and also show you how to use a python wrapper script that uses ``IRAF/PyRAF``
routines to reduce engineering data. The final product will be 1d spectra, mainly useful
to estimate throughput performance.

.. _photo_obs:
.. toctree::
   :maxdepth: 2
   :caption: Standard Star Observations:

   fourstar_obs

.. _photo_red:
.. toctree::
   :maxdepth: 2
   :caption: Standard Star Reduction:

   fourstar_red


.. toctree::
   :maxdepth: 2
   :caption: About the Cookbook
   :titlesonly:

   about


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
