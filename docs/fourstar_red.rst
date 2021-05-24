
Observation
===========

An extensive cookbook for most the the LDSS3 modes is available in the LCO web
page (`LDSS3 Observing Cookbook
<http://www.lco.cl/?epkb_post_type_1=observing-cookbooks-2>`_). There you will
find a detailed description on how to use the instrument GUI, how to obtain
a direct image, how to align a mask (long-slit and multi-object) and everything
else related to data acquisition. Here, we just provide some guidelines on
instrument setup and parameters required to observe a spectrophotometric
standard star that we use to keep tracking of the instrument performance.

Which Star or Field?
--------------------



Photometric Setup
-----------------

Unless a different configuration is requested, always use the setup as
defined below. For photometric zero points measurements, this is the
configuration we usually keep track of instrument performance.

Requested Setup:
  | **Binning**: 1x1
  | **Readout Speed**: Fast
  | **CCD Gain**: Low
  | **ROI**: Full Array (Amplifier 1)

Note that exposure times should not be lower than 7 seconds, since the shutter
has some limitations regarding linearity with shorter exposures. Note also that
we only need data from amplifier number one (c1) for now. The main reason
for that is explained in the next paragraph.

Since, for now, we are measuring the zero points by using a single star from the
field, make sure the star you are observing is not saturated. If the sky is
photometric, most probably, you will only be able to achieve that if the
telescope is defocused. For that, ask the telescope operator to defocus the
telescope by around 150-250 units. You will notice that the stars will look
like a donuts instead of having a 2d-gaussian shape, which is fine for our
purposes.

.. important::

   As we are measuring the zero points with a single star, we request
   that you do the *ltoslit* step described in the
   `LDSS Direct Imaging
   <http://www.lco.cl/?epkb_post_type_1=ldss-direct-imaging>`_
   section of the
   `LDSS3 Observing Cookbook
   <http://www.lco.cl/?epkb_post_type_1=observing-cookbooks-2>`_
   by using ``xcoord = 1500.0`` and ``ycoord = 2000.0`` as input to that IRAF
   routine.


We chose to put our star in the position stated above (which put the star
to fall in the amplifier #1) mainly because we want to use the same CCD
region that has been used for a long time in the observatory and also because
we are not doing any calibration before measuring the zero points. That is the
way it has been done and we don't want to change that until we implement
a way to measure the zero points by using all stars from the field.

Observing sequence
------------------

After proper field alignment, we recommend a sequence of three standard star
exposures for each filter installed on the instrument:

    | 3 x standard star (for each filter)

As we are not doing any data reduction before obtaining the zero points, we
do not require any calibration frame. Three images in each filter is everything
we need to obtain zero points.

Estimating the Zero Points
==========================

Overview
--------

If you access the GitLab repository associated with this cookbook, you will
see a folder there that is called ``tools``, which contains a python module
called ``zpt.py``. This module has a couple of methods that we use to 1)
identify sources in the image, 2) compute the aperture photometry of identified
sources, 3) display the results and 4) estimate the zero points of a single
source in the field. Those methods require some user intervention, although
we tried to keep it as minimum as possible. We made available inside the
``tools`` folder in the gitlab repository a jupyter notebook that explains
and guide the user to the methods made available in the ``zpt.py``.

If you have downloaded the entire repository to your computer, then you already
have the  ``zpt.ipynb`` jupyter  notebook in your system.
`Here <http://gitlab.lco.cl/dsanmartim/ldss/blob/master/zpt.ipynb>`_ (inside
the LCO VPN) you can see it and download it again if you want. So, assuming you
have a python3 environment with all the
`required packages <http://gitlab.lco.cl/dsanmartim/ldss/tree/master#requirements>`_,
and that you have activate your python 3 environment, open a jupyter notebook
session:

.. code:: bash

   $ (my_environment) jupyter notebook

Below, we reproduce the content of the jupyter notebook, so you can have easy
access to its content. However, the content shown here is the same that you will
find in the notebook file. Therefore, you can skip the sections below and go
directly to the ``zpt.ipynb`` notebook. Note that all the commands we show below
should be executed inside a jupyter notebook, although it should also work in
a normal (i)python session.

.. note::

   The documentation to the methods available in the ``zpt.py`` module can be
   found in the :ref:`api_doc` section. There, you will find some more details
   about each method and its input parameters.

