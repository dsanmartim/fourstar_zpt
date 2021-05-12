import numpy as np
from importlib_metadata import FileNotFoundError
import matplotlib.pyplot as plt
import pandas as pd

import seaborn as sns

from matplotlib.patches import Circle, Ellipse
from scipy.stats import norm


from astropy.io import fits
from astropy import wcs
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from astropy.coordinates import SkyCoord
from astroquery.irsa import Irsa
from photutils import DAOStarFinder
from photutils import CircularAperture, CircularAnnulus, aperture_photometry

import astropy.units as u
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_hdulist(filename):
    """
    It reads an fits file an return it as an astropy.io.fits.hdu.hdulist.HDUList object.

    Args:
        filename (str): The input is given as an string, but the file itself must contain a fits frame with a Primaray
        HDU and the corresponding data. Any file can be read, although a 2d image is expected.

    Returns (astropy.io.fits.hdu.hdulist.HDUList object): It returns an astropy HDUList.

    """
    hdu_list = fits.open(filename)

    ndim = hdu_list[0].data.ndim
    if ndim is not 2:
        raise TypeError('Data contains {} dimensions while it was '
                        'expected 2 dimensions.'.format(ndim))

    return hdu_list


def get_image(filename, ext=0):
    """
    Read a given extension of a mosaiced 2D FourStar image. It returns a 2D numpy array that contains the data and the
    corresponding image header.

    Args:

        filename (str): the file name of a FITS file.

        ext (int): extension of the header that will be read. As default it read the Primary HDU (ext=0).

    Returns:

        header (astropy.io.fits.Header) : An astropy header of the input fits image.

        data (numpy.ndarray) : A 2D numpy array that contains the data.

    """
    hdu = fits.open(filename)
    header = hdu[ext].header
    data = hdu[ext].data

    if data.ndim is not 2:
        raise TypeError('Data contains {} dimensions while it was '
                        'expected 2 dimensions.'.format(data.ndim))

    return header, data


def display(data, sources, saturated=None, radio=18.0, stretch='linear', percent=99.5, color_map='Greys'):
    """
    Display the input sources over plotted against the data array.

    Args:

        data (numpy.ndarray): A 2D numpy array that contains the data to be displayed.

        sources (pandas.DataFrame): a data frame with the source information

        saturated (padands.DataFrame or None): If an data frame with the saturated source is given as in input,
            then they are displayed. If None, saturated sources are not displayed.

        radio (int or float): radio of of the circles overplotted around the sources.

        stretch (str, optional): {'linear', 'sqrt', 'power', log', 'asinh'}. The stretch function to apply to the
            image. The default is 'linear'.

        percent (float): The percent scale. It works more or less as the ds9 percent display.

    """
    global ax

    '''
    data : `~numpy.ndarray`
        The image array.

    stretch : {'linear', 'sqrt', 'power', log', 'asinh'}, optional
        The stretch function to apply to the image.  The default is 'linear'
    '''

    if not isinstance(data, np.ndarray):
        raise TypeError('Please, use a np.array as input')

    if data.ndim is not 2:
        raise TypeError('Data contains {} dimensions while it was '
                        'expected 2 dimensions.'.format(data.ndim))

    if len(sources) < 1:
        print('No sources provided as input. The source input data frame is empty.')

    if len(sources) >= 1:

        fig, ax = plt.subplots(1, 1, figsize=(9, 9))

        # Getting image normalization
        norm_image = simple_norm(data, stretch=stretch, percent=percent)

        # Plotting image
        im = ax.imshow(data, origin='lower', norm=norm_image, cmap=color_map, interpolation='nearest')  #  cmap='viridis', Greys

        # Plotting circles and lables around sources
        positions = np.transpose((sources.xcentroid.values, sources.ycentroid.values))
        apertures = CircularAperture(positions, r=radio)
        apertures.plot(color='blue', lw=1.5, alpha=0.3)

        # Labels: using matplotlib
        xc, yc = list(sources.xcentroid), list(sources.ycentroid)

        for i, source_id in enumerate(list(sources['id'])):
            ax.annotate(source_id, (xc[i] + radio + 5.0, yc[i] + radio + 5.0), fontsize=9, color='blue')

        # create an axes on the right side of ax. The width of cax will be 4%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.06)

        # Define the colorbar and fix the labels
        cbar = plt.colorbar(im, cax=cax)

        # Define labels
        cbar.set_label(r'Arbitraty Flux Count Rate')

        ax.set_xlabel('X (pixels)')
        ax.set_ylabel('Y (pixels)')
        ax.set_title('Identified Sources')

    if saturated is not None:

        sat_pos = np.transpose((saturated.xcentroid.values, saturated.ycentroid.values))
        sat_ap = CircularAperture(sat_pos, r=radio)
        sat_ap.plot(color='red', lw=1.5, alpha=0.3)

        # Labels: using matplotlib
        sat_xc, sat_yc = list(saturated.xcentroid), list(saturated.ycentroid)

        for i, sat_id in enumerate(list(saturated['id'])):
            ax.annotate(sat_id, (sat_xc[i] + 25.0, sat_yc[i] + 25.0), fontsize=9, color='red')


def find_stars(data, fwhm=22, threshold=5.0, saturation=65000.0, round_low=-0.2, round_high=0.2, show=None):
    """
    It runs the DAOStarFinder ovef the data to identify star-like sources.

    Args:

        data (numpy.ndarray): A 2D numpy array that contains the data to be displayed.

        fwhm (float): The full-width half-maximum (FWHM) of the major axis a Gaussian kernel in
            units of pixels. If stars in the data array are defocused, the number of identified sources can
            change a lot depending of the fwhm input;

        threshold (float): Detection intensity threshold. Only stars above the sigma clipped (sigma=3) stdev of the data
            array will be detected;

        saturation (float): Sources with flux peak above the saturation level will be stored in a bad_sources dataframe
            and ruled out from the good_sources data_frame;

        show (str or None): If the input is not None, the identified sources will be display over the data image.
            If 'None' is given, no display is shown;

    Returns:

        df (pandas.Dataframe): a dataframe with all information of all potential sources identified.

        good_sources (pandas.Dataframe): a dataframe with only the good sources, which will exclude source with a
        a peak flux greater than the saturation level or with negative value (which is a fake source).

        good_sources (pandas.Dataframe): a dataframe with the saturated sources (stars) and with negative values.
    """

    if not isinstance(data, np.ndarray):
        raise TypeError('Please, use a np.array as input')

    if data.ndim is not 2:
        raise TypeError('Data contains {} dimensions while it was '
                        'expected 2 dimensions.'.format(data.ndim))

    # Getting sigma clipped statistics
    mean, median, std = sigma_clipped_stats(data, sigma=3.0, cenfunc='median', stdfunc='std')

    # Running DAOStarFinder
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold * std, peakmax=None, exclude_border=True,
                            roundlo=round_low, roundhi=round_high)

    sources = daofind(data - median)

    if len(sources) < 1:

        print('No source was found. Try the change the fwhm input or the threshold')

    else:

        for col in sources.colnames:
            sources[col].info.format = '%.8g'  # for consistent table output

        # Full data frame
        df = sources.to_pandas()

        # Get good sources and saturated ones
        good_sources = df[(df.peak < float(saturation)) & (df.peak > 0.0)]
        bad_sources = df[(df.peak >= float(saturation)) | (df.peak < 0.0)]

        if len(bad_sources) < 1:
            message = 'Good news. No saturated sources!'.format('%s')
            print(message)

        if show is not None:

            if len(bad_sources) >= 1:

                display(data, good_sources, bad_sources, radio=18.0, stretch='linear', percent=99.5, color_map='Greys')

            elif len(bad_sources) < 1:

                display(data, good_sources, saturated=None, radio=18.0, stretch='linear', percent=99.5, color_map='Greys')

        return df, good_sources, bad_sources


def load_wcs_from_file(filename, ext=0):
    """

    Args:
        ext:
        filename:

    Returns:

    """

    # Load the FITS hdulist using astropy.io.fits
    hdu_list = fits.open(filename)

    # Parse the WCS keywords in the primary HDU
    w = wcs.WCS(hdu_list[ext].header)

    return w


def get_wcs(hdu_list, ext=0):
    """

    Args:
        hdu_list:
        ext:

    Returns:

    """
    world = wcs.WCS(hdu_list[ext])

    return world


def sky2world(hdu_list, sources, max_error=1e-6):
    """

    Args:
        hdu_list:
        sources:
        max_error:

    Returns:

    """
    # Getting wcs from header
    w = get_wcs(hdu_list, ext=0)

    # Copying sources to df data frame
    df = sources

    # Three pixel coordinates of interest.
    # Note we've silently assumed an NAXIS=2 image here.
    # The pixel coordinates are pairs of [X, Y].
    # The "origin" argument indicates whether the input coordinates
    # are 0-based (as in Numpy arrays) or
    # 1-based (as in the FITS convention, for example coordinates
    # coming from DS9).
    pix = np.transpose([df.xcentroid, df.ycentroid])

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 0-based (Numpy-like) coordinates.
    world = w.wcs_pix2world(pix, 0)

    # Convert the same coordinates back to pixel coordinates.
    pix2 = w.wcs_world2pix(world, 0)

    # These should be the same as the original pixel coordinates, modulo
    # some floating-point error.
    message = 'Errors are greater than  {}!'.format(max_error)
    assert np.max(np.abs(pix - pix2)) < np.float(max_error), message

    df['ra_deg'] = np.array(world[:, 0])
    df['dec_deg'] = np.array(world[:, 1])

    return df

def world2sky(hdu_list, sources, catalog, max_error=1e-6, xcol_name='xcentroid', ycol_name='ycentroid'):
    """

    Args:
        hdu_list:
        sources:
        catalog:
        max_error:
        xcol_name:
        ycol_name:

    Returns:

    """

    # Getting wcs from header
    w = get_wcs(hdu_list, ext=0)

    # Three pixel coordinates of interest.
    # Note we've silently assumed an NAXIS=2 image here.
    # The pixel coordinates are pairs of [X, Y].
    # The "origin" argument indicates whether the input coordinates
    # are 0-based (as in Numpy arrays) or
    # 1-based (as in the FITS convention, for example coordinates
    # coming from DS9).
    world = np.transpose([catalog.ra, catalog.dec])

    # Convert pixel coordinates to world coordinates
    # The second argument is "origin" -- in this case we're declaring we
    # have 0-based (Numpy-like) coordinates.
    pix = w.wcs_world2pix(world, 0)

    # Convert the same coordinates back to world coordinates.
    world2 = w.wcs_pix2world(pix, 0)

    # These should be the same as the original pixel coordinates, modulo
    # some floating-point error.
    message = 'Errors are greater than  {}!'.format(max_error)
    assert np.max(np.abs(world - world2)) < np.float(max_error), message

    modified_catalog = catalog

    modified_catalog[xcol_name] = np.array(pix[:, 0])
    modified_catalog[ycol_name] = np.array(pix[:, 1])
    modified_catalog['id'] = sources.id

    return modified_catalog


def get_distance_from_center(hdu_list, sources, center=None, scale=0.16):
    """

    Args:
        hdu_list:
        sources:
        center:
        scale:

    Returns:

    """
    updated_sources = sources

    if not isinstance(scale, (float, int, type(None))):
        raise TypeError('Please, the scale parameter muste a float, int or None type.')

    header = hdu_list[0].header

    # Getting scale from header or from explicit input
    dxy = np.float(scale) if scale is not None else np.float(header['scale'])

    # Getting xc, yc coordintes (center field) from header or explicitly
    xc = np.float(center[0]) if center is not None else np.float(header['CRPIX1'])
    yc = np.float(center[1]) if center is not None else np.float(header['CRPIX2'])

    # Computing radio in image
    radio = np.sqrt((sources.xcentroid - xc) ** 2 + (sources.ycentroid - yc) ** 2) * (dxy / 60.0)

    updated_sources['r_arcmin'] = radio

    return updated_sources

def get_coordinates(hdu_list, ra_deg=330.0, dec_deg=29.3, frame='icrs'):
    """
    hdu_list (astropy.io.fits.hdu.hdulist.HDUList or None):

    ra_deg (None or float): The RA can be either a float number (in degrees) or a keyword from header. Note
    that the RA in header must be in degree also.

    dec_deg (None or float): The RA can be either a float number (in degrees) or a keyword from header. Note
    that the RA in header must be in degree also.

    """
    global hdr

    if hdu_list is not None:
        # Getting header
        hdr = hdu_list[0].header
        data = hdu_list[0].data

        if data.ndim is not 2:
            raise TypeError('Data contains {} dimensions while it was '
                            'expected 2 dimensions.'.format(data.ndim))

    elif hdu_list is None:

        if not isinstance(ra_deg, (float, type(None))):
            raise TypeError("RA might be a None or a float number given in degrees (e.g. 333.333)'")

        if not isinstance(dec_deg, (float, type(None))):
            raise TypeError("DEC might be None or a float number given in degrees (e.g. -29.666)'")

    alpha = np.float(hdr['CRVAL1']) if hdu_list is not None else np.float(ra_deg)
    delta = np.float(hdr['CRVAL2']) if hdu_list is not None else np.float(dec_deg)

    # Getting coordintates object
    c = SkyCoord(ra=alpha * u.degree, dec=delta * u.degree, frame=frame)

    coordinates = c.to_string('hmsdms')

    return coordinates


def get_two_mass_sources(coordinates, spatial='Cone', radius=7.5 * u.arcmin):
    """
    Doc >> https://astroquery.readthedocs.io/en/latest/irsa/irsa.html

    #Irsa.print_catalogs()

    coordinates (str, `astropy.coordinates` object):

        Gives the position of the center of the cone or box if
        performing a cone or box search. The string can give coordinates
        in various coordinate systems, or the name of a source that will
        be resolved on the server (see `here
        <https://irsa.ipac.caltech.edu/search_help.html>`_ for more
        details). Required if spatial is ``'Cone'`` or ``'Box'``. Optional
        if spatial is ``'Polygon'``.

    spatial (str):

        Type of spatial query: ``'Cone'``, ``'Box'``, ``'Polygon'``, and
        ``'All-Sky'``. If missing then defaults to ``'Cone'``.

    radius (str or `~astropy.units.Quantity` object) [optional for spatial is ``'Cone'``]
        The string must be parsable by `~astropy.coordinates.Angle`. The
        appropriate `~astropy.units.Quantity` object from
        `astropy.units` may also be used. Defaults to 10 arcsec.

    """
    # Catalog description: https://old.ipac.caltech.edu/2mass/releases/allsky/doc/sec2_2a.html#gal_contam
    table = Irsa.query_region(coordinates, catalog='fp_psc', spatial=spatial, radius=radius)

    catalog_df = table.to_pandas()

    return catalog_df


def match_sources(sources, catalog, hdu_list, max_separation=5 * u.arcsec):
    """

    Args:
        sources:
        catalog:
        hdu_list:
        max_separation:

    Returns:

    """

    # Getting header information
    header = hdu_list[0].header
    scale = header['scale']  # arcsec / pixel

    # Getting max separation
    max_sep = (max_separation / u.arcsec) / scale

    # Getting wcs
    _wcs = get_wcs(hdu_list, ext=0)

    _df = pd.DataFrame([])
    for i, star in catalog.iterrows():

        # Getting (x,y) coordinates for stars from catalog
        star_radec = SkyCoord(star.ra, star.dec, unit=(u.degree, u.degree))
        star_xy = SkyCoord.to_pixel(star_radec, wcs=_wcs, origin=0)

        # (x,y) positinos
        x, y = star_xy[0].item(0), star_xy[1].item(0)

        for j, source in sources.iterrows():

            separation = np.sqrt((source.xcentroid - x) ** 2 + (source.ycentroid - y) ** 2)

            if separation <= max_sep:
                star['id'] = '{:3d}'.format(int(source.id))
                star['xcentroid'] = source.xcentroid
                star['ycentroid'] = source.ycentroid
                star['sharpness'] = source.sharpness
                star['roundness1'] = source.roundness1
                star['roundness2'] = source.roundness2
                star['npix'] = source.npix
                star['sky'] = source.sky
                star['peak'] = source.peak
                star['flux'] = source.flux
                star['mag'] = source.mag

                _df = _df.append(star)

    columns = np.concatenate([sources.columns.values, catalog.columns.values])
    _df = _df[columns]

    # Sorting data frame based on source id
    sdf = _df.sort_values(by='id', axis=0, ascending=True)

    if len(_df.values) > 0:
        print('{} 2MASS Point Sources macth with sources from image.'.format(len(_df.values)))
    else:
        print('No star from 2MASS Point Sources Catalog match with sources from image.')

    return sdf


def get_photo(data, sources, aperture=25.0, annulus_width=6.0, show=True, join_frames=True):
    """
    It compute the aperture photometry for the input data. The user is asked to provide a few input parameters
    regarding the source finding sources code and two other regarding the photometry itself. The parameters related
    to the source finding are named starting with a 'f' prefix and the ones related to the photometry, starting with
    a 'p' prefix.

    Args:
  
        sources (pandas.Dataframe): a pandas dataframe that contains the (x,y) position of the sources identified in the
            data image. Note that sources must have already been identified. This routine will only calculate the 
            aperture photometry. It is mandatory that the columns with x and y position are named as 'xcentroid' and
            'ycentroid', respectively. 
        
        data (numpy.ndarray): A 2D numpy array that contains the data to be displayed;

        f_fwhm (float): The full-width half-maximum (FWHM) of the major axis a Gaussian kernel in units of pixels. If
            stars in the data array are defocused, the number of identified sources can
            change a lot depending of the fwhm input.

        f_threshold (float): Detection intensity threshold. Only stars above the sigma clipped (sigma=3) stdev of the
            data array will be detected.

        f_saturation (float): Sources with flux peak above the saturation level will be stored in a bad_sources data
            frame and ruled out from the good_sources data frame.

        aperture (float): The radius of the circular aperture in pixels.

        annulus_width (int): The length of the annulus, whose internal radio is by default r_in = p_aperture + 10
            pixels. The external radio will be then r_out = r_in + p_annulus_width [pixels].

        show (str, bool or None): If the input is not None, the identified sources will be display over the data image.
            If None is given, no display is shown.

    Returns:

        photo_df (pandas.DataFrame): A data frame with the resulting aperture photometry. Optionally, if show is True,
        it shows the identified sources with the aperture and annulus from from photometry over plotted against the
        data.

    """

    if not isinstance(data, np.ndarray):
        raise TypeError('Please, use a np.array as input')

    if data.ndim is not 2:
        raise TypeError('Data contains {} dimensions while it was '
                        'expected 2 dimensions.'.format(data.ndim))

    if not isinstance(join_frames, (bool, type(None))):
        raise TypeError('Please, join_frames parameter muste be a bool or None type.')

    s = sources

    positions = np.transpose((s.xcentroid.values, s.ycentroid.values))
    apertures = CircularAperture(positions, r=aperture)

    r1 = aperture + 10.0
    r2 = r1 + annulus_width

    annulus_aperture = CircularAnnulus(positions, r_in=r1, r_out=r2)
    annulus_masks = annulus_aperture.to_mask(method='center')

    bkg_median = []
    for mask in annulus_masks:
        annulus_data = mask.multiply(data)
        annulus_data_1d = annulus_data[mask.data > 0]
        _, median_sigclip, _ = sigma_clipped_stats(annulus_data_1d)
        bkg_median.append(median_sigclip)

    bkg_median = np.array(bkg_median)

    photo = aperture_photometry(data, apertures)
    photo['id'] = s.id.values
    photo['annulus_median'] = bkg_median
    photo['aper_bkg'] = bkg_median * apertures.area
    photo['aper_sum_bkgsub'] = photo['aperture_sum'] - photo['aper_bkg']

    # Passing photo to pandas
    photo_df = photo.to_pandas()
    photo_df.style.format("{:.8g}")     # for consistent table output

    if show is not None:

        norm = simple_norm(data, 'linear', percent=99.5)

        fig, ax = plt.subplots(1, 1, figsize=(9, 9))
        im = plt.imshow(data, norm=norm, origin='lower', cmap='Greys', interpolation='nearest')  # cmap='viridis', 'Greys'
        apertures.plot(color='red', lw=1)
        annulus_aperture.plot(color='green', lw=1)

        xc, yc = list(s.xcentroid), list(s.ycentroid)

        for i, source_id in enumerate(list(photo_df.id)):
            ax.annotate(source_id, (xc[i] + r2, yc[i] + r2), fontsize=9, color='blue')

        # create an axes on the right side of ax. The width of cax will be 4%
        # of ax and the padding between cax and ax will be fixed at 0.05 inch.
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.06)

        # Define the colorbar and fix the labels
        cbar = plt.colorbar(im, cax=cax)

        # Define labels
        cbar.set_label(r'Arbitraty Flux Count Rate')

    if join_frames:

        photo_df['id'] = photo_df['id'].astype(int)
        s['id'] = s['id'].astype(int)

        photo_df = pd.merge(s, photo_df, on=['id', 'id'])

    return photo_df


def get_header_info(hdu_list):
    """

    Args:
        hdu_list:

    Returns:

    """

    # Getting header information
    header = hdu_list[0].header

    # Getting gain from header or from explicit input
    _gain = np.float(header['gain'])

    # Getting another header
    # _filename = header['filename']
    _dateobs = str(header['date-obs'])
    _object = str(header['object']).split()[0]

    _filter = str(header['filter'])
    _airmass = np.float(header['airmass'])
    _exptime = np.float(header['exptime'])

    return _dateobs, _object, _filter, _airmass, _exptime, _gain


def get_sigma_clipped_statistic(array, sigma_clip=2.0):
    """

    Args:
        array:
        sigma_clip:

    Returns:

    """

    mean, median, stddev = sigma_clipped_stats(array, sigma=sigma_clip)

    low_std, upp_std = mean - sigma_clip * stddev, mean + sigma_clip * stddev

    return mean, median, stddev, low_std, upp_std


def add_zpt_column(hdu_list, sources, gain=1.0, show=True, sigma_clip=2.0):
    """

    Args:
        hdu_list:
        sources:
        gain:
        show:
        sigma_clip:

    Returns:

    """

    # Passing sources to updated_sources
    updated_sources = sources

    # Initializing filter parameters
    ext_coef = {'J': 0.1, 'H': 0.05, 'Ks': 0.05}

    # Filter dictionary from 2MASS vs Four Star filter
    filter_dict = {'J': 'j_m', 'H': 'h_m', 'Ks': 'k_m'}

    if not isinstance(gain, (float, int, type(None))):
        raise TypeError('Please, the gain parameter must be a float, int or None type.')

    if not isinstance(show, (bool, type(None))):
        raise TypeError('Please, the show parameter must be boolean (True or False) or None type.')

    _dateobs, _object, _filter, _airmass, _exptime, __gain = get_header_info(hdu_list)

    # Getting gain from header or from explicit input
    _gain = np.float(gain) if gain is not None else np.float(__gain)

    # Remove data with NaN mag from 2MASS catalog dependin filter_dict[_filter]
    sources.dropna(subset=[filter_dict[_filter]], inplace=True)

    # Getting mag. of reference stars
    mag_ref = sources[filter_dict[_filter]]

    # Getting aperture photometry
    aper_sum = np.asarray(sources['aper_sum_bkgsub'])

    # Calculatig zero points
    zpt = mag_ref + (2.5 * np.log10(_gain * aper_sum / _exptime)) + (_airmass - 1) * ext_coef[_filter]

    # Update sources data frame with:
    # 1. distance from center
    updated_sources = get_distance_from_center(hdu_list, sources, center=None, scale=None)

    # 2. Zero Points
    updated_sources['zpt'] = zpt

    # Fazer um plot com linha sobre a média do zpt +/- 1sigma
    if show:
        display_zpt_results(hdu_list, sources, sigma_clip=sigma_clip)

    return updated_sources


def display_zpt_results(hdu_list, sources, sigma_clip=2.0):
    """
    Args:
        hdu_list:
        sources:
        sigma_clip:

    Returns:

    """
    # Getting data
    data = hdu_list[0].data

    # Getting filter dictionary
    filter_dict = {'J': 'j_m', 'H': 'h_m', 'Ks': 'k_m'}

    # Getting header information
    _dateobs, _object, _filter, _airmass, _exptime, _gain = get_header_info(hdu_list)

    # Getting mag, radio and zpt from data frame sources
    zpt = sources.zpt
    r_arcmin = sources.r_arcmin
    mag = sources[filter_dict[_filter]]

    # Getting zpt statistics
    mean, median, stddev, low_std, upp_std = get_sigma_clipped_statistic(zpt, sigma_clip=sigma_clip)

    if len(sources) < 1:
        print('No sources provided as input. The source input data frame is empty.')

    elif len(sources) > 1:

        fig = plt.figure(figsize=(10, 6))

        gs = gridspec.GridSpec(3, 5)
        ax0 = plt.subplot(gs[0:3, 0:3])
        ax1 = plt.subplot(gs[0, 3:5])
        ax2 = plt.subplot(gs[1, 3:5])
        ax3 = plt.subplot(gs[2, 3:5])

        # 0. Plotting data + sources
        # Getting image normalization
        norm_image = simple_norm(data, stretch='linear', percent=99.5)

        # Plotting image
        ax0.imshow(data, origin='lower', norm=norm_image, cmap='Greys',
                   interpolation='nearest')  # cmap='viridis', Greys

        # Getting x,y positions to display circles around sources
        xc, yc = list(sources.xcentroid), list(sources.ycentroid)

        # Radio of circles in pixels
        radio = 25.0

        for i, source_id in enumerate(list(sources['id'])):
            circle = Circle((xc[i], yc[i]), radio, facecolor='None', edgecolor='blue', linewidth=1.5, alpha=0.3)
            ax0.add_patch(circle)

        ax0.set_xlabel('X (pixels)')
        ax0.set_ylabel('Y (pixels)')
        ax0.set_title('{} - {} - {} - Airmass = {}'.format(_dateobs, _object, _filter, _airmass))

        # 1. Zpt Histogram
        mu, sigma = median, stddev  # mean an standard deviation of distribution

        # Number of bins rule >> R(n ^ (1 / 3)) / (3.49σ)
        # Ref. https://stackoverflow.com/questions/33458566/how-to-choose-bins-in-matplotlib-histogram
        R = sources.zpt.max() - sources.zpt.min()
        n = len(sources.zpt)
        num_bins = int(R * (n ** (1 / 3.)) / (3.49 * sigma))

        # Numpy has some options to get optimal number of bins: 'auto', 'fd', etc
        #num_bins = len(np.histogram_bin_edges(sources.zpt, bins='fd'))

        # Label
        lab = 'N = {:d}'.format(len(zpt))

        # 1.1
        # mu, sigma = norm.fit(zpt)
        n, bins, patches = ax1.hist(zpt, num_bins, density=True, label=lab, alpha=0.7, color='green')

        best_fit_line = norm.pdf(bins, mu, sigma)
        ax1.plot(bins, best_fit_line, '--', color='orange')

        # 1.2
        # n, bins, patches = ax1.hist(zpt, num_bins, density=1, alpha=0.5, label=lab)

        # add a 'best fit' line
        # y = ((1 / (np.sqrt(2 * np.pi) * sigma)) * np.exp(-0.5 * (1 / sigma * (bins - mu))**2))
        # ax1.plot(bins, y, '--')

        # 1.3 Using data frame
        # sources.hist(column='zpt', ax=ax1, bins=num_bins, grid=False, label=lab, alpha=0.5)

        ax1.set_xlabel('Zpt (mag.)')
        ax1.set_ylabel('N (density)')
        ax1.legend(loc='best')

        # 2. Zpt vs R (arcmin)
        # Plotting scatter
        scatter = ax3.scatter(r_arcmin, zpt, c=mag, cmap='viridis')

        # Getting 1d fitting parameters
        pars = np.polyfit(sources.r_arcmin, sources.zpt, 1)

        # Setting label
        lab = 'd(Zpt)/dR = {:3.2f}'.format(pars[0])
        
        # Fit 1D model + 95% confidence interval (without scatter)
        sns.regplot(x="r_arcmin", y="zpt", data=sources, ci=95, ax=ax3, scatter=False, label=lab)

        # Setting colorbar axis
        divider = make_axes_locatable(ax3)
        cax = divider.append_axes("right", size="4%", pad=0.04)
        cb = fig.colorbar(scatter, ax=ax3, cax=cax)
        cb.set_label('{} mag.'.format(_filter), fontsize=10)

        # Axis limit + x and y labels
        # ax3.axis([xmin, xmax, ymin, ymax])
        ax3.set_ylabel('Zpt (mag.)')
        ax3.set_xlabel('R (arcmin)')
        ax3.legend(loc='best', fontsize=10)

        # 3. Zpt vs Mag
        scatter = ax2.scatter(mag, zpt, c=r_arcmin, cmap='cividis')

        # Setting colorbar
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="4%", pad=0.04)
        cb = fig.colorbar(scatter, ax=ax2, cax=cax)
        cb.set_label('R (arcmin)', fontsize=10)

        # Setting label
        r'$\sin (x)$'
        lab = r'${:4.2f} \pm {:3.2f}$'.format(median, stddev)

        # Showing median and low/high sigma_clipping lines
        ax2.axhline(y=mean, color='r', linestyle='-', label=lab)
        ax2.axhline(y=low_std, color='g', linestyle='-', alpha=0.5)
        ax2.axhline(y=upp_std, color='g', linestyle='-', alpha=0.5)
        # ax3.fill_between(np.arange(0, 9, 1), low_std, upp_std, alpha=0.3, color='green')

        # Axis limit + x and y labels
        # ax2.axis([xmin, xmax, ymin, ymax])
        ax2.set_ylabel('Zpt (mag.)')
        ax2.set_xlabel('{} (mag.)'.format(_filter))
        ax2.legend(loc='best', fontsize=10)

        plt.tight_layout()


def run_zpt(filename, f_fwhm=5, f_threshold=15.0, f_round_limit=0.2, f_saturation=26641.4,
            c_max_separation=5, z_gain=1.0, z_sigma_clip=2.0 ):

    """

    Args:
        filename:
        seeing:
        scale:
        threshold:
        saturation:

    Returns:

    """
    hdu = get_hdulist(filename)

    header, data = hdu[0].header, hdu[0].data

    # Finding Sources in the data array
    _, sources, _ = find_stars(data, fwhm=f_fwhm, threshold=f_threshold, saturation=f_saturation,
                               round_low=-f_round_limit, round_high=f_round_limit, show=None)

    # Getting coordinates
    coordinates = get_coordinates(hdu, ra_deg=None, dec_deg=None, frame='icrs')

    # Querying catalog
    catalog = get_two_mass_sources(coordinates, spatial='Cone', radius=7.5 * u.arcmin)

    # Getting main dataframe with matching sources
    dataframe = match_sources(sources, catalog, hdu, max_separation=c_max_separation * u.arcsec)

    # Selection a few keywords from dataframe
    keywods = ['id', 'xcentroid', 'ycentroid', 'sharpness', 'roundness1', 'roundness2', 'npix', 'sky', 'peak', 'flux',
               'mag', 'ra', 'dec', 'clon', 'clat', 'designation', 'j_m', 'j_cmsig', 'j_snr', 'h_m', 'h_cmsig', 'h_snr',
               'k_m', 'k_cmsig', 'k_snr']

    #
    dataframe = dataframe[keywods]

    dataframe = get_photo(data, dataframe, aperture=5 * f_fwhm, annulus_width=10.0, show=None, join_frames=True)

    dataframe = add_zpt_column(hdu, dataframe, gain=z_gain, show=None, sigma_clip=z_sigma_clip)

    display_zpt_results(hdu, dataframe, sigma_clip=z_sigma_clip)

    return dataframe


