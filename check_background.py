#!/usr/bin/env python

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import sep
import aplpy
from astropy import wcs
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clipped_stats
from photutils.aperture import ApertureStats
from photutils.aperture import SkyCircularAnnulus

from hostphot._constants import font_family
from hostphot.utils import suppress_stdout

import warnings
from astropy.utils.exceptions import AstropyWarning

def extract_image(file):
    """Obtains the data and other information from a FITS file.

    Parameters
    ----------
    file: str
        Name of the FITS file.

    Returns
    -------
    data: ndarray
        Image data/counts.
    header: ~fits.header
        Image header.
    img_wcs: ~astropy.wcs
        Image WCS.
    hdu: ~fits.hdu
        Image Header Data Unit.
    """
    hdu = fits.open(file)
    
    header = hdu[0].header
    data = hdu[0].data
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        img_wcs = wcs.WCS(header, naxis=2)

    data = data.astype(np.float64)
    
    return data, header, img_wcs, hdu

def get_sep_stats(data, **sep_kwargs):
    """Obtains the background mean and rms from an array
    using SEP (SExtractor).

    Parameters
    ----------
    data: ndarray
        Image data/counts.

    Returns
    -------
    bkg_mean: float
        Background mean.
    bkg_rms: float
        Background root-mean-square.
    """
    bkg = sep.Background(data, **sep_kwargs)
    bkg_mean = bkg.globalback
    bkg_rms = bkg.globalrms
    
    return bkg_mean, bkg_rms

def get_astropy_stats(data, sigma=3.0):
    """Obtains the background mean, median and std from an array
    using Astropy's sigma-clipping stats.

    Parameters
    ----------
    data: ndarray
        Image data/counts.
    sigma: float, default '3.0'
        Sigma used for the sigma clipping.

    Returns
    -------
    mean: float
        Background mean.
    mean: float
        Background median.
    std: float
        Background standard deviation.
    """
    mean, median, std = sigma_clipped_stats(data, sigma=sigma)
    
    return mean, median, std

def get_target_stats(data, img_wcs, ra, dec, r_in, r_out):
    """Obtains the background mean, median and std around
    the given coordinates using an annulus.

    Parameters
    ----------
    data: ndarray
        Image data/counts.
    img_wcs: ~astropy.wcs
        Image WCS.
    ra: float
        Right ascension.
    dec: float
        Declination.
    r_in: float
        Inner radius of the annulus.
    r_out: float
        Outer radius of the annulus.

    Returns
    -------
    aperture: ~astropy.aperture
        Annulus aperture used for the background.
    mean: float
        Background mean.
    mean: float
        Background median.
    std: float
        Background standard deviation.
    """
    coords = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame="icrs")

    r_in = r_in * u.arcsec
    r_out = r_out * u.arcsec
    aperture = SkyCircularAnnulus(coords, r_in=r_in, r_out=r_out)
    aperstats = ApertureStats(data, aperture, wcs=img_wcs) 

    return aperture, aperstats.mean, aperstats.median, aperstats.std
            
def plot_target(
    hdu,
    ra=None,
    dec=None,
    aperture=None,
    size=1.0,
    info_dict=None,
):
    """Plots the objects extracted with :func:`sep.extract()``.

    Parameters
    ----------
    hdu: ~fits.hdu
        Image Header Data Unit.
    ra: float, default 'None'
       Right ascension of an object, in degrees. Used for plotting the position of the object.
    dec: float, default 'None'
    aperture: ~astropy.aperture, default 'None'
        Annulus aperture used for the background.
    size: float, default '1.0'
        Size of the image to be plotted, in arcminutes.
    info_dict: dict, default 'None'
        Dictionary with background statistics.
    """
    figure = plt.figure(figsize=(10, 10))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", AstropyWarning)
        fig = aplpy.FITSFigure(hdu, figure=figure)

    with suppress_stdout():
        fig.show_grayscale(stretch="arcsinh")

    if (ra is not None) and (dec is not None):
        # resize image to show only around the coordinates
        size_arcmin = size*u.arcmin
        size_degree = size_arcmin.to(u.degree)
        fig.recenter(ra, dec, size_degree.value)
        
        if aperture is not None:
            # plot annulus
            fig.show_circles(
                ra,
                dec,
                aperture.r_in.to(u.degree),
                linewidth=2,
                edgecolor="r",
                label='Annulus',
                layer='r_in'
            )
            fig.show_circles(
                ra,
                dec,
                aperture.r_out.to(u.degree),
                linewidth=2,
                edgecolor="r",
                layer='r_out'
            )

    sep_mean, sep_std = info_dict['sep']
    astro_mean, astro_median, astro_std = info_dict['astro']
    target_mean, target_median, target_std = info_dict['target']
    sep_diff = info_dict['sep_diff']
    astro_diff = info_dict['astro_diff']
    
    text = 'Background stats\n'
    text += f'SEP: mean={sep_mean:.2f}, std={sep_std:.2f}\n'
    text += f'ASTROPY: mean={astro_mean:.2f}, median={astro_median:.2f}, std={astro_std:.2f}\n'
    text += f'Annulus: mean={target_mean:.2f}, median={target_median:.2f}, std={target_std:.2f}\n'
    text += f'$\Delta$(SEP): {sep_diff:.2f}$\sigma$, $\Delta$(ASTROPY): {astro_diff:.2f}$\sigma$'
    
    fig.add_label(0.04, 0.11, text, relative=True, **{"family": font_family, 
                                                      "size": 18, 
                                                      "weight":"bold",
                                                      "horizontalalignment":"left",
                                                      "bbox":{"boxstyle":"round", "facecolor":"white", "alpha":0.7},
                                                      #"alpha":0.6,
                                                     })
    
    # ticks
    fig.tick_labels.set_font(**{"family": font_family, "size": 18})
    fig.tick_labels.set_xformat("dd.dd")
    fig.tick_labels.set_yformat("dd.dd")
    fig.ticks.set_length(6)

    fig.axis_labels.set_font(**{"family": font_family, "size": 18})
    fig.set_theme("publication")

    plt.show()
    
def check_background(file, ra, dec, r_in=3, r_out=6, method='mean', show_plot=True, size=1.0):
    """Calculates the difference, in sigmas, between an image global
    background and the background around the given coordinates.
    
    diff = np.abs(bkg_mean1 - bkg_mean2)/np.sqrt(bkg_std1**2 + bkg_std2**2)
    
    Parameters
    ----------
    file: str
        Name of the FITS file.
    ra: float
        Right ascension.
    dec: float
        Declination.
    r_in: float, default '3'
        Inner radius of the annulus (in arcsec).
    r_out: float, default '6'
        Outer radius of the annulus (in arcsec).
    method: str, default 'mean'
        Method used to estimate the difference in background.
        Either 'mean' or 'median'. SEP only uses 'mean'.
    show_plot: bool default 'True'
        Whether to show the image with the annulus used.
    size: float, default '1.0'
        Size of the image to be plotted, in arcminutes.
    """
    data, header, img_wcs, hdu = extract_image(file)
    
    # background statistics
    sep_mean, sep_std = get_sep_stats(data)
    astro_mean, astro_median, astro_std = get_astropy_stats(data)
    aperture, target_mean, target_median, target_std = get_target_stats(data, 
                                                                        img_wcs,
                                                                        ra, dec, 
                                                                        r_in, r_out)
    
    # calculate difference in background level, in units of sigmas
    assert method in ['mean', 'median'], "Not a valid method!"
    
    if method=='mean':
        sep_diff = np.abs(target_mean-sep_mean)/np.sqrt(sep_std**2 + target_std)
        astro_diff = np.abs(target_mean-astro_mean)/np.sqrt(astro_std**2 + target_std)
    elif method=='median':
        sep_diff = np.abs(target_median-sep_mean)/np.sqrt(sep_std**2 + target_std)
        astro_diff = np.abs(target_median-astro_median)/np.sqrt(astro_std**2 + target_std)
    
    print(f'SEP: {np.round(sep_diff, 2)} sigmas')
    print(f'ASTROPY: {np.round(astro_diff, 2)} sigmas')
    
    if show_plot is True:
        info_dict = {'sep':[sep_mean, sep_std],
                     'astro':[astro_mean, astro_median, astro_std],
                     'target':[target_mean, target_median, target_std],
                     'sep_diff':sep_diff,
                     'astro_diff':astro_diff,
                    }
        plot_target(hdu, ra, dec, aperture, size, info_dict)
     
    # save output into a file
    out_dict = {'file':[file],
                'ra':[ra],
                'dec':[dec],
                'r_in':[r_in],
                'r_out':[r_out],
                'method':[method],
                'sep_diff':[sep_diff],
                'astro_diff':[astro_diff],
                'sep_mean':[sep_mean],
                'sep_std':[sep_std],
                'astro_mean':[astro_mean],
                'astro_median':[astro_median],
                'astro_std':[astro_std],
                'target_mean':[target_mean],
                'target_median':[target_median],
                'target_std':[target_std]
               }
    
    df = pd.DataFrame(out_dict)
    outfile = os.path.basename(file).replace('.fits', '')
    outfile = 'bkg_' + outfile + '.csv'
    df.to_csv(outfile, index=False)
        
        
def main(args=None):
    description = f"Checks image background to identify the need of templates for image subtraction"
    usage = "check_background file ra dec [options]"
    
    if not args:
        args = sys.argv[1:] if sys.argv[1:] else ["--help"]
        
    parser = argparse.ArgumentParser(prog='check_background',
                                     usage=usage,
                                     description=description
                                     )
    parser.add_argument("file",
                        type=str,
                        help="Name of the FITS file."
                        )
    parser.add_argument("ra",
                        type=float,
                        help="Right ascension."
                        )
    parser.add_argument("dec",
                        type=float,
                        help="Declination."
                        )
    parser.add_argument("--r_in",
                        dest="r_in",
                        action="store",
                        default=3,
                        type=float,
                        help="Inner radius of the annulus."
                        )
    parser.add_argument("--r_out",
                        dest="r_out",
                        action="store",
                        default=6,
                        type=float,
                        help="Outer radius of the annulus."
                        )
    parser.add_argument("-m"
                        "--method",
                        dest="method",
                        action="store",
                        default="mean",
                        choices=["mean", "median"],
                        type=str,
                        help=("Method used to estimate the difference in background."
                              "Either 'mean' or 'median'. SEP only uses 'mean'.")
                        )
    parser.add_argument("-p"
                        "--show_plot",
                        dest="show_plot",
                        action="store",
                        default=True,
                        choices=[True, False],
                        type=bool,
                        help="Whether to show the image with the annulus used."
                        )
    parser.add_argument("-s"
                        "--size",
                        dest="size",
                        action="store",
                        default=1,
                        type=float,
                        help="Size of the image to be plotted, in arcminutes."
                        )
    
    args = parser.parse_args(args)
    check_background(args.file, args.ra, args.dec, args.r_in, args.r_out, 
                     args.method, args.show_plot, args.size)

if __name__ == "__main__":
    main(sys.argv[1:])