import warnings
from typing import Optional, cast

import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io.fits import PrimaryHDU
from astropy.table import Table
from astroquery.simbad import Simbad
from astroquery.skyview import SkyView

from .utils import StrEnum, numeric


def query_simbad(
    coo_centre: SkyCoord, radius: numeric = u.Quantity(24, u.arcmin)
) -> tuple[Optional[Table], Optional[SkyCoord]]:
    """
    Query SIMBAD using cone-search around the position using astroquery
    Parameters:
        coo_centre (:class:`astropy.coordinates.SkyCoord`): Coordinates of centre of search cone.
        radius (float, optional):
    Returns:
        list: Astropy Table with SIMBAD information.
    .. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
    .. codeauthor:: Emir Karamehmetoglu <emir.k@phys.au.dk>
    """

    s = Simbad()
    s.ROW_LIMIT = 0
    s.remove_votable_fields("coordinates")
    s.add_votable_fields("ra(d;A;ICRS;J2000)", "dec(d;D;ICRS;2000)", "pmra", "pmdec")
    s.add_votable_fields("otype")
    s.add_votable_fields("flux(B)", "flux(V)", "flux(R)", "flux(I)", "flux(J)", "flux(H)", "flux(K)")
    s.add_votable_fields("flux(u)", "flux(g)", "flux(r)", "flux(i)", "flux(z)")

    rad = Angle(radius).arcmin
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning)
        results = s.query_criteria(
            f"region(circle, icrs, {coo_centre.icrs.ra.deg:.10f} {coo_centre.icrs.dec.deg:+.10f}, {rad}m)",
            otypes="Star",
        )

    if not results:
        return None, None

    # Rename columns:
    results.rename_column("MAIN_ID", "main_id")
    results.rename_column("RA_d_A_ICRS_J2000", "ra")
    results.rename_column("DEC_d_D_ICRS_2000", "dec")
    results.rename_column("PMRA", "pmra")
    results.rename_column("PMDEC", "pmdec")
    results.rename_column("FLUX_B", "B_mag")
    results.rename_column("FLUX_V", "V_mag")
    results.rename_column("FLUX_R", "R_mag")
    results.rename_column("FLUX_I", "I_mag")
    results.rename_column("FLUX_J", "J_mag")
    results.rename_column("FLUX_H", "H_mag")
    results.rename_column("FLUX_K", "K_mag")
    results.rename_column("FLUX_u", "u_mag")
    results.rename_column("FLUX_g", "g_mag")
    results.rename_column("FLUX_r", "r_mag")
    results.rename_column("FLUX_i", "i_mag")
    results.rename_column("FLUX_z", "z_mag")
    results.rename_column("OTYPE", "otype")
    results.remove_column("SCRIPT_NUMBER_ID")
    results.sort(["V_mag", "B_mag", "H_mag"])

    # Filter out object types which shouldn'r really be in there anyway:
    indx = (results["otype"] == "Galaxy") | (results["otype"] == "LINER") | (results["otype"] == "SN")
    results = results[~indx]

    if len(results) == 0:
        return None, None

    # Build sky coordinates object:
    simbad = SkyCoord(
        ra=results["ra"],
        dec=results["dec"],
        pm_ra_cosdec=results["pmra"],
        pm_dec=results["pmdec"],
        frame="icrs",
        obstime="J2000",
    )

    return results, simbad


class SkyViewSurveys(StrEnum):
    TWO_MASS_H = "2MASS-H"
    DSS = "DSS"
    DSS_B = "DSS2-B"
    DSS_R = "DSS2-R"
    SDSS_G = "SDSS-G"
    SDSS_R = "SDSS-R"


def query_skyview(
    ra: float,
    dec: float,
    pixels: int = 2500,
    radius: numeric = 20,
    scale: str = "linear",
    survey: str = SkyViewSurveys.DSS,
) -> PrimaryHDU:
    """
    Query SkyView using astroquery
    Parameters:
        ra (float): Right Ascension of centre of search cone.
        dec (float): Declination of centre of search cone.
        radius (float, optional):
    Returns:
        dict: Dictionary of HDU objects.
    """
    qradius: u.Quantity = radius << u.arcmin  # type: ignore # Convert to astropy units
    out = SkyView.get_images(
        position="{}, {}".format(ra, dec),
        survey=survey,
        pixels=str(pixels),
        coordinates="J2000",
        scaling=scale,
        radius=qradius,
    )
    hdul = out.pop()
    return cast(PrimaryHDU, hdul.pop())


def query_2mass_image(
    ra: float,
    dec: float,
    pixels: int = 2500,
    radius: numeric = 50,
    scale: str = "linear",
    survey: str = SkyViewSurveys.TWO_MASS_H,
) -> PrimaryHDU:  # ImageHDU
    qradius: u.Quantity = radius << u.arcmin  # type: ignore # Convert to astropy units
    out = SkyView.get_images(
        position="{}, {}".format(ra, dec),
        survey=survey,
        pixels=str(pixels),
        coordinates="J2000",
        scaling=scale,
        radius=qradius,
    )
    hdul = out.pop()
    return cast(PrimaryHDU, hdul.pop())
