import warnings
from astropy.coordinates import SkyCoord, Angle
from astroquery.simbad import Simbad
import astropy.units as u
from astroquery.skyview import SkyView
from astropy.table import Table
from typing import Optional


def query_simbad(coo_centre, radius=24*u.arcmin) -> tuple[Optional[Table], Optional[SkyCoord]]:
    """
    Query SIMBAD using cone-search around the position using astroquery.
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
    s.remove_votable_fields('coordinates')
    s.add_votable_fields('ra(d;A;ICRS;J2000)', 'dec(d;D;ICRS;2000)', 'pmra', 'pmdec')
    s.add_votable_fields('otype')
    s.add_votable_fields('flux(B)', 'flux(V)', 'flux(R)', 'flux(I)', 'flux(J)', 'flux(H)', 'flux(K)')
    s.add_votable_fields('flux(u)', 'flux(g)', 'flux(r)', 'flux(i)', 'flux(z)')

    rad = Angle(radius).arcmin
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', category=UserWarning)
        results = s.query_criteria((f'region(circle, icrs, {coo_centre.icrs.ra.deg:.10f} '
                                    f'{coo_centre.icrs.dec.deg:+.10f}, {rad}m)'), otypes='Star')

    if not results:
        return None, None

    # Rename columns:
    results.rename_column('MAIN_ID', 'main_id')
    results.rename_column('RA_d_A_ICRS_J2000', 'ra')
    results.rename_column('DEC_d_D_ICRS_2000', 'dec')
    results.rename_column('PMRA', 'pmra')
    results.rename_column('PMDEC', 'pmdec')
    results.rename_column('FLUX_B', 'B_mag')
    results.rename_column('FLUX_V', 'V_mag')
    results.rename_column('FLUX_R', 'R_mag')
    results.rename_column('FLUX_I', 'I_mag')
    results.rename_column('FLUX_J', 'J_mag')
    results.rename_column('FLUX_H', 'H_mag')
    results.rename_column('FLUX_K', 'K_mag')
    results.rename_column('FLUX_u', 'u_mag')
    results.rename_column('FLUX_g', 'g_mag')
    results.rename_column('FLUX_r', 'r_mag')
    results.rename_column('FLUX_i', 'i_mag')
    results.rename_column('FLUX_z', 'z_mag')
    results.rename_column('OTYPE', 'otype')
    results.remove_column('SCRIPT_NUMBER_ID')
    results.sort(['V_mag', 'B_mag', 'H_mag'])

    # Filter out object types which shouldn'r really be in there anyway:
    indx = (results['otype'] == 'Galaxy') | (results['otype'] == 'LINER') | (results['otype'] == 'SN')
    results = results[~indx]

    if len(results) == 0:
        return None, None

    # Build sky coordinates object:
    simbad = SkyCoord(
        ra=results['ra'],
        dec=results['dec'],
        pm_ra_cosdec=results['pmra'],
        pm_dec=results['pmdec'],
        frame='icrs',
        obstime='J2000')

    return results, simbad


def query_2mass_image(ra, dec, pixels=2500, radius=50):
    out = SkyView.get_images(position='{}, {}'.format(ra, dec),
                             survey='2MASS-H', pixels=str(pixels),
                             coordinates='J2000', scaling='Linear', radius=radius * u.arcmin)
    return out[0][0]
