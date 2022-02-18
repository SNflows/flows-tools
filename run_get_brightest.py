from __future__ import annotations

from tendrils import api
from astropy.time import Time
from astropy.table import Table
from catalogs import query_simbad, query_2mass_image
from astropy.coordinates import SkyCoord, ICRS
import astropy.units as u
import argparse
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from astropy.wcs import WCS
from astropy.wcs.utils import celestial_frame_to_wcs
import regions
from plots import plot_image
from astropy.visualization import ZScaleInterval
from dataclasses import dataclass, field
from typing import Union, Optional
numeric = Union[float, int, u.Quantity]
from reproject import reproject_interp
from typing import Union

def parse():
    '''parse command line input to get target, position angle (rotate), alpha and delta offsets (shifta, shiftd)'''
    parser = argparse.ArgumentParser(description='Calculate Brightest Star')
    parser.add_argument('-t', '--target', help="calculate for this targetname or targetid",
                        type=str, default='None', action='store')
    parser.add_argument('-r', '--rotate', help='rotation angle in degrees',
                        type=float, default=0.0, action='store')
    parser.add_argument('-a', '--shifta', help='shift alpha in arcsec',
                        type=float, default=0.0, action='store')
    parser.add_argument('-d', '--shiftd', help='shift delta in arcsec',
                        type=float, default=0.0, action='store')
    args = parser.parse_args()
    if args.target == 'None': parser.error('target id or name not provided, use -t <targetid> or <targetname>')
    elif args.target.isnumeric(): args.target = int(args.target)
    return args.rotate, args.target, args.shifta, args.shiftd


def get_corners(ra, dec, radius):
    '''Get corners of a rectangle for a given ra dec and side length'''
    return [(ra - radius, dec - radius), (ra - radius, dec + radius),
            (ra + radius, dec + radius), (ra + radius, dec - radius)]


def process_corner(corners):
    '''Given corners of a rectangle defined as astropy Quantity objects, return it as an np array of floats'''
    _points = [u.quantity.Quantity(corner) for corner in corners]
    return np.array(_points)

def get_flows_info(tid):
    """This function is used for getting information from the flows database and catalog via the API."""

    # Get refcat2 catalog
    c = api.get_catalog(tid)
    ref = c['references']
    refcoords = SkyCoord(ref['ra'], ref['decl'])

    # Get target info
    target_info = api.get_target(tid)
    tar = c['target']
    ra_tar = tar['ra'].data[0]
    dec_tar = tar['decl'].data[0]

    # Get simbad catalog
    target_info['skycoord'] = SkyCoord(target_info['ra'] * u.deg, target_info['decl'] * u.deg)
    simbad = query_simbad(target_info['skycoord'])

    # propagte Simbad catalog coords to 2mass reference year
    simbad_coords = simbad[1].apply_space_motion(new_obstime=Time(2000, format='decimalyear'))

    return ref, refcoords, tar, target_info, simbad, simbad_coords, ra_tar, dec_tar

class Instrument:
    """Instrument class"""
    # Rotation, tid, alpha, delta,
    #rotation = 0.0
    #rotate = False
    #alpha = 0.0
    #delta = 0.0
    #skip_shift = True

    def __init__(self):
        #self.rotation = rotation
        #self.alpha = alpha
        #self.delta = delta


    def point(self, rotation: float | u.Quantity, alpha: float | u.Quantity, delta: float | u.Quantity):
        """point telescope to rot=rotation in degrees, alpha and delta offset in arcseconds"""
        self.rotation = rotation
        self.alpha = alpha
        self.delta = delta

    def default_point(self):
        self.point()

    def offset(self):
        self.alpha *u.arcsecond + default_offset




    #def get_wcs(self,image):
    #		return WCS(image.header)

class Hawki(Instrument):
    """Class for storing the hardcoded chip and detector information of a VLT HAWKI pointing"""
    chip1_center = 3.75 * u.arcmin + 15 * u.arcsecond  # This is the distance to Chip 1 center from Field center.
    field_hw = 7.5 * u.arcminute  # full field height width
    chip1_hw = 3.5 * u.arcminute  # chip1 field height width
    chip1_offset = 112.5 * u.arcsecond # default pointing offset of Hawki.


    def __init__(self, coords: SkyCoord, rotation: Optional[u.Quantity] = 0.0 *u.deg):
        super().__init__()  # is this needed?

        # Coordinates
        self.rotation = rotation
        self.coords = coords # initial (target) coord
        self.center_coords = self.default_point()  # Hawki center coord
        self.chip1_center_coords = self.point_chip1()  # CCD4 center coord

        # Regions
        self.hawki_region, self.chip1_region = self.get_regions()

        # Position angle and separation between from new center to target (init) coords.
        # self.position_angle, self.separation = self._get_pa_sep(self.center_coords, self.coords)

    def point(self, plan: Plan) -> None:
        """point telescope to rot=rotation in degrees, alpha and delta offset in arcseconds"""
        self.coords = plan.target.coords # Assume unchanged
        if plan.shift:
            self.center_coords = self.default_point(plan.alpha, plan.delta)
            self.chip1_center_coords = self.point_chip1()
        if plan.rotate:
            self.rotation = plan.rotation
        self.hawki_region, self.chip1_region = self.get_regions()


    def offset(self, shifta: Optional[numeric] = 0.0, shiftd: Optional[numeric] = 0.0) -> SkyCoord:
        shifta = shifta << u.arcsecond
        shiftd = shiftd << u.arcsecond
        return self.coords.spherical_offsets_by(shifta, shiftd)

    def get_regions(self) -> tuple[regions.RectangleSkyRegion, regions.RectangleSkyRegion]:
        hawki_region = self.make_region(self.center_coords, self.field_hw, self.field_hw, self.rotation)
        chip1_region = self.make_region(self.chip1_center_coords, self.chip1_hw, self.chip1_hw, self.rotation)
        return hawki_region, chip1_region

    @staticmethod
    def make_region(coords:SkyCoord, width:u.Quantity, height:u.Quantity, angle:u.Quantity):
        return regions.RectangleSkyRegion(coords, width=width, height=height, angle=angle)

    def default_point(self, alpha: u.Quantity = 0.0*u.arcsec, delta: u.Quantity = 0.0*u.arcsec ) -> SkyCoord:
        return self.offset(-self.chip1_offset + alpha, -self.chip1_offset + delta)  # Hawki_center_coord

    def point_chip1(self) -> SkyCoord:
        """returns the CCD4 (chip1) ra and dec center given a pointing in ra dec and rotation"""
        return self.center_coords.spherical_offsets_by(self.chip1_center, self.chip1_center)

    @staticmethod
    def _get_pa_sep(coord1: SkyCoord, coord2: SkyCoord) -> tuple[u.Quantity, u.Quantity]:
        pa = coord1.position_angle(coord2)
        sep = coord1.separation(coord2)
        return pa, sep




@dataclass
class Target:
    tid: int | str  # target id
    ra: numeric
    dec: numeric
    info: dict[str, ...] = field(default_factory=dict)
    coords: Optional[SkyCoord] = None

    def __post_init__(self):
        if self.coords is None:
            self.coords = SkyCoord(self.ra<<u.deg, self.dec<<u.deg, frame='icrs')


@dataclass()
class Plan:
    target: Target
    rotation: numeric = 0.0
    alpha: numeric = 0.0
    delta: numeric = 0.0
    rotate: bool = False
    shift: bool = False

    def __post_init__(self):
        self.shift = self.set_shift()
        self.rotate = self.set_rotation()

    def set_rotation(self) -> bool:
        if self.rotation == 0.0:
            return False
        elif self.rotation != 0.0:
            return True

    def set_shift(self) -> bool:
        shift = (np.array((self.alpha, self.delta)) == 0.0).all()  # skip shift if alpha and delta 0
        return not shift


class Observer:

    def __init__(self, instrument: Hawki, target: Target, plan: Plan):
        self.tar = target  # target info
        if not isinstance(self.tar.coords, SkyCoord):
            raise TypeError('Target.coords must be an instance of SkyCoord representing target coordinates.')
        self.plan = plan  # Store offsets and rotations.
        self.ins = instrument(target.coords)  # store instrument specific info
        self.ins.point(plan)

        self.refcat_coords = self._make_refcat_catalog()
        self.simbad_coords = self._make_simbad_catalog()
        print((f'Simbad and refcat skyframes are '
               f'equivalent: {self.refcat_coords.frame.is_equivalent_frame(self.simbad_coords.frame)}'
               f'. If not, there might be a slight mismatch in alignment in current year.'))
        self.wcs = self.get_wcs(self.refcat_coords.frame)


    def _make_refcat_catalog(self) -> SkyCoord:
        refcat = api.get_catalog(self.target.tid)
        refcat_coords = SkyCoord(ra=refcat['references']['ra'],
                                 dec=refcat['references']['decl'],
                                 pm_ra_cosdec=refcat['references']['pm_ra'],
                                 pm_dec=refcat['references']['pm_dec'],
                                 obstime=Time(2015.5, format='decimalyear'))
        return refcat_coords

    def _make_simbad_catalog(self) -> SkyCoord:
        simbad = query_simbad(self.tar.coords)
        # propagate Simbad catalog coords to refcat reference year
        simbad_coords = simbad[1].apply_space_motion(new_obstime=Time(2015.5, format='decimalyear'))
        return simbad_coords

    @staticmethod
    def get_wcs(self, frame: Optional[object] = None, header: Optional[object] = None) -> WCS:
        if frame is not None:
            return celestial_frame_to_wcs(frame)
        elif header is not None:
            return WCS(header)
        else:
            raise AttributeError("Both frame and header were None, cannot calculate WCS")

    def get_image(self, pixels=2500, radius=50):
        return query_2mass_image(self.tar.ra, self.tar.dec, pixels, radius)


if __name__ == '__main__':
    # Parse input
    rot, tid, shifta, shiftd = parse()
    target_info = api.get_target(tid)
    # define useful values and get target info and catalog.
    #today = Time(datetime.datetime.today())

    # Create Observer
    target = Target(tid, target_info['ra'], target_info['decl'], target_info['coords'], target_info)
    plan = Plan(target, rot, shifta, shiftd)
    obs = Observer(Hawki, target, plan)

    obs

    #ref, refcoords, tar, target_info, simbad, simbad_coords, ra_tar, dec_tar = get_flows_info(tid)

    #hawki = Hawki(ra_tar, dec_tar)
    ra, dec = hawki.default_point_chip1() # Get default pointing roughly centered at chip1 center.
    chip1_center, field_hw, chip1_hw = hawki.chip1_center, hawki.field_hw, hawki.chip1_hw

    # Get corners of Chip 1 and full-field centered, where full field is centered at ra,dec.
    #corners = get_corners(ra, dec, field_hw)
    #chip1 = get_corners(ra - chip1_center,
    #					dec + chip1_center, chip1_hw)

    # Process to be easily plottable
    #chip1 = process_corner(chip1)
    #corners = process_corner(corners)

    wcs_H = hawki.get_image_wcs()
    ref_pix = wcs_H.all_world2pix(ref['ra'], ref['decl'], 0)

        # wcs_H = WCS(image.header)


    # Rotation of image.
    # This is disabled for now, we rotate the region boxes instead!
    #wcs_H2 = wcs_H.deepcopy()
    #if rot:
    #	rot_angle = rot
    #else:
    #	rot_angle = 0.0
    #theta = np.radians(0.0)
    #rot_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    #wcs_H2.wcs.pc = rot_matrix
    #array, footprint = reproject_interp((image.data, wcs_H), wcs_H2, shape_out=image.data.shape)
    #wcs_H = wcs_H2

    # Regions
    offset_hawki = 112.5 * u.arcsecond
    #if skip_shift:
    #	offset_alpha, offset_delta = offset_hawki, offset_hawki
    #	new_center_coord = target_info['skycoord']
    #	new_center_coords = self.tar_coords
    #else:
    #if not hawki.skip_shift and hawki.rotate:
    hawki.point(rot,shifta,shifta)

    offset_alpha, offset_delta = offset_hawki + shifta * u.arcsecond, offset_hawki + shiftd * u.arcsecond
    hawki.point(rot, shifta, shiftd)
    new_center_coord = target_info['skycoord'].spherical_offsets_by(-shifta * u.arcsecond, -shiftd * u.arcsecond)

    chip1_center_true = 3.75 / 2 * u.arcminute + 15 * u.arcsecond  # This is the distance to Chip 1 center from Field center.
    field_hw = 7.5 * u.arcminute  # full field height width
    chip1_hw = 3.5 * u.arcminute  # chip1 field height width

    Hawki_center_coord = target_info['skycoord'].spherical_offsets_by(-offset_alpha, -offset_delta)

    sep = Hawki_center_coord.spherical_offsets_to(target_info['skycoord'])  # unused
    field_angle = rot_angle * u.deg

    chip1_center_coord = Hawki_center_coord.directional_offset_by(
        Hawki_center_coord.position_angle(new_center_coord) + field_angle,
        Hawki_center_coord.separation(new_center_coord) + 15 * u.arcsecond
    )

    chip1_r = regions.RectangleSkyRegion(chip1_center_coord, width=chip1_hw, height=chip1_hw, angle=field_angle)
    Hawki_r = regions.RectangleSkyRegion(Hawki_center_coord, width=field_hw, height=field_hw, angle=field_angle)

    Hawki_pixel_region = Hawki_r.to_pixel(wcs_H)
    chip1_pixel_region = chip1_r.to_pixel(wcs_H)

    mask = ref['H_mag'] <= 14.0
    ref_pix_masked = wcs_H.all_world2pix(ref['ra'][mask], ref['decl'][mask], 0)

    corners_pix = wcs_H.all_world2pix(corners[:, 0] * u.deg, corners[:, 1] * u.deg, 0)
    tar_pix = wcs_H.all_world2pix(ra_tar, dec_tar, 0)
    chip1_pix = wcs_H.all_world2pix(chip1[:, 0], chip1[:, 1], 0)
    simbad_stars_pix = wcs_H.all_world2pix(simbad_coords.ra, simbad_coords.dec, 0)

    # Check bright stars
    simbad_stars = Hawki_r.contains(simbad_coords, wcs_H)
    catalog_stars = Hawki_r.contains(refcoords, wcs_H)

    bright_stars = np.hstack((simbad[0][simbad_stars]['H_mag'].data, ref[catalog_stars]['H_mag'].data))
    if np.ma.isMaskedArray(bright_stars): bright_stars = bright_stars.data[~np.isnan(bright_stars.data)]

    print('Brightest star has H-mag = {0:.1f}'.format(np.round(bright_stars.min(), 1)))

    # Plot
    # Set style using seaborn for prettier plots
    sns.set()
    zscale = ZScaleInterval()
    vmin, vmax = zscale.get_limits(image.data.flat)
    sns.set_theme(style='ticks', font_scale=2.5)
    fig, ax = plt.subplots(figsize=(20, 20), subplot_kw={'projection': wcs_H}, dpi=200)

    Hawki_pixel_region.plot(ax=ax, edgecolor='cyan', linestyle='-.', label='HAWK-I')
    chip1_pixel_region.plot(ax=ax, edgecolor='black', linestyle='-.', label='chip1')
    plot_image(array, ax=ax, cmap='viridis', scale='linear', vmin=vmin, vmax=vmax)

    ref_pix_masked = wcs_H.all_world2pix(ref['ra'][mask], ref['decl'][mask], 0)
    ax.scatter(ref_pix_masked[0], ref_pix_masked[1],
               facecolors='none', edgecolors='red', zorder=5, alpha=3.0, marker='o', s=200, label='bright stars')

    corners_pix = wcs_H.all_world2pix(corners[:, 0] * u.deg, corners[:, 1] * u.deg, 0)
    tar_pix = wcs_H.all_world2pix(ra_tar, dec_tar, 0)
    ax.scatter(tar_pix[0], tar_pix[1], marker='*', s=450, label='SN', color='orange')

    chip1_pix = wcs_H.all_world2pix(chip1[:, 0], chip1[:, 1], 0)
    simbad_stars = wcs_H.all_world2pix(simbad_coords.ra, simbad_coords.dec, 0)
    ax.scatter(simbad_stars[0], simbad_stars[1],
               facecolors='none', edgecolors='orange', zorder=5, alpha=3.0, marker='s', s=200, label='simbad stars')

    ax.legend(fontsize=18)
    fig.savefig('{}_HAWKI_H_FC.png'.format(target_info['target_name']))
