from __future__ import annotations
from abc import ABC, abstractmethod
from tendrils import api
from astropy.time import Time
from astropy.table import Table
from .catalogs import query_simbad, query_2mass_image
from astropy.coordinates import SkyCoord, ICRS
import astropy.units as u
from erfa import ErfaWarning
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astropy.wcs.utils import celestial_frame_to_wcs
import regions
from .plots import plot_image
from astropy.visualization import ZScaleInterval
from dataclasses import dataclass, field
from typing import Union, Optional
import warnings

numeric = Union[float, int, u.Quantity]

# Most useless warnings ever spammed for every operation by this package!
warnings.filterwarnings('ignore', category=ErfaWarning, append=True)


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
    parser.add_argument('-p', '--plot', help='whether to query images and plot',
                        type=bool, default=False, action='store')

    args = parser.parse_args()
    if args.target == 'None':
        parser.error('target id or name not provided, use -t <targetid> or <targetname>')
    elif args.target.isnumeric():
        args.target = int(args.target)
    return args.rotate, args.target, args.shifta, args.shiftd


def process_corner(corners):
    '''Given corners of a rectangle defined as astropy Quantity objects, return it as an np array of floats'''
    _points = [u.quantity.Quantity(corner) for corner in corners]
    return np.array(_points)


class Instrument(ABC):
    """Instrument class"""

    @abstractmethod
    def __init__(self, coords: SkyCoord, rotation: Optional[u.Quantity] = 0.0 * u.deg):
        pass

    @abstractmethod
    def point(self, target: Target, plan: Plan) -> list[regions.RectangleSkyRegion, Optional[regions.RectangleSkyRegion]]:
        pass

    @abstractmethod
    def default_point(self):
        pass

    @abstractmethod
    def offset(self):
        pass

    @abstractmethod
    def get_regions(self):
        pass

    @abstractmethod
    def get_corners(self):
        pass


class Corner:

    def __init__(self, x: u.Quantity, y: u.Quantity, hw: u.Quantity):
        self.x = x
        self.y = y
        self.hw = hw
        self.corner_xy = self.set_corners()
        self.corners = self.process_corner()

    def set_corners(self):
        """Get corners of a rectangle for a given ra dec and side length
        :return: list[corner 1, corner 2, ..]
        """
        return [(self.x - self.hw, self.y - self.hw), (self.x - self.hw, self.y + self.hw),
                (self.x + self.hw, self.y + self.hw), (self.x + self.hw, self.y - self.hw)]

    def process_corner(self):
        '''Given corners of a rectangle defined as astropy Quantity objects, return it as an np array of floats'''
        _points = [u.quantity.Quantity(corner) for corner in self.corner_xy]
        return np.array(_points)


class Hawki(Instrument):
    """Class for storing the hardcoded chip and detector information of a VLT HAWKI pointing"""
    chip1_center = 3.75 * u.arcmin + 15 * u.arcsecond  # This is the distance to Chip 1 center from Field center.
    field_hw = 7.5 * u.arcminute  # full field height width
    chip1_hw = 3.5 * u.arcminute  # chip1 field height width
    chip1_offset = 112.5 * u.arcsecond  # default pointing offset of Hawki.

    def __init__(self, coords: SkyCoord, rotation: Optional[u.Quantity] = 0.0 * u.deg):
        # Coordinates
        self.rotation = rotation
        self.coords = coords  # initial (target) coord
        self.center_coords = self.default_point()  # Hawki center coord
        self.chip1_center_coords = self.point_chip1()  # CCD4 center coord

        # Regions
        self.hawki_region, self.chip1_region = self.get_regions()

    def point(self, target: Target, plan: Plan) -> list[
        regions.RectangleSkyRegion, Optional[regions.RectangleSkyRegion]]:
        """point telescope to rot=rotation in degrees, alpha and delta offset in arcseconds"""
        self.coords = target.coords  # Assume unchanged
        if plan.shift:
            self.center_coords = self.default_point(plan.alpha, plan.delta)
            self.chip1_center_coords = self.point_chip1()
        if plan.rotate:
            self.rotation = plan.rotation
        self.hawki_region, self.chip1_region = self.get_regions()
        return [self.hawki_region, self.chip1_region]

    def offset(self, shifta: Optional[numeric] = 0.0, shiftd: Optional[numeric] = 0.0) -> SkyCoord:
        shifta = shifta << u.arcsecond
        shiftd = shiftd << u.arcsecond
        return self.coords.spherical_offsets_by(shifta, shiftd)

    def get_regions(self) -> tuple[regions.RectangleSkyRegion, regions.RectangleSkyRegion]:
        hawki_region = self.make_region(self.center_coords, self.field_hw, self.field_hw, self.rotation)
        chip1_region = self.make_region(self.chip1_center_coords, self.chip1_hw, self.chip1_hw, self.rotation)
        return hawki_region, chip1_region

    def get_corners(self) -> list[Corner, Corner]:
        chip1_corners = Corner(self.chip1_center_coords.ra, self.chip1_center_coords.dec, self.chip1_hw)
        field_corners = Corner(self.center_coords.ra, self.center_coords.dec, self.field_hw)
        return [field_corners, chip1_corners]

    @staticmethod
    def make_region(coords: SkyCoord, width: u.Quantity, height: u.Quantity, angle: u.Quantity):
        return regions.RectangleSkyRegion(coords, width=width, height=height, angle=angle)

    def default_point(self, alpha: u.Quantity = 0.0 * u.arcsec, delta: u.Quantity = 0.0 * u.arcsec) -> SkyCoord:
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
    coords: Optional[SkyCoord] = None
    info: dict[str, ...] = field(default_factory=dict)

    def __post_init__(self):
        if self.coords is None:
            self.coords = SkyCoord(self.ra << u.deg, self.dec << u.deg, frame='icrs')


@dataclass
class Plan:
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

    def __init__(self, instrument: Instrument, target: Target, plan: Plan, verbose: bool = False):
        self.target = target  # target info
        if not isinstance(self.target.coords, SkyCoord):
            raise TypeError('Target.coords must be an instance of SkyCoord representing target coordinates.')
        self.plan = plan  # Store offsets and rotations.
        self.ins = instrument  # store instrument specific info
        self.regions = self.ins.point(target, plan)
        self.nregions = len(self.regions)  # How many sub-regions of important passed by Instrument.?

        self.refcat_coords, self.refcat = self._make_refcat_catalog(mask=True, mask_dict={'H_mag': 14.0})
        self.simbad_coords, self.simbad = self._make_simbad_catalog()

        if not verbose:
            print((f'Simbad and refcat skyframes are '
                   f'equivalent: {self.refcat_coords.frame.is_equivalent_frame(self.simbad_coords.frame)}'
                   f'. If not, there might be a slight mismatch in alignment in current year.'))
        self.wcs = self.get_wcs(self.refcat_coords.frame)

    def _make_refcat_catalog(self, mask: bool = True, mask_dict: Optional[dict[str, float]] = None) \
            -> tuple[SkyCoord, Table]:
        refcat = api.get_catalog(self.target.tid)
        if mask and mask_dict is not None:
            masks = []
            keys_values = mask_dict.items()
            for key, value in keys_values:
                masks.append(refcat['references'][str(key)] <= float(value))
            if len(masks) == 1:
                refcat['references'] = refcat['references'][masks[0]]
            else:
                for msk in masks:
                    refcat['references'] = refcat['references'][msk]
        refcat_coords = SkyCoord(ra=refcat['references']['ra'],
                                 dec=refcat['references']['decl'],
                                 pm_ra_cosdec=refcat['references']['pm_ra'],
                                 pm_dec=refcat['references']['pm_dec'],
                                 obstime=Time(2015.5, format='decimalyear'))
        return refcat_coords, refcat['references']

    def _make_simbad_catalog(self) -> tuple[SkyCoord, Table]:
        simbad = query_simbad(self.target.coords)
        # propagate Simbad catalog coords to refcat reference year
        simbad_coords = simbad[1].apply_space_motion(new_obstime=Time(2015.5, format='decimalyear'))
        return simbad_coords, simbad[0]


    def get_wcs(self, frame: Optional[object] = None, header: Optional[object] = None) -> WCS:
        if frame is not None:
            wcs = celestial_frame_to_wcs(frame)
            wcs.wcs.crval = np.array((self.target.ra, self.target.dec))
            return wcs
        elif header is not None:
            return WCS(header)
        else:
            raise AttributeError("Both frame and header were None, cannot calculate WCS")

    def get_image(self, pixels=2500, radius=50):
        return query_2mass_image(self.tar.ra, self.tar.dec, pixels, radius)

    def regions_to_physical(self) -> list[regions.RectangleSkyRegion]:
        pixel_regions = []
        for region in self.regions:
            pixel_regions.append(region.to_pixel(self.wcs))
        return pixel_regions

    def check_bright_stars(self, region: regions.RectangleSkyRegion, wcs: Optional[WCS] = None) -> np.ndarray:
        wcs = self.wcs if wcs is None else wcs
        # Check bright stars
        simbad_stars = region.contains(self.simbad_coords, wcs)
        catalog_stars = region.contains(self.refcat_coords, wcs)

        bright_stars = np.hstack((self.simbad[simbad_stars]['H_mag'].data, self.refcat[catalog_stars]['H_mag'].data))
        if np.ma.isMaskedArray(bright_stars):
            bright_stars = bright_stars.data[~np.isnan(bright_stars.data)]

        print('Brightest star has H-mag = {0:.1f}'.format(np.round(bright_stars.min(), 1)))
        return bright_stars


class Plotter:
    """
    Takes an observer with WCS, Target, Plan, Instrument, Regions, and Corners and makes a finding chart.
    """

    def __init__(self, obs: Observer):
        self.obs = obs

    def plot(self):
        """not implemented yet"""
        pass


def main():
    # Parse input
    rot, tid, shifta, shiftd = parse()

    # Get query flows database
    target_info = api.get_target(tid)
    target_info['skycoord'] = SkyCoord(target_info['ra'] * u.deg, target_info['decl'] * u.deg)

    # Create Observer
    target = Target(tid, target_info['ra'], target_info['decl'], target_info['skycoord'], target_info)
    plan = Plan(rot, shifta, shiftd)
    hawki = Hawki(target.coords)
    obs = Observer(hawki, target, plan)

    # Print brightest star in field and get list of bright stars.
    brightest_stars = obs.check_bright_stars(region=obs.regions[0])

if __name__ == '__main__':
    main()
