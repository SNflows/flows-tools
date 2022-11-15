from abc import ABC, abstractmethod
from typing import Optional
from .target import Target
from .plan import Plan
from .corner import Corner
from .utils import numeric

from astropy.coordinates import SkyCoord
import astropy.units as u
import regions


class Instrument(ABC):
    """Instrument class"""

    @abstractmethod
    def __init__(self, coords: SkyCoord, rotation: Optional[u.Quantity] = 0.0 * u.deg):
        pass

    @abstractmethod
    def point(self, target: Target, plan: Plan) -> list[
        regions.RectangleSkyRegion, Optional[regions.RectangleSkyRegion]]:
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

    @property
    @abstractmethod
    def region_names(self):
        pass

    @abstractmethod
    def get_corners(self):
        pass


class Hawki(Instrument):
    """Class for storing the hardcoded chip and detector information of a VLT HAWKI pointing"""
    chip1_center = 3.75/2 * u.arcmin + 15 * u.arcsecond  # This is the distance to Chip 1 center from Field center.
    field_hw = 7.5 * u.arcminute  # full field height width
    chip1_hw = 3.5 * u.arcminute  # chip1 field height width
    chip1_offset = 112.5 * u.arcsecond  # default pointing offset of Hawki.

    def __init__(self, coords: SkyCoord, rotation: Optional[u.Quantity] = 0.0 * u.deg):
        # Coordinates
        self.rotation = rotation
        self.coords = coords  # initial (target) coord
        self.center_coords = self.default_point()  # Hawki center coord
        self.chip1_center_coords = self.offset()  # CCD4 center coord

        # Regions
        self.hawki_region, self.chip1_region = self.get_regions()

    def point(self, target: Target, plan: Plan) -> list[
        regions.RectangleSkyRegion, Optional[regions.RectangleSkyRegion]]:
        """point telescope to rot=rotation in degrees, alpha and delta offset in arcseconds"""
        self.coords = target.coords  # Assume unchanged
        if plan.rotate:
            self.rotation = plan.rotation

        self.center_coords = self.default_point(plan.alpha, plan.delta)
        self.chip1_center_coords = self.point_chip1(plan.alpha, plan.delta)
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

    @property
    def region_names(self) -> tuple[str, str]:
        return 'HAWK-I', 'chip1'

    def get_corners(self) -> list[Corner, Corner]:
        chip1_corners = Corner(self.chip1_center_coords.ra, self.chip1_center_coords.dec, self.chip1_hw)
        field_corners = Corner(self.center_coords.ra, self.center_coords.dec, self.field_hw)
        return [field_corners, chip1_corners]

    @staticmethod
    def make_region(coords: SkyCoord, width: u.Quantity, height: u.Quantity, angle: u.Quantity):
        return regions.RectangleSkyRegion(coords, width=width, height=height, angle=angle)

    def default_point(self, alpha: u.Quantity = 0.0 * u.arcsec, delta: u.Quantity = 0.0 * u.arcsec) -> SkyCoord:
        return self.offset(-self.chip1_offset + alpha, -self.chip1_offset + delta)  # Hawki_center_coord

    def point_chip1(self, alpha: u.Quantity = 0.0 * u.arcsec, delta: u.Quantity = 0.0 * u.arcsec) -> SkyCoord:
        new_coords = self.chip1_center_coords.spherical_offsets_by(alpha + 15/2 * u.arcsecond, delta + 15/2 *
                                                                   u.arcsecond)
        """returns the CCD4 (chip1) ra and dec center given a pointing in ra dec and rotation"""
        return self.center_coords.directional_offset_by(
            self.center_coords.position_angle(new_coords) + self.rotation,
            self.center_coords.separation(new_coords))
        # return self.center_coords.spherical_offsets_by(self.chip1_center, self.chip1_center)

    @staticmethod
    def _get_pa_sep(coord1: SkyCoord, coord2: SkyCoord) -> tuple[u.Quantity, u.Quantity]:
        pa = coord1.position_angle(coord2)
        sep = coord1.separation(coord2)
        return pa, sep
