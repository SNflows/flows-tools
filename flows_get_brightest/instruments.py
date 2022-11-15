from abc import ABC, abstractmethod
from typing import Optional
from .target import Target
from .plan import Plan
from .corner import Corner
from .utils import numeric, is_quantity

from astropy.coordinates import SkyCoord
import astropy.units as u
import regions


class Instrument(ABC):
    """Instrument class"""

    @abstractmethod
    def __init__(self, coords: SkyCoord, rotation: Optional[u.Quantity] = u.Quantity(0.0, u.deg)) -> None:
        pass

    @abstractmethod
    def point(self, target: Target, plan: Plan) -> list[regions.RectangleSkyRegion]:
        pass

    @abstractmethod
    def default_point(self) -> SkyCoord:
        pass

    @abstractmethod
    def offset(self) -> SkyCoord:
        pass

    @abstractmethod
    def get_regions(self) -> list[regions.RectangleSkyRegion]:
        pass

    @property
    @abstractmethod
    def region_names(self) -> list[str]:
        pass
    
    @property
    @abstractmethod
    def nregions(self) -> int:
        pass

    @abstractmethod
    def get_corners(self) -> list[Corner]:
        pass


class Hawki(Instrument):
    """Class for storing the hardcoded chip and detector information of a VLT HAWKI pointing"""
    # This is the distance to Chip 1 center from Field center:
    chip1_center= u.Quantity(2.125, u.arcmin) #3.75 / 2, u.arcmin) + u.Quantity(15, u.arcsecond) 
    field_hw = u.Quantity(7.5, u.arcminute)  # full field height width
    chip1_hw = u.Quantity(3.5, u.arcminute)  # chip1 field height width
    chip1_offset = u.Quantity(112.5, u.arcsecond)  # default pointing offset of Hawki.

    def __init__(self, coords: SkyCoord, rotation: numeric = u.Quantity(0.0,u.deg)):
        # Coordinates
        self.rotation: u.Quantity = rotation << u.deg  # type: ignore
        self.coords = coords  # initial (target) coord
        self.center_coords = self.default_point()  # Hawki center coord
        self.chip1_center_coords = self.offset()  # CCD4 center coord

        # Regions
        self.hawki_region, self.chip1_region = self.get_regions()

    def point(self, target: Target, plan: Plan) -> list[regions.RectangleSkyRegion]:
        """point telescope to rot=rotation in degrees, alpha and delta offset in arcseconds"""
        self.coords = target.coords  # Assume unchanged
        if plan.rotate:
            self.rotation = u.Quantity(plan.rotation, u.deg)

        self.center_coords = self.default_point(plan.alpha, plan.delta)
        self.chip1_center_coords = self.point_chip1(plan.alpha, plan.delta)
        self.hawki_region, self.chip1_region = self.get_regions()
        return [self.hawki_region, self.chip1_region]

    def offset(self, shifta: Optional[numeric] = 0.0, shiftd: Optional[numeric] = 0.0) -> SkyCoord:
        shifta = shifta << u.arcsecond  # type: ignore
        shiftd = shiftd << u.arcsecond  # type: ignore
        return self.coords.spherical_offsets_by(shifta, shiftd)

    def get_regions(self) -> list[regions.RectangleSkyRegion]:
        hawki_region = self.make_region(self.center_coords, self.field_hw, self.field_hw, self.rotation)
        chip1_region = self.make_region(self.chip1_center_coords, self.chip1_hw, self.chip1_hw, self.rotation)
        return [hawki_region, chip1_region]

    @property
    def region_names(self) -> list[str]:
        return ['HAWK-I', 'chip1']
    
    @property
    def nregions(self) -> int:
        return 2

    def get_corners(self) -> list[Corner]:
        chip1ra, chip1dec = self.chip1_center_coords.ra, self.chip1_center_coords.dec
        fieldra, fielddec = self.center_coords.ra, self.center_coords.dec
        quants = [chip1ra, chip1dec, fieldra, fielddec]
        if not is_quantity(quants):
            raise ValueError(f'chip1_center_coords and center_coords must be SkyCoord with Quantity ra and dec.')
        chip1_corners = Corner(quants[0], quants[1], self.chip1_hw)
        field_corners = Corner(quants[2], quants[3], self.field_hw)
        return [field_corners, chip1_corners]

    @staticmethod
    def make_region(coords: SkyCoord, width: u.Quantity, height: u.Quantity, angle: u.Quantity):
        return regions.RectangleSkyRegion(coords, width=width, height=height, angle=angle)

    def default_point(self, alpha: u.Quantity = u.Quantity(0.0, u.arcsec), delta: u.Quantity = u.Quantity(0.0, u.arcsec)) -> SkyCoord:
        return self.offset(-self.chip1_offset + alpha, -self.chip1_offset + delta) # type: ignore # Hawki_center_coord  

    def point_chip1(self, alpha: u.Quantity = u.Quantity(0.0, u.arcsec), delta: u.Quantity = u.Quantity(0.0, u.arcsec)) -> SkyCoord:
        new_coords = self.chip1_center_coords.spherical_offsets_by(
            alpha + 15/2 * u.arcsecond, delta + 15/2 * u.arcsecond)  # type: ignore
        """returns the CCD4 (chip1) ra and dec center given a pointing in ra dec and rotation"""
        return self.center_coords.directional_offset_by(
            self.center_coords.position_angle(new_coords) + self.rotation, # type: ignore
            self.center_coords.separation(new_coords))
        # return self.center_coords.spherical_offsets_by(self.chip1_center, self.chip1_center)

    @staticmethod
    def _get_pa_sep(coord1: SkyCoord, coord2: SkyCoord) -> tuple[u.Quantity, u.Quantity]:
        pa = coord1.position_angle(coord2)
        if pa is None: raise ValueError('Position angle is None')
        sep = coord1.separation(coord2)
        return pa, sep
