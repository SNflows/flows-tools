from typing import Optional, Union

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
from astropy.wcs import WCS
from astropy.wcs.utils import celestial_frame_to_wcs
from astropy.io.fits import PrimaryHDU
import regions
import numpy as np
from numpy.typing import NDArray
from tendrils import api

from .instruments import Instrument, Hawki
from .target import Target
from .utils import numeric, tabular
from .catalogs import query_2mass_image, query_simbad
from .plan import Plan


class Observer:
    def __init__(self, instrument: Instrument, target: Target, plan: Plan, verbose: bool = False):
        self.target = target  # target info
        self.plan = plan  # Store offsets and rotations.
        self.ins = instrument  # store instrument specific info
        self.regions = self.ins.point(target, self.plan)
        self.nregions = len(self.regions)  # How many sub-regions of important passed by Instrument.?

        self.refcat_coords, self.refcat = self._make_refcat_catalog(mask=True, mask_dict={"H_mag": 14.0})
        self.simbad_coords, self.simbad = self._make_simbad_catalog()

        if not verbose:
            print(
                "Simbad and refcat skyframes are "
                f"equivalent: {self.refcat_coords.frame.is_equivalent_frame(self.simbad_coords.frame)}"
                ". If not, there might be a slight mismatch in alignment in current year."
            )
        self.wcs = self.get_wcs(self.refcat_coords.frame)

    def _make_refcat_catalog(self, mask: bool = True, mask_dict: Optional[dict[str, float]] = None) -> tabular:
        refcat = api.get_catalog(self.target.tid)
        if mask and mask_dict is not None:
            masks = []
            keys_values = mask_dict.items()
            for key, value in keys_values:
                masks.append(refcat["references"][str(key)] <= float(value))
            if len(masks) == 1:
                refcat["references"] = refcat["references"][masks[0]]
            else:
                for msk in masks:
                    refcat["references"] = refcat["references"][msk]
        refcat_coords = SkyCoord(
            ra=refcat["references"]["ra"],
            dec=refcat["references"]["decl"],
            pm_ra_cosdec=refcat["references"]["pm_ra"],
            pm_dec=refcat["references"]["pm_dec"],
            obstime=Time(2015.5, format="decimalyear"),
        )
        return refcat_coords, refcat["references"]

    def _make_simbad_catalog(self) -> tabular:
        simbad, coords = query_simbad(self.target.coords)
        if simbad is None or coords is None:
            raise ValueError("No Simbad catalog found for target.")
        # propagate Simbad catalog coords to refcat reference year
        simbad_coords = coords.apply_space_motion(new_obstime=Time(2015.5, format="decimalyear"))
        return simbad_coords, simbad

    def get_wcs(self, frame: Optional[object] = None, header: Optional[object] = None) -> WCS:
        if frame is not None:
            wcs = celestial_frame_to_wcs(frame)
            wcs.wcs.crval = np.array((self.target.ra, self.target.dec))
            return wcs
        elif header is not None:
            return WCS(header)
        else:
            raise AttributeError("Both frame and header were None, cannot calculate WCS")

    def get_image(self, pixels: int = 2500, radius: numeric = 50) -> PrimaryHDU:
        return query_2mass_image(self.target.ra, self.target.dec, pixels, radius)

    def regions_to_physical(self) -> list[regions.RectangleSkyRegion]:
        pixel_regions = []
        for region in self.regions:
            pixel_regions.append(region.to_pixel(self.wcs))
        return pixel_regions

    def check_bright_stars(self, region: regions.RectangleSkyRegion, wcs: Optional[WCS] = None) -> np.ndarray:
        wcs = self.wcs if wcs is None else wcs
        # Check bright stars
        simbad_stars: NDArray[np.bool_] = region.contains(self.simbad_coords, wcs)
        catalog_stars: NDArray[np.bool_] = region.contains(self.refcat_coords, wcs)

        simbad_mag: NDArray[np.float_] = self.simbad[simbad_stars]["H_mag"].data  # type: ignore
        refcat_mag: NDArray[np.float_] = self.refcat[catalog_stars]["H_mag"].data  # type: ignore
        bright_stars = np.hstack((simbad_mag, refcat_mag))
        if np.ma.isMaskedArray(bright_stars):
            bright_stars: NDArray[np.float_] = bright_stars.data[~np.isnan(bright_stars.data)]  # type: ignore

        print("Brightest star has H-mag = {0:.1f}".format(np.round(bright_stars.min(), 1)))
        return bright_stars


def get_flows_observer(
    rot: numeric, tid: Union[int, str], shifta: numeric, shiftd: numeric, instrument: type[Instrument] = Hawki
) -> Observer:
    """
    Returns the H-mag of the brightest star in the given region.
    """
    # Get query flows database
    target_info = api.get_target(tid)
    target_info["skycoord"] = SkyCoord(target_info["ra"] * u.deg, target_info["decl"] * u.deg)

    # Create Observer
    target = Target(tid, target_info["ra"], target_info["decl"], target_info["skycoord"], target_info)
    plan = Plan.from_numeric(rot, shifta, shiftd)
    hawki = instrument(target.coords)
    return Observer(hawki, target, plan)
