from dataclasses import dataclass, field
from typing import Optional
from .utils import numeric
from astropy.coordinates import SkyCoord
import astropy.units as u


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
