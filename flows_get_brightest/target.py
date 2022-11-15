from dataclasses import dataclass, field
from typing import Any
from .utils import MISSING_COORDS
from astropy.coordinates import SkyCoord

@dataclass
class Target:
    tid: int | str  # target id
    ra: float
    dec: float
    coords: SkyCoord = MISSING_COORDS
    info: dict[str, Any] = field(default_factory=dict)
    name: str = "UNKNOWN"
    loaded: bool = field(init=False, default=False)

    def __post_init__(self) -> None:
        if self.coords is MISSING_COORDS or self.coords is None:
            c = SkyCoord(self.ra << u.deg, self.dec << u.deg, frame='icrs')  # type: ignore
            self.coords = c  
        if self.name == "UNKNOWN":
            self.name = self.info.get('target_name', f"Target {self.tid}")
        if isinstance(self.tid, str):
            if self.tid.isnumeric():
                self.tid = int(self.tid)
            else:
                tid = self.info.get('target_id', None)
                if tid:
                    self.tid = int(tid)
                
        self.ra = float(self.coords.ra.deg)  # type: ignore
        self.deg = float(self.coords.dec.deg)  # type: ignore
        self.loaded = True
        