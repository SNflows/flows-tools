from dataclasses import dataclass
import numpy as np
import astropy.units as u
from .utils import numeric

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
        self._sanitize_quanity_input()

    def set_rotation(self) -> bool:
        if self.rotation == 0.0:
            return False
        elif self.rotation != 0.0:
            return True

    def set_shift(self) -> bool:
        shift = (np.array((self.alpha, self.delta)) == 0.0).all()  # skip shift if alpha and delta 0
        return not shift

    def _sanitize_quanity_input(self):
        self.alpha = self.alpha << u.arcsec
        self.delta = self.delta << u.arcsec
        self.rotation = self.rotation << u.deg