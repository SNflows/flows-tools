from dataclasses import dataclass
from typing import Any, TypeGuard
import numpy as np
import astropy.units as u
from .utils import numeric

@dataclass
class Plan:
    rotation: u.Quantity = u.Quantity(0, u.deg)
    alpha: u.Quantity = u.Quantity(0, u.deg)
    delta: u.Quantity = u.Quantity(0, u.deg)
    rotate: bool = False
    shift: bool = False

    def __post_init__(self) -> None:
        self.shift = self.set_shift()
        self.rotate = self.set_rotation()
        # self._sanitize_quanity_input()

    def set_rotation(self) -> bool:
        if self.rotation == u.Quantity(0, u.deg):
            return False
        return True

    def set_shift(self) -> bool:
        shift = (np.array((self.alpha.value, self.delta.value)) == 0.0).all()  # skip shift if alpha and delta 0
        return not shift

    # def _sanitize_quanity_input(self) -> None:
    #     self.alpha = self.alpha << u.arcsec  # type: ignore
    #     self.delta = self.delta << u.arcsec  # type: ignore
    #     self.rotation = self.rotation << u.deg  # type: ignore
        
    @classmethod
    def from_numeric(cls, rotation: numeric, alpha: numeric, delta: numeric) -> "Plan":
        return cls(rotation=rotation << u.deg, alpha=alpha << u.arcsec, delta=delta << u.arcsec) # type: ignore
        
    # def finalize(self) -> FinalizedPlan:
    #     if self._finalized(self):
    #         return FinalizedPlan(self.rotation, self.alpha, self.delta, self.rotate, self.shift)
    #     raise ValueError('Plan not finalized, attributes must be quantites.')
    
    # @staticmethod
    # def _finalized(p: Any) -> TypeGuard[FinalizedPlan]:
    #     if isinstance(p.rotation, u.Quantity) and isinstance(p.alpha, u.Quantity) and isinstance(p.delta, u.Quantity):
    #         return True
    #     return False
    