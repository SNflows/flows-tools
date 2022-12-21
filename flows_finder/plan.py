from dataclasses import dataclass
from typing import Optional

import astropy.units as u
import numpy as np

from .catalogs import SkyViewSurveys
from .utils import numeric


@dataclass
class Plan:
    rotation: u.Quantity = u.Quantity(0, u.deg)
    alpha: u.Quantity = u.Quantity(0, u.deg)
    delta: u.Quantity = u.Quantity(0, u.deg)
    rotate: bool = False
    shift: bool = False
    survey: str = SkyViewSurveys.TWO_MASS_H.value  # Planned filter for survey image.
    image_scale: str = "linear"
    local_image: Optional[str] = None

    def __post_init__(self) -> None:
        self.shift = self.set_shift()
        self.rotate = self.set_rotation()

    def plan_finder_chart(
        self, image: Optional[str] = None, survey: Optional[str] = None, scale: str = "linear"
    ) -> None:
        self.local_image = image
        self.image_scale = scale
        if survey is not None:
            self.survey = survey

    def set_rotation(self) -> bool:
        if self.rotation == u.Quantity(0, u.deg):
            return False
        return True

    def set_shift(self) -> bool:
        shift = (np.array((self.alpha.value, self.delta.value)) == 0.0).all()  # skip shift if alpha and delta 0
        return not shift

    @classmethod
    def from_numeric(cls, rotation: numeric, alpha: numeric, delta: numeric) -> "Plan":
        return cls(rotation=rotation << u.deg, alpha=alpha << u.arcsec, delta=delta << u.arcsec)  # type: ignore


def make_plan(
    rotation: numeric = 0.0,
    alpha: numeric = 0.0,
    delta: numeric = 0.0,
    image: Optional[str] = None,
    survey: Optional[str] = SkyViewSurveys.TWO_MASS_H.value,
    scale: str = "linear",
) -> Plan:
    plan = Plan.from_numeric(rotation, alpha, delta)
    plan.plan_finder_chart(image=image, survey=survey, scale=scale)
    return plan
