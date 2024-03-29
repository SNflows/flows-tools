class PixCoord:
    def __init__(self, x, y) -> None: ...
    def copy(self): ...
    @property
    def isscalar(self): ...
    def __len__(self) -> int: ...
    def __iter__(self): ...
    def __getitem__(self, key): ...
    def __add__(self, other): ...
    def __sub__(self, other): ...
    def __eq__(self, other): ...
    def to_sky(self, wcs, origin=..., mode=...): ...
    @classmethod
    def from_sky(cls, skycoord, wcs, origin=..., mode=...): ...
    def separation(self, other): ...
    @property
    def xy(self): ...
    def rotate(self, center, angle): ...
