from _typeshed import Incomplete

class Regions:
    regions: Incomplete
    def __init__(self, regions) -> None: ...
    def __getitem__(self, index): ...
    def __len__(self) -> int: ...
    def append(self, region) -> None: ...
    def extend(self, regions) -> None: ...
    def insert(self, index, region) -> None: ...
    def reverse(self) -> None: ...
    def pop(self, index: int = ...): ...
    def copy(self): ...
    def plot(self): ...
    @classmethod
    def get_formats(cls): ...
    @classmethod
    def read(cls, filename, format: Incomplete | None = ..., cache: bool = ..., **kwargs): ...
    @classmethod
    def parse(cls, data, format: Incomplete | None = ..., **kwargs): ...
    def write(self, filename, format: Incomplete | None = ..., overwrite: bool = ..., **kwargs): ...
    def serialize(self, format: Incomplete | None = ..., **kwargs): ...
