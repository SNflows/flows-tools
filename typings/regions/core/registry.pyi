from _typeshed import Incomplete

class IORegistryError(Exception): ...

class RegionsRegistry:
    registry: Incomplete
    @classmethod
    def register(cls, classobj, methodname, filetype): ...
    @classmethod
    def get_identifiers(cls, classobj): ...
    @classmethod
    def identify_format(cls, filename, classobj, methodname): ...
    @classmethod
    def read(cls, filename, classobj, format: Incomplete | None = ..., **kwargs): ...
    @classmethod
    def parse(cls, data, classobj, format: Incomplete | None = ..., **kwargs): ...
    @classmethod
    def write(cls, regions, filename, classobj, format: Incomplete | None = ..., **kwargs): ...
    @classmethod
    def serialize(cls, regions, classobj, format: Incomplete | None = ..., **kwargs): ...
    @classmethod
    def get_formats(cls, classobj): ...