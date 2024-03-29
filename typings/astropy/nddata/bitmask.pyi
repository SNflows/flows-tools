"""
This type stub file was generated by pyright.
"""

import numpy as np

"""
A module that provides functions for manipulating bit masks and data quality
(DQ) arrays.

"""
__all__ = ['bitfield_to_boolean_mask', 'interpret_bit_flags', 'BitFlagNameMap', 'extend_bit_flag_map', 'InvalidBitFlag']
_ENABLE_BITFLAG_CACHING = ...
_MAX_UINT_TYPE = np.maximum_sctype(np.uint)
_SUPPORTED_FLAGS = ...
class InvalidBitFlag(ValueError):
    """ Indicates that a value is not an integer that is a power of 2. """
    ...


class BitFlag(int):
    """ Bit flags: integer values that are powers of 2. """
    def __new__(cls, val, doc=...): # -> Self@BitFlag:
        ...
    


class BitFlagNameMeta(type):
    def __new__(mcls, name, bases, members): # -> Self@BitFlagNameMeta:
        ...
    
    def __setattr__(cls, name, val): # -> None:
        ...
    
    def __getattr__(cls, name): # -> Any:
        ...
    
    def __getitem__(cls, key): # -> Any:
        ...
    
    def __add__(cls, items): # -> BitFlagNameMeta:
        ...
    
    def __iadd__(cls, other):
        ...
    
    def __delattr__(cls, name):
        ...
    
    def __delitem__(cls, name):
        ...
    
    def __repr__(cls): # -> str:
        ...
    


class BitFlagNameMap(metaclass=BitFlagNameMeta):
    """
    A base class for bit flag name maps used to describe data quality (DQ)
    flags of images by provinding a mapping from a mnemonic flag name to a flag
    value.

    Mapping for a specific instrument should subclass this class.
    Subclasses should define flags as class attributes with integer values
    that are powers of 2. Each bit flag may also contain a string
    comment following the flag value.

    Examples
    --------

        >>> from astropy.nddata.bitmask import BitFlagNameMap
        >>> class ST_DQ(BitFlagNameMap):
        ...     __version__ = '1.0.0'  # optional
        ...     CR = 1, 'Cosmic Ray'
        ...     CLOUDY = 4  # no docstring comment
        ...     RAINY = 8, 'Dome closed'
        ...
        >>> class ST_CAM1_DQ(ST_DQ):
        ...     HOT = 16
        ...     DEAD = 32

    """
    ...


def extend_bit_flag_map(cls_name, base_cls=..., **kwargs): # -> BitFlagNameMeta:
    """
    A convenience function for creating bit flags maps by subclassing an
    existing map and adding additional flags supplied as keyword arguments.

    Parameters
    ----------
    cls_name : str
        Class name of the bit flag map to be created.

    base_cls : BitFlagNameMap, optional
        Base class for the new bit flag map.

    **kwargs : int
        Each supplied keyword argument will be used to define bit flag
        names in the new map. In addition to bit flag names, ``__version__`` is
        allowed to indicate the version of the newly created map.

    Examples
    --------

        >>> from astropy.nddata.bitmask import extend_bit_flag_map
        >>> ST_DQ = extend_bit_flag_map('ST_DQ', __version__='1.0.0', CR=1, CLOUDY=4, RAINY=8)
        >>> ST_CAM1_DQ = extend_bit_flag_map('ST_CAM1_DQ', ST_DQ, HOT=16, DEAD=32)
        >>> ST_CAM1_DQ['HOT']  # <-- Access flags as dictionary keys
        16
        >>> ST_CAM1_DQ.HOT  # <-- Access flags as class attributes
        16

    """
    ...

def interpret_bit_flags(bit_flags, flip_bits=..., flag_name_map=...): # -> int | None:
    """
    Converts input bit flags to a single integer value (bit mask) or `None`.

    When input is a list of flags (either a Python list of integer flags or a
    string of comma-, ``'|'``-, or ``'+'``-separated list of flags),
    the returned bit mask is obtained by summing input flags.

    .. note::
        In order to flip the bits of the returned bit mask,
        for input of `str` type, prepend '~' to the input string. '~' must
        be prepended to the *entire string* and not to each bit flag! For
        input that is already a bit mask or a Python list of bit flags, set
        ``flip_bits`` for `True` in order to flip the bits of the returned
        bit mask.

    Parameters
    ----------
    bit_flags : int, str, list, None
        An integer bit mask or flag, `None`, a string of comma-, ``'|'``- or
        ``'+'``-separated list of integer bit flags or mnemonic flag names,
        or a Python list of integer bit flags. If ``bit_flags`` is a `str`
        and if it is prepended with '~', then the output bit mask will have
        its bits flipped (compared to simple sum of input flags).
        For input ``bit_flags`` that is already a bit mask or a Python list
        of bit flags, bit-flipping can be controlled through ``flip_bits``
        parameter.

        .. note::
            When ``bit_flags`` is a list of flag names, the ``flag_name_map``
            parameter must be provided.

        .. note::
            Only one flag separator is supported at a time. ``bit_flags``
            string should not mix ``','``, ``'+'``, and ``'|'`` separators.

    flip_bits : bool, None
        Indicates whether or not to flip the bits of the returned bit mask
        obtained from input bit flags. This parameter must be set to `None`
        when input ``bit_flags`` is either `None` or a Python list of flags.

    flag_name_map : BitFlagNameMap
         A `BitFlagNameMap` object that provides mapping from mnemonic
         bit flag names to integer bit values in order to translate mnemonic
         flags to numeric values when ``bit_flags`` that are comma- or
         '+'-separated list of menmonic bit flag names.

    Returns
    -------
    bitmask : int or None
        Returns an integer bit mask formed from the input bit value or `None`
        if input ``bit_flags`` parameter is `None` or an empty string.
        If input string value was prepended with '~' (or ``flip_bits`` was set
        to `True`), then returned value will have its bits flipped
        (inverse mask).

    Examples
    --------

        >>> from astropy.nddata.bitmask import interpret_bit_flags, extend_bit_flag_map
        >>> ST_DQ = extend_bit_flag_map('ST_DQ', CR=1, CLOUDY=4, RAINY=8, HOT=16, DEAD=32)
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags(28))
        '0000000000011100'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('4,8,16'))
        '0000000000011100'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('CLOUDY,RAINY,HOT', flag_name_map=ST_DQ))
        '0000000000011100'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('~4,8,16'))
        '1111111111100011'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('~(4+8+16)'))
        '1111111111100011'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags('~(CLOUDY+RAINY+HOT)',
        ... flag_name_map=ST_DQ))
        '1111111111100011'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags([4, 8, 16]))
        '0000000000011100'
        >>> "{0:016b}".format(0xFFFF & interpret_bit_flags([4, 8, 16], flip_bits=True))
        '1111111111100011'

    """
    ...

def bitfield_to_boolean_mask(bitfield, ignore_flags=..., flip_bits=..., good_mask_value=..., dtype=..., flag_name_map=...): # -> NDArray[bool_]:
    """
    bitfield_to_boolean_mask(bitfield, ignore_flags=None, flip_bits=None, \
good_mask_value=False, dtype=numpy.bool_)
    Converts an array of bit fields to a boolean (or integer) mask array
    according to a bit mask constructed from the supplied bit flags (see
    ``ignore_flags`` parameter).

    This function is particularly useful to convert data quality arrays to
    boolean masks with selective filtering of DQ flags.

    Parameters
    ----------
    bitfield : ndarray
        An array of bit flags. By default, values different from zero are
        interpreted as "bad" values and values equal to zero are considered
        as "good" values. However, see ``ignore_flags`` parameter on how to
        selectively ignore some bits in the ``bitfield`` array data.

    ignore_flags : int, str, list, None (default = 0)
        An integer bit mask, `None`, a Python list of bit flags, a comma-,
        or ``'|'``-separated, ``'+'``-separated string list of integer
        bit flags or mnemonic flag names that indicate what bits in the input
        ``bitfield`` should be *ignored* (i.e., zeroed), or `None`.

        .. note::
            When ``bit_flags`` is a list of flag names, the ``flag_name_map``
            parameter must be provided.

        | Setting ``ignore_flags`` to `None` effectively will make
          `bitfield_to_boolean_mask` interpret all ``bitfield`` elements
          as "good" regardless of their value.

        | When ``ignore_flags`` argument is an integer bit mask, it will be
          combined using bitwise-NOT and bitwise-AND with each element of the
          input ``bitfield`` array (``~ignore_flags & bitfield``). If the
          resultant bitfield element is non-zero, that element will be
          interpreted as a "bad" in the output boolean mask and it will be
          interpreted as "good" otherwise. ``flip_bits`` parameter may be used
          to flip the bits (``bitwise-NOT``) of the bit mask thus effectively
          changing the meaning of the ``ignore_flags`` parameter from "ignore"
          to "use only" these flags.

        .. note::

            Setting ``ignore_flags`` to 0 effectively will assume that all
            non-zero elements in the input ``bitfield`` array are to be
            interpreted as "bad".

        | When ``ignore_flags`` argument is a Python list of integer bit
          flags, these flags are added together to create an integer bit mask.
          Each item in the list must be a flag, i.e., an integer that is an
          integer power of 2. In order to flip the bits of the resultant
          bit mask, use ``flip_bits`` parameter.

        | Alternatively, ``ignore_flags`` may be a string of comma- or
          ``'+'``(or ``'|'``)-separated list of integer bit flags that should
          be added (bitwise OR) together to create an integer bit mask.
          For example, both ``'4,8'``, ``'4|8'``, and ``'4+8'`` are equivalent
          and indicate that bit flags 4 and 8 in the input ``bitfield``
          array should be ignored when generating boolean mask.

        .. note::

            ``'None'``, ``'INDEF'``, and empty (or all white space) strings
            are special values of string ``ignore_flags`` that are
            interpreted as `None`.

        .. note::

            Each item in the list must be a flag, i.e., an integer that is an
            integer power of 2. In addition, for convenience, an arbitrary
            **single** integer is allowed and it will be interpreted as an
            integer bit mask. For example, instead of ``'4,8'`` one could
            simply provide string ``'12'``.

        .. note::
            Only one flag separator is supported at a time. ``ignore_flags``
            string should not mix ``','``, ``'+'``, and ``'|'`` separators.

        .. note::

            When ``ignore_flags`` is a `str` and when it is prepended with
            '~', then the meaning of ``ignore_flags`` parameters will be
            reversed: now it will be interpreted as a list of bit flags to be
            *used* (or *not ignored*) when deciding which elements of the
            input ``bitfield`` array are "bad". Following this convention,
            an ``ignore_flags`` string value of ``'~0'`` would be equivalent
            to setting ``ignore_flags=None``.

        .. warning::

            Because prepending '~' to a string ``ignore_flags`` is equivalent
            to setting ``flip_bits`` to `True`, ``flip_bits`` cannot be used
            with string ``ignore_flags`` and it must be set to `None`.

    flip_bits : bool, None (default = None)
        Specifies whether or not to invert the bits of the bit mask either
        supplied directly through ``ignore_flags`` parameter or built from the
        bit flags passed through ``ignore_flags`` (only when bit flags are
        passed as Python lists of integer bit flags). Occasionally, it may be
        useful to *consider only specific bit flags* in the ``bitfield``
        array when creating a boolean mask as opposed to *ignoring* specific
        bit flags as ``ignore_flags`` behaves by default. This can be achieved
        by inverting/flipping the bits of the bit mask created from
        ``ignore_flags`` flags which effectively changes the meaning of the
        ``ignore_flags`` parameter from "ignore" to "use only" these flags.
        Setting ``flip_bits`` to `None` means that no bit flipping will be
        performed. Bit flipping for string lists of bit flags must be
        specified by prepending '~' to string bit flag lists
        (see documentation for ``ignore_flags`` for more details).

        .. warning::
            This parameter can be set to either `True` or `False` **ONLY** when
            ``ignore_flags`` is either an integer bit mask or a Python
            list of integer bit flags. When ``ignore_flags`` is either
            `None` or a string list of flags, ``flip_bits`` **MUST** be set
            to `None`.

    good_mask_value : int, bool (default = False)
        This parameter is used to derive the values that will be assigned to
        the elements in the output boolean mask array that correspond to the
        "good" bit fields (that are 0 after zeroing bits specified by
        ``ignore_flags``) in the input ``bitfield`` array. When
        ``good_mask_value`` is non-zero or ``numpy.True_`` then values in the
        output boolean mask array corresponding to "good" bit fields in
        ``bitfield`` will be ``numpy.True_`` (if ``dtype`` is ``numpy.bool_``)
        or 1 (if ``dtype`` is of numerical type) and values of corresponding
        to "bad" flags will be ``numpy.False_`` (or 0). When
        ``good_mask_value`` is zero or ``numpy.False_`` then the values
        in the output boolean mask array corresponding to "good" bit fields
        in ``bitfield`` will be ``numpy.False_`` (if ``dtype`` is
        ``numpy.bool_``) or 0 (if ``dtype`` is of numerical type) and values
        of corresponding to "bad" flags will be ``numpy.True_`` (or 1).

    dtype : data-type (default = ``numpy.bool_``)
        The desired data-type for the output binary mask array.

    flag_name_map : BitFlagNameMap
         A `BitFlagNameMap` object that provides mapping from mnemonic
         bit flag names to integer bit values in order to translate mnemonic
         flags to numeric values when ``bit_flags`` that are comma- or
         '+'-separated list of menmonic bit flag names.

    Returns
    -------
    mask : ndarray
        Returns an array of the same dimensionality as the input ``bitfield``
        array whose elements can have two possible values,
        e.g., ``numpy.True_`` or ``numpy.False_`` (or 1 or 0 for integer
        ``dtype``) according to values of to the input ``bitfield`` elements,
        ``ignore_flags`` parameter, and the ``good_mask_value`` parameter.

    Examples
    --------

        >>> from astropy.nddata import bitmask
        >>> import numpy as np
        >>> dqarr = np.asarray([[0, 0, 1, 2, 0, 8, 12, 0],
        ...                     [10, 4, 0, 0, 0, 16, 6, 0]])
        >>> flag_map = bitmask.extend_bit_flag_map(
        ...     'ST_DQ', CR=2, CLOUDY=4, RAINY=8, HOT=16, DEAD=32
        ... )
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags=0,
        ...                                  dtype=int)
        array([[0, 0, 1, 1, 0, 1, 1, 0],
               [1, 1, 0, 0, 0, 1, 1, 0]])
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags=0,
        ...                                  dtype=bool)
        array([[False, False,  True,  True, False,  True,  True, False],
               [ True,  True, False, False, False,  True,  True, False]]...)
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags=6,
        ...                                  good_mask_value=0, dtype=int)
        array([[0, 0, 1, 0, 0, 1, 1, 0],
               [1, 0, 0, 0, 0, 1, 0, 0]])
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags=~6,
        ...                                  good_mask_value=0, dtype=int)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags=6, dtype=int,
        ...                                  flip_bits=True, good_mask_value=0)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags='~(2+4)',
        ...                                  good_mask_value=0, dtype=int)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags=[2, 4],
        ...                                  flip_bits=True, good_mask_value=0,
        ...                                  dtype=int)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags='~(CR,CLOUDY)',
        ...                                  good_mask_value=0, dtype=int,
        ...                                  flag_name_map=flag_map)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])
        >>> bitmask.bitfield_to_boolean_mask(dqarr, ignore_flags='~(CR+CLOUDY)',
        ...                                  good_mask_value=0, dtype=int,
        ...                                  flag_name_map=flag_map)
        array([[0, 0, 0, 1, 0, 0, 1, 0],
               [1, 1, 0, 0, 0, 0, 1, 0]])

    """
    ...

