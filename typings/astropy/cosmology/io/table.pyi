"""
This type stub file was generated by pyright.
"""

def from_table(table, index=..., *, move_to_meta=..., cosmology=...):
    """Instantiate a `~astropy.cosmology.Cosmology` from a |QTable|.

    Parameters
    ----------
    table : `~astropy.table.Table`
        The object to parse into a |Cosmology|.
    index : int, str, or None, optional
        Needed to select the row in tables with multiple rows. ``index`` can be
        an integer for the row number or, if the table is indexed by a column,
        the value of that column. If the table is not indexed and ``index``
        is a string, the "name" column is used as the indexing column.

    move_to_meta : bool (optional, keyword-only)
        Whether to move keyword arguments that are not in the Cosmology class'
        signature to the Cosmology's metadata. This will only be applied if the
        Cosmology does NOT have a keyword-only argument (e.g. ``**kwargs``).
        Arguments moved to the metadata will be merged with existing metadata,
        preferring specified metadata in the case of a merge conflict
        (e.g. for ``Cosmology(meta={'key':10}, key=42)``, the ``Cosmology.meta``
        will be ``{'key': 10}``).

    cosmology : str, `~astropy.cosmology.Cosmology` class, or None (optional, keyword-only)
        The cosmology class (or string name thereof) to use when constructing
        the cosmology instance. The class also provides default parameter values,
        filling in any non-mandatory arguments missing in 'table'.

    Returns
    -------
    `~astropy.cosmology.Cosmology` subclass instance

    Examples
    --------
    To see loading a `~astropy.cosmology.Cosmology` from a Table with
    ``from_table``, we will first make a |QTable| using
    :func:`~astropy.cosmology.Cosmology.to_format`.

        >>> from astropy.cosmology import Cosmology, Planck18
        >>> ct = Planck18.to_format("astropy.table")
        >>> ct
        <QTable length=1>
          name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                 km / (Mpc s)            K                 eV
          str8     float64    float64 float64 float64  float64[3] float64
        -------- ------------ ------- ------- ------- ----------- -------
        Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    Now this table can be used to load a new cosmological instance identical
    to the ``Planck18`` cosmology from which it was generated.

        >>> cosmo = Cosmology.from_format(ct, format="astropy.table")
        >>> cosmo
        FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                      Tcmb0=2.7255 K, Neff=3.046, m_nu=[0. 0. 0.06] eV, Ob0=0.04897)

    Specific cosmology classes can be used to parse the data. The class'
    default parameter values are used to fill in any information missing in the
    data.

        >>> from astropy.cosmology import FlatLambdaCDM
        >>> del ct["Tcmb0"]  # show FlatLambdaCDM provides default
        >>> FlatLambdaCDM.from_format(ct)
        FlatLambdaCDM(name="Planck18", H0=67.66 km / (Mpc s), Om0=0.30966,
                      Tcmb0=0.0 K, Neff=3.046, m_nu=None, Ob0=0.04897)

    For tables with multiple rows of cosmological parameters, the ``index``
    argument is needed to select the correct row. The index can be an integer
    for the row number or, if the table is indexed by a column, the value of
    that column. If the table is not indexed and ``index`` is a string, the
    "name" column is used as the indexing column.

    Here is an example where ``index`` is needed and can be either an integer
    (for the row number) or the name of one of the cosmologies, e.g. 'Planck15'.

        >>> from astropy.cosmology import Planck13, Planck15, Planck18
        >>> from astropy.table import vstack
        >>> cts = vstack([c.to_format("astropy.table")
        ...               for c in (Planck13, Planck15, Planck18)],
        ...              metadata_conflicts='silent')
        >>> cts
        <QTable length=3>
          name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                 km / (Mpc s)            K                 eV
          str8     float64    float64 float64 float64  float64[3]  float64
        -------- ------------ ------- ------- ------- ----------- --------
        Planck13        67.77 0.30712  2.7255   3.046 0.0 .. 0.06 0.048252
        Planck15        67.74  0.3075  2.7255   3.046 0.0 .. 0.06   0.0486
        Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06  0.04897

        >>> cosmo = Cosmology.from_format(cts, index=1, format="astropy.table")
        >>> cosmo == Planck15
        True

    For further examples, see :doc:`astropy:cosmology/io`.
    """
    ...

def to_table(cosmology, *args, cls=..., cosmology_in_meta=...): # -> QTable:
    """Serialize the cosmology into a `~astropy.table.QTable`.

    Parameters
    ----------
    cosmology : `~astropy.cosmology.Cosmology` subclass instance
    *args
        Not used. Needed for compatibility with
        `~astropy.io.registry.UnifiedReadWriteMethod`
    cls : type (optional, keyword-only)
        Astropy :class:`~astropy.table.Table` class or subclass type to return.
        Default is :class:`~astropy.table.QTable`.
    cosmology_in_meta : bool
        Whether to put the cosmology class in the Table metadata (if `True`,
        default) or as the first column (if `False`).

    Returns
    -------
    `~astropy.table.QTable`
        With columns for the cosmology parameters, and metadata and
        cosmology class name in the Table's ``meta`` attribute

    Raises
    ------
    TypeError
        If kwarg (optional) 'cls' is not a subclass of `astropy.table.Table`

    Examples
    --------
    A Cosmology as a `~astropy.table.QTable` will have the cosmology's name and
    parameters as columns.

        >>> from astropy.cosmology import Planck18
        >>> ct = Planck18.to_format("astropy.table")
        >>> ct
        <QTable length=1>
          name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                 km / (Mpc s)            K                 eV
          str8     float64    float64 float64 float64  float64[3] float64
        -------- ------------ ------- ------- ------- ----------- -------
        Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    The cosmological class and other metadata, e.g. a paper reference, are in
    the Table's metadata.

        >>> ct.meta
        OrderedDict([..., ('cosmology', 'FlatLambdaCDM')])

    To move the cosmology class from the metadata to a Table row, set the
    ``cosmology_in_meta`` argument to `False`:

        >>> Planck18.to_format("astropy.table", cosmology_in_meta=False)
        <QTable length=1>
          cosmology     name        H0        Om0    Tcmb0    Neff      m_nu      Ob0
                               km / (Mpc s)            K                 eV
            str13       str8     float64    float64 float64 float64  float64[3] float64
        ------------- -------- ------------ ------- ------- ------- ----------- -------
        FlatLambdaCDM Planck18        67.66 0.30966  2.7255   3.046 0.0 .. 0.06 0.04897

    Astropy recommends `~astropy.table.QTable` for tables with
    `~astropy.units.Quantity` columns. However the returned type may be
    overridden using the ``cls`` argument:

        >>> from astropy.table import Table
        >>> Planck18.to_format("astropy.table", cls=Table)
        <Table length=1>
        ...
    """
    ...

def table_identify(origin, format, *args, **kwargs): # -> bool:
    """Identify if object uses the Table format.

    Returns
    -------
    bool
    """
    ...
