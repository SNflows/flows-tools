"""
This type stub file was generated by pyright.
"""

from ..query import BaseQuery
from ..utils import async_to_sync, prepend_docstr_nosections

__doctest_skip__ = ...
@async_to_sync
class SkyViewClass(BaseQuery):
    URL = ...
    def __init__(self) -> None:
        ...
    
    def get_images(self, position, survey, coordinates=..., projection=..., pixels=..., scaling=..., sampler=..., resolver=..., deedger=..., lut=..., grid=..., gridlabels=..., radius=..., height=..., width=..., cache=..., show_progress=...): # -> list[Unknown]:
        """
        Query the SkyView service, download the FITS file that will be
        found and return a generator over the local paths to the
        downloaded FITS files.

        Note that the files will be downloaded when the generator will be
        exhausted, i.e. just calling this method alone without iterating
        over the result won't issue a connection to the SkyView server.

        Parameters
        ----------
        position : str
            Determines the center of the field to be retrieved. Both
            coordinates (also equatorial ones) and object names are
            supported. Object names are converted to coordinates via the
            SIMBAD or NED name resolver. See the reference for more info
            on the supported syntax for coordinates.
        survey : str or list of str
            Select data from one or more surveys. The number of surveys
            determines the number of resulting file downloads. Passing a
            list with just one string has the same effect as passing this
            string directly.
        coordinates : str
            Choose among common equatorial, galactic and ecliptic
            coordinate systems (``"J2000"``, ``"B1950"``, ``"Galactic"``,
            ``"E2000"``, ``"ICRS"``) or pass a custom string.
        projection : str
            Choose among the map projections (the value in parentheses
            denotes the string to be passed):

            Gnomonic (Tan), default value
                good for small regions
            Rectangular (Car)
                simplest projection
            Aitoff (Ait)
                Hammer-Aitoff, equal area projection good for all sky maps
            Orthographic (Sin)
                Projection often used in interferometry
            Zenith Equal Area (Zea)
                equal area, azimuthal projection
            COBE Spherical Cube (Csc)
                Used in COBE data
            Arc (Arc)
                Similar to Zea but not equal-area
        pixels : str
            Selects the pixel dimensions of the image to be produced. A
            scalar value or a pair of values separated by comma may be
            given. If the value is a scalar the number of width and height
            of the image will be the same. By default a 300x300 image is
            produced.
        scaling : str
            Selects the transformation between pixel intensity and
            intensity on the displayed image. The supported values are:
            ``"Log"``, ``"Sqrt"``, ``"Linear"``, ``"HistEq"``,
            ``"LogLog"``.
        sampler : str
            The sampling algorithm determines how the data requested will
            be resampled so that it can be displayed.
        resolver : str
            The name resolver allows to choose a name resolver to use when
            looking up a name which was passed in the ``position`` parameter
            (as opposed to a numeric coordinate value). The default choice
            is to call the SIMBAD name resolver first and then the NED
            name resolver if the SIMBAD search fails.
        deedger : str
            When multiple input images with different backgrounds are
            resampled the edges between the images may be apparent because
            of the background shift. This parameter makes it possible to
            attempt to minimize these edges by applying a de-edging
            algorithm. The user can elect to choose the default given for
            that survey, to turn de-edging off, or to use the default
            de-edging algorithm. The supported values are: ``"_skip_"`` to
            use the survey default, ``"skyview.process.Deedger"`` (for
            enabling de-edging), and ``"null"`` to disable.
        lut : str
            Choose from the color table selections to display the data in
            false color.
        grid : bool
            overlay a coordinate grid on the image if True
        gridlabels : bool
            annotate the grid with coordinates positions if True
        radius : `~astropy.units.Quantity` or None
            The radius of the specified field.  Overrides width and height.
        width : `~astropy.units.Quantity` or None
            The width of the specified field.  Must be specified
            with ``height``.
        height : `~astropy.units.Quantity` or None
            The height of the specified field.  Must be specified
            with ``width``.

        References
        ----------
        .. [1] http://skyview.gsfc.nasa.gov/current/help/fields.html

        Examples
        --------
        >>> sv = SkyView()
        >>> paths = sv.get_images(position='Eta Carinae',
        ...                       survey=['Fermi 5', 'HRI', 'DSS'])
        >>> for path in paths:
        ...     print('\tnew file:', path)

        Returns
        -------
        A list of `~astropy.io.fits.HDUList` objects.

        """
        ...
    
    @prepend_docstr_nosections(get_images.__doc__)
    def get_images_async(self, position, survey, coordinates=..., projection=..., pixels=..., scaling=..., sampler=..., resolver=..., deedger=..., lut=..., grid=..., gridlabels=..., radius=..., height=..., width=..., cache=..., show_progress=...): # -> list[FileContainer]:
        """
        Returns
        -------
        A list of context-managers that yield readable file-like objects
        """
        ...
    
    @prepend_docstr_nosections(get_images.__doc__, sections=['Returns', 'Examples'])
    def get_image_list(self, position, survey, coordinates=..., projection=..., pixels=..., scaling=..., sampler=..., resolver=..., deedger=..., lut=..., grid=..., gridlabels=..., radius=..., width=..., height=..., cache=...): # -> list[Unknown]:
        """
        Returns
        -------
        list of image urls

        Examples
        --------
        >>> SkyView().get_image_list(position='Eta Carinae',
        ...                          survey=['Fermi 5', 'HRI', 'DSS'])
        [u'http://skyview.gsfc.nasa.gov/tempspace/fits/skv6183161285798_1.fits',
         u'http://skyview.gsfc.nasa.gov/tempspace/fits/skv6183161285798_2.fits',
         u'http://skyview.gsfc.nasa.gov/tempspace/fits/skv6183161285798_3.fits']
        """
        ...
    
    @property
    def survey_dict(self): # -> dict[Any, list[Any]]:
        ...
    
    def list_surveys(self): # -> None:
        """
        Print out a formatted version of the survey dict
        """
        ...
    


def parse_coordinates(position):
    ...

SkyView = SkyViewClass()
