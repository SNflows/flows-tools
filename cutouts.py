import astropy.io.fits as fits
# noinspection PyProtectedMember
from astropy.io.fits.hdu.image import _ImageBaseHDU
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from pathlib import Path
import argparse
from astropy.visualization.interval import BaseInterval
import numpy as np
from abc import ABC, abstractmethod

"""Expects fits compatible images saved with a single file extension (no .fits.gz), 
Makes a cutout based on values of y and x and recomputes wcs,
Saves the new fits file in output dir, 
preserves all extensions and applies the same cut to all image extensions
optionally plots the primary image in cut and uncut form.

code_author: emirkmo@github.com.
"""


class AbstractImage(ABC):

    @abstractmethod
    def __init__(self, image: _ImageBaseHDU):
        self.image = image
        self.ImType = type(image)
        self.wcs_init = None
        self.cut_image = None

    @abstractmethod
    def cutout(self, y: int, x: int):
        pass

    @abstractmethod
    def make_cutout(self, y: int, x: int):
        pass


# noinspection PyAttributeOutsideInit
class Image(AbstractImage):

    def __init__(self, image: _ImageBaseHDU,
                 primary: bool = False):
        self.image = image
        self.data = image.data
        self.header = image.header
        self.wcs_init = WCS(self.image.header)
        self.isprimary = primary
        self.ImType = type(image)
        if self.isprimary:
            self.wcs = self.wcs_init

    # def get_class(self):
    #    if isinstance(self.image,fits.ImageHDU):
    #        self.ImType = fits.ImageHDU
    #        self.isprimary = False
    #    elif isinstance(self.image,fits.PrimaryHDU):
    #        self.ImType = fits.PrimaryHDU
    #        self.isprimary = True
    #    else: raise(TypeError)

    def cutout(self, y, x):
        self.cut = Cutout2D(
            self.data,
            position=self.wcs.wcs.crpix,
            wcs=self.wcs,
            size=(y, x)
        )

    # noinspection PyCallingNonCallable
    def make_fits_from_cutout(self):
        self.cut_image = self.ImType(self.cut.data, self.cut.wcs.to_header())

    def find_missing_keys(self):
        self.missing_keys = set(self.header.keys()) - set(self.cut_image.header.keys())

    def fill_cut_image_header(self):
        for key in self.missing_keys:
            try:
                self.cut_image.header[key] = self.header[key]
            except ValueError:
                print(f'Bad fits key: {key}')

    def make_cutout(self, y: int, x: int):
        """convenience function for making a cutout image and forward filling header"""
        self.cutout(y, x)
        self.make_fits_from_cutout()
        self.find_missing_keys()
        self.fill_cut_image_header()


class Plotter:

    def __init__(self,
                 figsize: tuple = (8, 8),
                 zs: BaseInterval = ZScaleInterval(1000),
                 subplots: int = 1,
                 dpi: int = 125
                 ):
        self.figsize = figsize
        self.interval_finder = zs
        self.fig = plt.figure(figsize=self.figsize, dpi=dpi)
        self.subplots = subplots
        self.active_axis = 1
        self.axes = []
        self.ax = None

    # noinspection SpellCheckingInspection
    def plot_img(self, image: Image, wcs: WCS = 'None'):
        if wcs == 'None':
            wcs = image.wcs
        self.ax = self.fig.add_subplot(1, self.subplots, self.active_axis, projection=wcs)
        vmin, vmax = self.interval_finder.get_limits(image.data)
        self.ax.imshow(image.data, vmin=vmin, vmax=vmax)
        self.ax.grid()
        self.active_axis += 1
        self.axes.append(self.ax)


def parse():
    """Parse input
    :rtype: object
    """
    parser = argparse.ArgumentParser(prog='cut_image',
                                     description='Cut image and all image extensions')

    parser.add_argument("file", help="Specify the original file here")

    parser.add_argument("--output_dir", help="Give output directory for saving", default='./')
    parser.add_argument("-x", help='output suffix', type=int, default=4096)
    parser.add_argument("-y", help='output suffix', type=int, default=4030)
    parser.add_argument("-o", "--overwrite", help='overwrite output if exists',
                        action='store_true', default=False)
    parser.add_argument("-p", "--plot", help='overwrite output if exists',
                        action='store_true', default=False)

    return parser.parse_args()


def find_primary(hdul):
    # Find primary HDU and use its wcs
    primary_index = np.argwhere([isinstance(i, fits.PrimaryHDU) for i in hdul]).flatten()
    if primary_index.size <= 0 or primary_index.size > 1:
        raise ValueError('WCS to use could not be set as there were less than or greater than one Primary HDUs \
                         This should never be the case, please check your fits image!')
    return primary_index[0]


def get_primary(hdul) -> AbstractImage:
    return Image(hdul[find_primary(hdul)], primary=True)


def make_img_cutout(i, wcs, y, x) -> AbstractImage:
    img = Image(i, primary=False)
    img.wcs = wcs
    img.make_cutout(y, x)
    return img


def main():
    args = parse()
    curr = Path(args.file)
    save_path = Path.joinpath(Path(args.output_dir), curr.stem + '_cutout' + curr.suffix).resolve()

    hdul = fits.open(curr)
    new_hdul = fits.HDUList()

    primary_image = get_primary(hdul)
    primary_image.make_cutout(args.y, args.x)
    new_hdul.append(primary_image.cut_image)

    for i in hdul:
        if isinstance(i, fits.PrimaryHDU):
            continue
        elif isinstance(i, fits.ImageHDU):
            img = make_img_cutout(i, primary_image.wcs_init, args.y, args.x)
            new_hdul.append(img.cut_image)
        else:  # Mostly expecting extension to be a BinaryFitsTable.
            new_hdul.append(i)

    if args.plot:
        plotter = Plotter(subplots=2, dpi=200)
        plotter.plot_img(Image(new_hdul[0], primary=True))
        plotter.ax.set_title('cutout')
        plotter.plot_img(Image(hdul[0], primary=True))
        plotter.ax.set_title('uncut')
        plt.show(block=True)

    # Save
    new_hdul.writeto(str(save_path), output_verify='exception',
                     overwrite=args.overwrite, checksum=False)
    print(f'saved to: {save_path}')


if __name__ == '__main__':
    main()
