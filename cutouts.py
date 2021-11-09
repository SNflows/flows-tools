import astropy.io.fits as fits
from astropy.nddata import Cutout2D
from astropy import units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from pathlib import Path
import argparse
import typing
from astropy.visualization.interval import BaseInterval
        
"""Expects fits compatible images saved with a single file extension (no .fits.gz), 
Makes a cutout based on values of y and x and recomputes wcs,
Saves the new fits file in output dir, 
preserves all extensions and applies the same cut to all image extensions
optionally plots the primary image in cut and uncut form.

code_author: emirkmo@github.com.
"""

class Image:
    
    def __init__(self, image: typing.Union[fits.PrimaryHDU, fits.ImageHDU]):
        self.image = image
        self.get_class()
        self.data = image.data
        self.header = image.header
        self.wcs_init = self.update_wcs()
        
    def get_class(self):
        if isinstance(self.image,fits.ImageHDU):
            self.ImType = fits.ImageHDU
        elif isinstance(self.image,fits.PrimaryHDU):
            self.ImType = fits.PrimaryHDU
        else: raise(TypeError)
        
    def update_wcs(self):
        wcs = WCS(self.image.header)
        self.wcs = wcs
        return wcs
    
    #  @TODO: Move cutouts to new class.
    def cutout(self,y,x):
        self.cut = Cutout2D(
            self.data,
            position=self.wcs_init.wcs.crpix,
            wcs=self.wcs_init,
            size=(y,x)
        )
        
    def make_fits_from_cutout(self):
            self.cut_image = self.ImType(self.cut.data,self.cut.wcs.to_header())
        
    def find_missing_keys(self):
        self.missing_keys = set(self.header.keys())-set(self.cut_image.header.keys())
        
    def fill_cut_image_header(self):
        for key in self.missing_keys:
            self.cut_image.header[key] = self.header[key] 
           
    def make_cutout(self, y:int, x:int):
        """convenience function for making a cutout image and forward filling header"""
        self.cutout(y,x)
        self.make_fits_from_cutout()
        self.find_missing_keys()
        self.fill_cut_image_header()

class Plotter:
    
    def __init__(self, 
                 figsize:tuple = (8,8), 
                 zs: BaseInterval = ZScaleInterval(1000), 
                 subplots: int=1,
                 dpi: int=125
                ):
        self.figsize = figsize
        self.interval_finder = zs
        self.fig = plt.figure(figsize=self.figsize, dpi=dpi)
        self.subplots = subplots
        self.active_axis = 1
        self.axes = []
        
    def plot_img(self,image: Image, wcs:WCS = 'None'):
        if wcs=='None':
            wcs = image.wcs
        self.ax = self.fig.add_subplot(1,self.subplots,self.active_axis, projection=wcs)
        vmin,vmax = self.interval_finder.get_limits(image.data)
        self.ax.imshow(image.data,vmin=vmin,vmax=vmax)
        self.ax.grid()
        self.active_axis += 1
        self.axes.append(self.ax)

def main():
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

    args=parser.parse_args()
    
    curr = Path(args.file)
    savepath = Path.joinpath(Path(args.output_dir) /
                             'cutout',curr.stem+'_cutout'+curr.suffix).resolve()
    
    hdul = fits.open(curr)
    new_hdul = fits.HDUList()
    
    for i in hdul:
        if isinstance(i,fits.ImageHDU) or isinstance(i,fits.PrimaryHDU):
            img = Image(i)
            img.make_cutout(args.y,args.x)
            new_hdul.append(img.cut_image)
        else: #  Mostly expecting extension to be a binaryfitstable.  
            new_hdul.append(i)
    
    if args.plot:
        plotter = Plotter(subplots=2, dpi=200)
        plotter.plot_img(Image(new_hdul[0]))
        plotter.ax.set_title('cutout')
        plotter.plot_img(Image(hdul[0]))
        plotter.ax.set_title('uncut')
        plt.show(block=True)
    
    # Save
    new_hdul.writeto(str(savepath), output_verify='exception', 
                     overwrite=args.overwrite, checksum=False)
    print(f'saved to: {savepath}')
    
if __name__== '__main__':
    main()
    
   

    