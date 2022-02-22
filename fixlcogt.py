import argparse
import astropy.io.fits as fits
from astropy.wcs import WCS
import numpy as np
from pathlib import Path


def read_hdul(filepath, header_index=0):
    hdul = fits.open(filepath)
    wcs = WCS(hdul[header_index].header)
    return hdul, wcs


def swap_wcs(hdul_orig, hdul_diff):
    new_hdul = fits.HDUList()

    new_hdul.append(hdul_diff[0])
    new_hdul.append(hdul_orig[1])
    new_hdul.append(hdul_orig[2])
    return new_hdul


def get_basename_with_suffix(filename: Path, suffix=''):
    return Path(
        str(filename)[: str(filename).rfind(''.join(filename.suffixes))]
        + suffix
        + ''.join(filename.suffixes)
    )


def parse():
    parser = argparse.ArgumentParser(prog='FixTemplateSubHeaderforLCOGT', description='Attach header and extra dimensions from second \
                                     file to the first, while checking to make sure WCS is identical.\
                                     If not, uses WCS offirst file.')

    parser.add_argument("diff_file", help="Specify the diff file here")
    parser.add_argument("orig_file", help="Specify the original file here")
    parser.add_argument("--suffix", help='output suffix', type=str, default="_hfix")
    parser.add_argument("-o", "--overwrite", help='overwrite output if exists', action='store_true', default=False)

    return parser.parse_args()


def main():
    args = parse()

    diff_path = Path(args.diff_file)
    orig_path = Path(args.orig_file)

    # Get HDU list and WCS
    hdul_diff, wcs_diff = read_hdul(diff_path)
    hdul_orig, wcs_orig = read_hdul(orig_path)

    if wcs_diff.low_level_wcs == wcs_diff.low_level_wcs:
        save_file = get_basename_with_suffix(diff_path, args.suffix)
        new_hdul = swap_wcs(hdul_orig, hdul_diff)
        new_hdul.writeto(save_file, overwrite=args.overwrite)
        print(f'saved to: {save_file}')


if __name__ == '__main__':
    main()
