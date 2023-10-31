#!/usr/bin/env python

import os
import sys
import glob
import argparse
import numpy as np
import pandas as pd

def get_metadata(file):
    """Obtains info from an image reduced with the FLOWS pipeline.
    
    Parameters
    ----------
    file: str
        'photometry.ecsv' file in an image directory.
        
    Returns
    -------
    filt: str
        Filter used.
    mjd: float
        Time of observation.
    """
    filt = mjd = None
    
    with open(file, 'r') as info_file:
        for line in info_file.readlines():
            if 'photfilter' in line:
                filt = line.split()[-1].split('}')[0]
            if 'obstime-bmjd' in line:
                mjd = line.split()[-1].split('}')[0]

    return filt, mjd

def get_magnitude(file):
    """Extracts photometry from an image reduced with the FLOWS pipeline.
    
    Parameters
    ----------
    file: str
        'photometry.ecsv' file in an image directory.
        
    Returns
    -------
    mag: float
        SN magnitude.
    mag_err: float
        SN magnitude error.
    sub: int
        Whether the image was template subtracted:
        '-1' is template subtracted and '0' is no subtraction.
    """
    init_df = pd.read_csv(file, comment='#')
    df = init_df[init_df.starid==-1]
    sub = -1
    if len(df)==0:
        df = init_df[init_df.starid==0]
        sub = 0
    mag, mag_err = df.mag.values[0], df.mag_error.values[0]
    
    return mag, mag_err, sub

def flows_photometry(target, subtraction=1):
    """Extracts FLOWS photometry.
    
    Parameters
    ----------
    target: str
        Name of the target. Used for the output photometry file.
    subtraction: int, default '1'
        Whether to extract the photometry from template-subtracted
        images (-1), unsubtracted images (0) or both (1).
    """
    img_directories = glob.glob(f'{target}/*')

    phot_dir = {'mjd':[], 'mag':[], 'mag_err':[], 'filt':[], 'subtraction':[]}
    for directory in img_directories:
        if os.path.isdir(directory) is False:
            continue  # not a directory, so skip

        # get all the files in an image directory and get the photometry file
        dir_files = glob.glob(f'{directory}/*')
        try:
            phot_file = [file for file in dir_files if file.endswith('photometry.ecsv')][0]
        except:
            print(f'Skipping {directory} - no "photometry.ecsv" file')
            continue  # skip this image 

        filt, mjd = get_metadata(phot_file)
        mag, mag_err, sub = get_magnitude(phot_file)
        mag = np.round(mag, 2)
        mag_err = np.round(mag_err, 2)
        mjd = np.round(float(mjd), 2)

        phot_dir['mjd'].append(mjd)
        phot_dir['mag'].append(mag)
        phot_dir['mag_err'].append(mag_err)
        phot_dir['filt'].append(filt)
        phot_dir['subtraction'].append(sub)
        
    
    phot_df = pd.DataFrame(phot_dir)
    
    # sort by filter, then by mjd
    sorter_dict = {'B':0, 'gp':1, 'V':2, 'rp':3, 'ip':4}
    phot_df['filt_num'] = [sorter_dict[filt] for filt in phot_df.filt.values]
    phot_df.sort_values(['filt_num', 'mjd'], inplace=True, ignore_index=True)
    
    # mask by subtraction
    if subtraction in [-1, 0]:
        phot_df = phot_df[phot_df.subtraction==subtraction]
    
    phot_df.to_csv(f'{target}_phot.csv', index=False)


def create_snoopy_file(target, z, ra, dec):
    """Creates a snpy file from FLOWS photometry

    Parameters
    ----------
    target : str
        Name of the target. Used for the output photometry file.
    z : float
        Target's redshift.
    ra : float
        Target's right ascension.
    dec : float
        Target's declination.
    """
    phot_df = pd.read_csv(f'{target}_phot.csv')

    with open(f'{target}_snpy.dat', 'w') as outfile:
        # metadata
        outfile.write(f'{target} {z} {ra} {dec}\n')

        for filt in phot_df.filt.unique():
            # add each filter
            outfile.write(f'filter {filt}\n')

            filt_df = phot_df[phot_df.filt==filt]
            for t, m, me in zip(filt_df.mjd.values, 
                                filt_df.mag.values, 
                                filt_df.mag_err.values):
                outfile.write(f'{t}\t{m}\t{me}\n')
            

def main(args=None):
    description = f"FLOWS photometry"
    usage = "flows_photometry <directory> [options]"
    
    if not args:
        args = sys.argv[1:] if sys.argv[1:] else ["--help"]
        
    parser = argparse.ArgumentParser(prog='flows_photometry',
                                     usage=usage,
                                     description=description
                                     )
    parser.add_argument("directory",
                        help="name of the directory"
                        )
    parser.add_argument("-s",
                        "--subtraction",
                        dest="subtraction",
                        action="store",
                        default=1,
                        choices=[-1, 0, 1],
                        type=int,
                        help=("Whether to extract the photometry from template-subtracted "
                              "images (-1), unsubtracted images (0) or both (1; deafault).")
                        )
    parser.add_argument("--snpy",
                        dest="snpy",
                        nargs='+',
                        type=str,
                        required=False,
                        help=("Creates a snpy file. Redshift, right ascension and declination "
                              "must be given (in that order), e.g. '--snpy 0.0532 12.340 0.344'.")
                        )
    
    args = parser.parse_args(args)
    flows_photometry(args.directory, args.subtraction)

    if args.snpy is not None:
        z, ra, dec = args.snpy
        create_snoopy_file(args.directory, z, ra, dec)
    
if __name__ == "__main__":
    main(sys.argv[1:])