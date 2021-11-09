# flows-tools# flows-tools
Observer tools for the FLOWS pipeline

Available scripts:
 - get_brightest.py The HAWK-I pointing script that will create a finding chart for HAWK-I H-band images given a flows target, and optional rotation (position angle) and telescope offsets.
 - fixlcogt.py A script for fixing LCOGT headers not preserved by other sripts such as imagematch.
 - cutouts.py A flexible script for making image cutouts while preserving original headers and extension and applying the same cutout to all image extensions present in the fits file.

Available notebooks:
 - Subtraction.ipynb Example notebook for doing manual subtractions following @sholmbo's method.
 - image_cutout.ipynb Example notebook that was the inspiration for cutouts.py.
