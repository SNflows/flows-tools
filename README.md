[![CI Tests](https://github.com/SNflows/flows-tools/actions/workflows/python-package.yml/badge.svg)](https://github.com/SNflows/flows-tools/actions/workflows/python-package.yml)
[![Test coverage](https://img.shields.io/badge/Tests%20Coverage-73%25-yellow.svg)](https://github.com/SNflows/flows-tools/pull/27#issuecomment-1317015483)

# flows-tools   
Observer tools for the FLOWS pipeline. Contains python packages, scripts, and notebooks for the following tasks:
 - Pointing, OB generation helper with VLT HAWK-I, finder chart maker.
 - Image cutouts
 - Template subtraction
 - Fixing fits image headers of multi-extension fits files
 - Hubble Diagrams
 - Lightcurve fitting of Type Ia SN lightcurves

---
## Available Packages:
 - **flows_get_brightest**: The VLT HAWK-I pointing script that will create a finding chart for HAWK-I 
H-band images given a flows target, and optional rotation (position angle) and telescope offsets.
Uses `tendrils` to query objects. Optionally plots a finding chart. 
Easily extensible to use with other instruments.

### install with:
``pip install flows_get_brightest``   
**NOTE**: Installation requires `python >= 3.10` due to `tendrils` dependency.  
### Usage:
Installs a user script `get_brightest` that can be run from the command line.
``get_brightest --help``.

### test with: 
``pip install flows_get_brightest[test]``    
``pytest -v --cov``
``mypy flows_get_brightest/ typings/tendrils/ --config-file=pyproject.toml --check-untyped-defs``
``flake8 . --count --exit-zero --max-complexity=10 --max-line-length=120 --statistics --extend-ignore=F821 --exclude .git,.idea,.mypy_cache,Notebooks/``
``black --check --verbose flows_get_brightest/``
Note: user script only tested to work on Linux and Mac OS.

---
## Available scripts:
 - fixlcogt.py A script for fixing LCOGT headers not preserved by other scripts such as imagematch.
 - cutouts.py A flexible script for making image cutouts while preserving original headers and extension and applying 
the same cutout to all image extensions present in the fits file.

## Available notebooks:
 - Subtraction.ipynb Example notebook for doing manual subtractions following @sholmbo's method.
 - image_cutout.ipynb Example notebook that was the inspiration for cutouts.py.
 - hubble diagram notebook.
 - lightcurve fitting notebook.
 
