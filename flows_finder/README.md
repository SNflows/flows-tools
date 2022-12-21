https://pypi.org/project/flows_get_brightest
# Finder chart maker and observational planner
[![Python package](https://github.com/SNflows/flows-tools/actions/workflows/python-package.yml/badge.svg)](https://github.com/SNflows/flows-tools/actions/workflows/python-package.yml)

### install with:
``pip install flows_get_brightest``   
**NOTE**: Installation requires `python >= 3.10` due to `tendrils` dependency.  
### Usage:
Installs a user script `get_brightest` that can be run from the command line.
``get_brightest --help``.

### test with: 
``pip install flows_get_brightest[test]``    
``pytest``

## Description
**Note:While currently only vlt HAWK-I is implemented, the code is abstracted to work
generically with any Instrument. See "How to extend?" below.**

The pointing script that will create a finding chart for (currently) HAWK-I 
H-band images given a flows target, and optional rotation (position angle) and telescope offsets.
Uses `tendrils` to query objects. Optionally plots a finding chart. 
Easily extensible to use with other instruments.

## How to extend?
Simply add `Instrument` classes for new telescope and instrument combos in `instruments.py`.

`Plotter` module/class "Takes an observer with `WCS`, `Target`, `Plan`, `Instrument`, `Regions`, and `Corners` and makes a finding chart."
`Regions` and `Target` are queried from transient name via `tendrils`, and `Plan` is specified based on offsets, rotations, time, etc.

Meaning that new instrument needs to define its observing area (which will be converted `Corners` class instances) and used for
making finder charts, or calculating the bright star(s) in a field.

The code could be extended to automatically schedule the observation by intelligently
avoiding or including stars based on catalog values, which will be needed in a few years
time when FLOWS project needs to grow 10x in size. The API for this is already there, 
and the only missing piece is to come up with scheduling criteria and write an algorithm 
for using the exposed interfaces of the `Observer`, `Plan`, and `Target` classes.
