# `DIAGNOSIS` v0.1 (Low Latency HET Observation Planner for GW events) 

# Table of Contents

[Overview](https://github.com/Majoburo/GWHET-obs/edit/master/README.md#Overview)

[Running Diagnosis]

[Data Products]

[Examples]

Listener

Submiter
HETDEX observation tools for a GW event

## Overview
`DIAGNOSIS` is a low latency HET observation planner for GW events. `DIAGNOSIS` will:
- Continuosly listen to the [Gamma-ray Coordinates Network/Transient Astronomy Network (GCN/TAN)](https://gcn.gsfc.nasa.gov/) for alerts on gravitational wave events.
- When triggered, if the event is likely to be a binary neutron star merger (BNS) or a black hole-neutron star merger (BHNS), `DIAGNOSIS` will download the asociated skymap, identify *if* and *when* the 90% probability region falls within the HET pupil, and if so, inform the observers.
- `DIAGNOSIS` will also query the 2MASS catalog for galaxies within the observable 90% probability region, organize them by probability, and give their local sidereal times to start observations.
- Finally, if VIRUS where to locate the source of the graviational waves, `DIAGNOSIS` also provides a tool to report the coordinates back to the GCN/TAN.

This is a beta tool.

## Running `DIAGNOSIS`

### Requirements
`DIAGNOSIS` relies on a few third-party Python packages. These include:

- Astropy
- gcn
- healpy
- smtplib
- numpy
- matplotlib
- os
- scipy
- argparse
- astroquery
- pandas

The fastest way to install the dependencies is with pip, a package manager that comes with most Python distributions. To install these packages with pip, run the following command:
```
$ pip install astropy gcn healpy smtplib numpy matplotlib os scipy argparse astroquery pandas
```

### Download and Install
`DIAGNOSIS` is a library code that can be acquired simply with:
```
cd WHEREVER
git clone https://github.com/Majoburo/Diagnosis.git
```

### Running the code
### Data Products
### Examples
