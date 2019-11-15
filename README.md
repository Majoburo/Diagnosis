# `DIAGNOSIS` v0.1 (Low Latency HET Observation Planner for GW events) 

# Table of Contents

[Overview](https://github.com/Majoburo/GWHET-obs/edit/master/README.md#Overview)

[Running `DIAGNOSIS`](https://github.com/Majoburo/GWHET-obs/edit/master/README.md#Running-`DIAGNOSIS`)

[Data Products](https://github.com/Majoburo/GWHET-obs/edit/master/README.md#Data-Products)


## Overview
`DIAGNOSIS` is a low latency HET observation planner for GW events. `DIAGNOSIS` is a modular code that will:
- Continuosly listen to the [Gamma-ray Coordinates Network/Transient Astronomy Network (GCN/TAN)](https://gcn.gsfc.nasa.gov/) for alerts on gravitational wave events.
- When triggered, if the event is likely to be a binary neutron star merger (BNS) or a black hole-neutron star merger (BHNS), `DIAGNOSIS` will download the asociated skymap, identify *if* and *when* the 90% confidence region falls within the HET pupil, and if so, inform the observers.
- `DIAGNOSIS` will also query a galaxy catalog (for now GLADE) for galaxies within the observable 90% probability region, organize them by probability, give their local sidereal times (LSTs) to start observations and make a .tsl file for phase II in HET observations.
- Finally, if VIRUS were to locate the source of the gravitational waves, `DIAGNOSIS` also provides a tool to report the coordinates back to the GCN/TAN.

This is a beta tool.

## Running `DIAGNOSIS`

### Requirements

`DIAGNOSIS` relies on a few third-party Python packages. These include:

- astropy
- pygcn
- healpy
- numpy
- matplotlib
- scipy
- argparse
- astroquery
- pandas

The fastest way to install the dependencies is with pip, a package manager that comes with most Python distributions. To install these packages with pip, run the following command:
```
$ pip install astropy pygcn healpy numpy matplotlib scipy argparse astroquery pandas
```

### Download and Install
`DIAGNOSIS` can be acquired simply with:
```
cd WHEREVER
git clone https://github.com/Majoburo/Diagnosis.git
```

In order to alert the relevant parties, `DIAGNOSIS` requires sensitive information that cannot be saved in a public git repository. Before running, be sure to acquire `recipients.py` from:
```
cd WHEREVER
scp username@stampede2.tacc.utexas.edu:/work/03237/majoburo/stampede2/DIAGNOSIS/recipients.py
```

### Running the code

`DIAGNOSIS` is meant to be constantly running in the background in order to listen and process any GCN/TAN alerts.
`DIAGNOSIS` can also be run once to process a specific event within the gracedb database using the toggle `--graceid`.
Details respect to the inputs can be seen in the arguments below.
```
usage: diagnose.py [-h] [-g GRACEID] recipients

Receive and parse GCN alerts, alert observers and create observing tools.

positional arguments:
  recipients            Specify python file with list of recipients.

optional arguments:
  -h, --help            show this help message and exit
  -g GRACEID, --graceid GRACEID
                        graceID of event
```
If VIRUS were to locate the source of the gravitational waves, `DIAGNOSIS` also provides submit_gracedb.py, a tool to report the coordinates back to the GCN/TAN.

In order to submit GCN alerts you will have to sign up to the GCN network. Please follow steps 1 and 2 in the [LIGO-Virgo EM Follow-Up Tutorial](https://dcc.ligo.org/public/0118/G1500442/010/ligo-virgo-emfollowup-tutorial.html) to do so.

Details with respect to the input arguments of submit_gracedb.py are below.
```
usage: submit_gracedb.py [-h] [-c COMMENT]
                         group grace_id ralist declist rawithlist decwithlist
                         durationlist starttimelist

Submit event to GraceDB

positional arguments:
  group          MOU group responsible
  grace_id       Identifier for the GW event
  ralist         List of ra of centers of footprints (degrees)
  declist        List of ra of centers of footprints (degrees)
  rawithlist     List (or one for all) of footprint widths in ra (degrees)
  decwithlist    List (or one for all) of footprint widths in dec (degrees)
  durationlist   List (or one for all) of exposure times in sec
  starttimelist  List of beginnings of exposures (utc)

optional arguments:
  -h, --help     show this help message and exit
  -c COMMENT     Comments
```

### Data Products

#### {EVENT}.fits.gz

Healpy probability skymap of the GW event.

#### galaxies_{CATALOG}_{EVENT}.dat

List of all galaxies, within the 90% confidence region observable by HET, with at most 99% of the probability of the most probable galaxy .

#### LSTs_{EVENT}.out

List of at most 100 most probable observable galaxies with their corresponding LSTs.

#### {EVENT}.tsl

List of at most 100 most probable observable galaxies formated for phase II observations.

#### MOLL_GWHET_{EVENT}.pdf
![MOLL_GWHET_{EVENT}.pdf](https://github.com/Majoburo/Diagnosis/blob/master/MOLL_GWHET_MS191113v.png)

Mollweide projection plot of the sky for the time of the alert. It displays:
- PURPLE: HET pupil.
- GREEN: 18 degree circle around the sun.
- GREY: AirMass > 2.5.
- BLUE: 90% probability region.
- PINK: Bellow the horizon.

#### LSTs_{EVENT}.pdf
![LSTs_{EVENT}.pdf](https://github.com/Majoburo/Diagnosis/blob/master/LSTs_MS191113v.pdf)


Plot of at most 100 most probable observable galaxies with their corresponding LSTs.
