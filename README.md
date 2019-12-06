# SKReact
A reactor neutrino simulator for Super-Kamiokande, simulating production in
reactors, oscillation to SK, cross section for IBD on free protons, and detector
smearing. Built referencing M. Baldoncini *et al.*, *"Reference worldwide model
for antineutrinos from reactors"*,
[arXiv:1411.6475v2](https://arxiv.org/abs/1411.6475).


## Installation and Running
SKReact uses python 3(.7.0), it does not support python 2.
Install the required modules using pip: 

`$ pip install -r requirements.txt`.

If it is your first time running and you don't have a `reactors_main.pkl` file,
you can generate your own by downloading the `.xls` files containing reactor
information from the PRIS kindly compiled by INFN
[here](https://www.fe.infn.it/radioactivity/antineutrino/index.html#download).

Download the "Input database" files (you will need to request the download), and
place the files into `react_p/` in the SKReact directory, with their names unchanged.
Simply create the directory if it does not already exist, or change the
`REACT_DIR` variable in `params.py` if you want to pull the data from somewhere
else.

The first time running without a `reactors_main.pkl` file will be very slow to
start as the relevant information needs to be pulled and compiled from the
`.xls` files. SKReact tries to handle all inconcistencies across the files and
will print any changes in reactor information as it compiles, if you want
verbose import/import errors, change the `VERBOSE_IMPORT` and
`VERBOSE_IMPORT_ERR` values in `params.py`.

Running from then on with a `reactors_main.pkl` file will be much faster on
startup. There is a chance you will need to recompile from `.xls` between
SKReact releases.


## Features

### Reactor List