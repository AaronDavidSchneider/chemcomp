# chemcomp
Modeling the chemical composition of gas giants by accretion of pebbles, planetesimals and gas.

## install:
clone project:

`git clone git@github.com:chemcomppla/chemcomp.git`

move into directory:

`cd chemcomp`

Install conda environment:

`conda create -n chemcomp numpy "astropy>=4.0" scipy matplotlib pyyaml jupyter parse h5py pip pytables`

Install `chemcomp`:

`pip install -e .`

Adjust paths in `chemcomp/helper/main_helper` if you need different directories for output and config.

## Help:
Check the [WIKI](https://chemcomp.readthedocs.io/en/latest/ "wiki")
