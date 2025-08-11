# Toricity of vertically parametrized systems with applications to reaction network theory

This repoistory contains files for the forthcoming paper _Toricity of vertically parametrized systems with applications to reaction network theory_ by Elisenda Feliu and Oskar Henriksson.

## File descriptions
The repository contains the following files:
* A directory `src` that contains Julia functions for testing whether a network satisfies the various notions of toricity treated in the paper, as well as properties that can easilly be checked in the presence of toricity such as (local) ACR and multistationarity. These functions also come with tests in a directory `test`.
* A directory `results` that contains list of networks in [ODEbase](https://www.odebase.org/) (as of November 2, 2023) with various properties. All the results are summerized in a text file `report.txt` and the table `results.csv` (see the file `legend.md` for explanation of the column headings).

## Examples
See the notebook `IDH_example.ipynb`.

## Dependencies
The code is based on `Oscar v1.1.1` and `DataStructures v0.18.20`, `Graphs v1.12.0` and `HomotopyContinuation v2.15.0`. 

For exact dependencies, see the file `Manifest.toml`.

