# DuneWCana
Some scripts for analysing WC MC with Dune flux vectors

## DuneWCana.py
Loops through fiTQun / WCSim files, applies selection cuts and produces histograms.
If run without arguments it will loop through all files in all modes (FHC/RHC, nue/numu/nutau/nuebar/...)

Takes a single optional argument: a "header" python file that can be used to overwrite global variables in the script.

## Headers
This directory contains header files that can be used to run over only one horn mode / neutrino type combination.

## submitAllSBU.py
This script loops over all the header files in the Headers directory and submits a DuneWCana.py job for each of the headers

## setupSBU.sh
Sets up environment on the nngroup machines at SBU.
