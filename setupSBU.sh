#!/bin/bash

source /home/cvilela/nuPRISM/Software/Analysis/Source_At_Start_nuPRISM.sh
export PYTHONPATH=${ROOTSYS}/lib:${PYTHONPATH}
export PYTHONPATH=/storage/shared/cvilela/DuneWC/DuneWCana/guppy-0.1.10/build/lib.linux-x86_64-2.7:${PYTHONPATH}

unset WCSIMDIR
unset G4WORKDIR

export WCSIMDIR=/storage/shared/cvilela/WCSim_v151/WCSim
export G4WORKDIR=${WCSIMDIR}/bin/
