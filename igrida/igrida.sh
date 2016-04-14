#!/usr/bin/env bash
#OAR -l {mem_core < 2048}/nodes=1/core=3,walltime=168:00:00
#OAR -l {mem_core > 2048}/nodes=1/core=1,walltime=168:00:00
#OAR --array-param-file /udd/rtavenar/src/LR-DTW/igrida/params.data
#OAR -O /temp_dd/igrida-fs1/rtavenar/LR-DTW/output/%jobid%.output
#OAR -E /temp_dd/igrida-fs1/rtavenar/LR-DTW/output/%jobid%.error
# set -xv

SOURCEDIR=/udd/rtavenar/src/LR-DTW/c/ && WORKINGDIR=${SOURCEDIR} && EXECUTABLE=${SOURCEDIR}/test_1nn
cd ${WORKINGDIR} && ${EXECUTABLE} $*


