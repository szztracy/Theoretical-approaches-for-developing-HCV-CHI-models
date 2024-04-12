#! /usr/bin/env bash

set -eu

if [ "$#" -ne 2 ]; then
  script_name=$(basename $0)
  echo "Usage: ${script_name} exp_id cfg_file"
  exit 1
fi

export TURBINE_LOG=0 TURBINE_DEBUG=0 ADLB_DEBUG=0
# export TURBINE_STDOUT=out-%%r.txt
export TURBINE_STDOUT=
export ADLB_TRACE=0
export EMEWS_PROJECT_ROOT=$( cd $( dirname $0 )/.. ; /bin/pwd )
# source some utility functions used by EMEWS in this script                                                                                 
source "${EMEWS_PROJECT_ROOT}/etc/emews_utils.sh"

export EXPID=$1
export TURBINE_OUTPUT=$EMEWS_PROJECT_ROOT/experiments/$EXPID
check_directory_exists

CFG_FILE=$2
source $CFG_FILE

echo "--------------------------"
echo "WALLTIME:              $CFG_WALLTIME"
echo "PROCS:                 $CFG_PROCS"
echo "PPN:                   $CFG_PPN"
echo "TRIALS:                $CFG_TRIALS"
echo "MODEL:                 $CFG_MODEL_X"
echo "VLF:                   $CFG_VIRAL_LOAD_FILE"
echo "--------------------------"

export PROCS=$CFG_PROCS
export QUEUE=$CFG_QUEUE
export WALLTIME=$CFG_WALLTIME
export PPN=$CFG_PPN
export TURBINE_JOBNAME="${EXPID}_job"
export PROJECT=$CFG_PROJECT

# if R cannot be found, then these will need to be
# uncommented and set correctly.
# export R_HOME=/path/to/R
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$R_HOME/lib



# Resident task workers and ranks
export TURBINE_RESIDENT_WORK_WORKERS=1
export RESIDENT_WORK_RANKS=$(( PROCS - 2 ))

# EQ/R location
# EQR=/lcrc/project/EMEWS/bebop/repos/spack/opt/spack/linux-centos7-broadwell/gcc-7.1.0/eqr-1.0-5hb4aszbbtezlifks6fz4g24zldnkdbx
EQPY=$EMEWS_PROJECT_ROOT/ext/EQ-Py
export PYTHONPATH=$EMEWS_PROJECT_ROOT/python:$EQPY

export SITE=bebop


# set machine to your schedule type (e.g. pbs, slurm, cobalt etc.),
# or empty for an immediate non-queued unscheduled run
MACHINE="slurm"

if [ -n "$MACHINE" ]; then
  MACHINE="-m $MACHINE"
fi

mkdir -p $TURBINE_OUTPUT/instances
mkdir -p $TURBINE_OUTPUT/results
cp $CFG_FILE $TURBINE_OUTPUT/cfg.cfg

DEAP_CFG=$TURBINE_OUTPUT/deap.yaml
cp $EMEWS_PROJECT_ROOT/data/$CFG_DEAP_CFG_FILE $DEAP_CFG

VL_FILE=$TURBINE_OUTPUT/viral_load.csv
cp $EMEWS_PROJECT_ROOT/data/$CFG_VIRAL_LOAD_FILE $VL_FILE

cp $EMEWS_PROJECT_ROOT/data/$CFG_GA_PARAMS_FILE $TURBINE_OUTPUT/ga_params.json

CMD_LINE_ARGS="$* -trials=$CFG_TRIALS -model_x=$CFG_MODEL_X -viral_load_file=$VL_FILE "
CMD_LINE_ARGS+=" -deap_cfg=$DEAP_CFG"

# Add any script variables that you want to log as
# part of the experiment meta data to the USER_VARS array,
# for example, USER_VARS=("VAR_1" "VAR_2")
USER_VARS=("MODEL_DIR" "STOP_AT" "MODEL_PROPS" \
 "STOP_AT")

# Enable MPI multi-threading
export TURBINE_MPI_THREAD=1 

# Put RInside on the LD_LIBRARY_PATH
export RINSIDE_LIB=/lcrc/project/EMEWS/bebop/sfw/R-4.2.2/lib64/R/library/RInside/lib
export LD_LIBRARY_PATH=$RINSIDE_LIB:$LD_LIBRARY_PATH

# Use fast interconnect MPI fabric
export I_MPI_FABRICS=shm:ofi

# Launch with srun rather than mpirun
export TURBINE_LAUNCHER=srun

# log variables and script to to TURBINE_OUTPUT directory
log_script

# echo's anything following this standard out
# set -x

swift-t -n $PROCS $MACHINE -p \
    -r $EQPY -I $EQPY \
    -e EMEWS_PROJECT_ROOT \
    -e TURBINE_OUTPUT \
    -e TURBINE_LOG \
    -e TURBINE_DEBUG \
    -e ADLB_DEBUG \
    -e SITE \
    $EMEWS_PROJECT_ROOT/swift/ga_workflow.swift $CMD_LINE_ARGS

chmod g+rw $TURBINE_OUTPUT/*.tic
