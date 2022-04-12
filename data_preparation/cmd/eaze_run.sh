#!/bin/sh
#
# VERSION: 1997b 
# DATE:    Aug 1997
# 
#
# -------------------------------------------------------------------
# 
#         Eaze_Run
#         ==========
# 
#     This command file runs an eaze job.
# 
#   Turn on Verbose Mode
set -x
#
#  >>>>>>>>>>>>>>>>>>>> set  USERS STUFF <<<<<<<<<<<<<<<<<<<<<<<<<<<<<
# 
#... the root of your EAZE executeable
#
#MY_ROOT=$AFS_HOME/zeus/eaze_example
MY_ROOT=.
#
#MY_EXE =  $MY_ROOT/exe/$(ZARCH_TYPE)
MY_EXE=.
# 
#... The Program Name
# 
PROGRAM=jps_compton
#
LOGFILE=${PROGRAM}.log
ERRFILE=${PROGRAM}.err
#
EXE_TYPE=.exe
# 
PROGRAM_EXE=${PROGRAM}${EXE_TYPE}
#
# 
#... your  Programme run control cards
# 
\rm -f fort.7 
#
#ln -s  $MY_ROOT/extract_control.cards  fort.7
ln -s  compton_v1.cards  fort.7
#
#
#... your Input file 
#
\rm -f input 
#
# >>>>>>>>>>>>>>>>>>>>>>>> END USER STUFF <<<<<<<<<<<<<<<<<<<<<<<<<<
# 
#..  Magnetic field map
# 
\rm -f fort.22 
\rm -f for022.dat 
ln -s  $ZOF_ROOT/field/s3955c889y3000_sca.maps  fort.22
ln -s  $ZOF_ROOT/field/s3955c889y3000_sca.maps  for022.dat
# 
#
#  .. GAFs
# 
\rm -f gaf       
ln -s  $ZOF_ROOT/gaf  gaf
#
#
echo  "Running program $PROGRAM ........................."
#
time $MY_EXE/$PROGRAM_EXE > $LOGFILE 2> $ERRFILE
#
#
#   Check the status
#   ----------------
if [ $? != 0 ]
then
        echo "ERROR running executable for $PROGRAM."
        exit 1
fi
echo "Successfully ran executable for $PROGRAM."
echo ''
#
ls -l
#
echo ''
echo ' ------------------------------------------------------- '
echo ''
date
echo "$PROGRAM Done."
#
#   Turn off Verbose Mode
set -
#
