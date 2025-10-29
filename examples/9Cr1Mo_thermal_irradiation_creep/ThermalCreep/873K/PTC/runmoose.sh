#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=20
STRESS=(9 26 39 53 66)  # Stress values in MPa
TIME1=(5e5 5e5 5e5 5e5 5e5)     # End times for initial simulations
TIME2=(7e7 7e7 7e7 7e7 7e7)     # End times for recovered simulations
TEMP=873

# Export path to rhocp-opt
# provide appropriate path here
export RHOCPOPT=~/projects/rhocp/rhocp-opt
# change mpirun to appropriate mpi command

# Iterate over the array index to pair STRESS, TIME1, and TIME2 correctly
for i in "${!STRESS[@]}"
do
    stress=${STRESS[$i]}
    time1=${TIME1[$i]}
    time2=${TIME2[$i]}

    echo "Running initial simulation for STRESS=${stress} MPa, TEMP=${TEMP} K, end_time=${time1}"

    # Run the initial simulation
    mpirun -n "$NPROCS" $RHOCPOPT -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time1" \
        Executioner/TimeStepper/log_dt=0.10 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    echo "Initial run complete. Starting recovered simulation with end_time=${time2}"

    # Run the recovered simulation
    mpirun -n "$NPROCS" $RHOCPOPT -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time2" \
        Executioner/TimeStepper/log_dt=0.05 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"
done
