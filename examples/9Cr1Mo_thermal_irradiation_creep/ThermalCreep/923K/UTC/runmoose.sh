#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=20
STRESS=(55 65 75 88 100 110)
TIME1=(5e5 5e5 5e5 1e5 5e4 5e4)     # End times for initial simulations
TIME2=(7e7 3e7 3e7 5e6 3e6 3e6)     # End times for recovered simulations
TEMP=923

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


STRESS=(125 150 160 200)
TIME1=(5e4 1e3 1e3 1e2)     # End times for initial simulations
TIME2=(1e6 5e4 3e4 3e3)     # End times for recovered simulations
TEMP=923

Iterate over the array index to pair STRESS, TIME1, and TIME2 correctly
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
        Executioner/TimeStepper/log_dt=0.050 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    echo "Initial run complete. Starting recovered simulation with end_time=${time2}"

    # Run the recovered simulation
    mpirun -n "$NPROCS" $RHOCPOPT -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time2" \
        Executioner/TimeStepper/log_dt=0.01 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"
done


STRESS=(120 172) 
TIME1=(5e4 5e2)     # End times for initial simulations
TIME2=(1e6 2e4)     # End times for recovered simulations
TEMP=923

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
        Executioner/TimeStepper/log_dt=0.050 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    echo "Initial run complete. Starting recovered simulation with end_time=${time2}"

    # Run the recovered simulation
    mpirun -n "$NPROCS" $RHOCPOPT -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time="$time2" \
        Executioner/TimeStepper/log_dt=0.01 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"
done
