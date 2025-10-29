#!/bin/bash

# Define Simulation Parameters for Thermal Creep Simulation

NPROCS=24
STRESS=(121 87 52)
TEMP=823

# Export path to rhocp-opt
# provide appropriate path here
export RHOCPOPT=~/projects/rhocp/rhocp-opt
# change mpirun to appropriate mpi command

for stress in "${STRESS[@]}"
do
    # Run the simulation with the specified parameters
    mpirun -n "$NPROCS" $RHOCPOPT -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e5 \
        Executioner/TimeStepper/log_dt=0.10 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" > "log${stress}.run"

    # Run the simulation with a longer end time and smaller log_dt
    mpirun -n "$NPROCS" $RHOCPOPT -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e6 \
        Executioner/TimeStepper/log_dt=0.05 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r1.run"

        # Run the simulation with a longer end time and smaller log_dt
    mpirun -n "$NPROCS" $RHOCPOPT -i bcc_pxtal.i \
        sigmaV="$stress" \
        Materials/CPStressUpdate/temp="$TEMP" \
        Materials/CPStressUpdate/climbmodel=true \
        Executioner/end_time=5e7 \
        Executioner/TimeStepper/log_dt=0.005 \
        Outputs/file_base="out_Creep_${TEMP}K_${stress}MPa" --recover > "log${stress}r2.run"
done

