#!/bin/bash

# Export path to rhocp-opt
# provide appropriate path here
export RHOCPOPT=~/projects/rhocp/rhocp-opt
# change mpirun to appropriate mpi command

--- sigmaV = 80 MPa ---
echo "Running sigmaV=80 MPa, stage 1..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=80 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.10 \
  Outputs/file_base=out_ICreep_473K_80MPa \
  > log80.run

echo "Running sigmaV=80 MPa, stage 2 (recover)..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=80 \
  Executioner/end_time=3e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_ICreep_473K_80MPa \
  --recover \
  > log80r1.run


--- sigmaV = 170 MPa ---
echo "Running sigmaV=170 MPa, stage 1..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=170 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.10 \
  Outputs/file_base=out_ICreep_473K_170MPa \
  > log170.run

echo "Running sigmaV=170 MPa, stage 2 (recover)..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=170 \
  Executioner/end_time=3e7 \
  Executioner/TimeStepper/log_dt=0.01 \
  Outputs/file_base=out_ICreep_473K_170MPa \
  --recover \
  > log170r2.run


# --- sigmaV = 350 MPa ---
echo "Running sigmaV=350 MPa, stage 1..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=350 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.020 \
  Outputs/file_base=out_ICreep_473K_350MPa \
  > log350.run

echo "Running sigmaV=350 MPa, stage 2 (recover)..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=350 \
  Executioner/end_time=3e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_473K_350MPa \
  --recover \
  > log350r1.run



