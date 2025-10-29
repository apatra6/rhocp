#!/bin/bash

# Export path to rhocp-opt
# provide appropriate path here
export RHOCPOPT=~/projects/rhocp/rhocp-opt
# change mpirun to appropriate mpi command

# --- sigmaV = 125 MPa ---
echo "Running sigmaV=125 MPa, stage 1..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=125 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.10 \
  Outputs/file_base=out_ICreep_603K_125MPa \
  > log125.run

echo "Running sigmaV=125 MPa, stage 2 (recover)..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=125 \
  Executioner/end_time=3e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_ICreep_603K_125MPa \
  --recover \
  > log125r1.run


# # --- sigmaV = 250 MPa ---
echo "Running sigmaV=250 MPa, stage 1..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=250 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.10 \
  Outputs/file_base=out_ICreep_603K_250MPa \
  > log250.run

echo "Running sigmaV=250 MPa, stage 2 (recover)..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=250 \
  Executioner/end_time=3e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_ICreep_603K_250MPa \
  --recover \
  > log250r2.run


# --- sigmaV = 370 MPa ---
echo "Running sigmaV=370 MPa, stage 1..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=370 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.020 \
  Outputs/file_base=out_ICreep_603K_370MPa \
  > log370.run

echo "Running sigmaV=370 MPa, stage 2 (recover)..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=370 \
  Executioner/end_time=3e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_603K_370MPa \
  --recover \
  > log370r1.run


# --- sigmaV = 490 MPa ---
echo "Running sigmaV=490 MPa, stage 1..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=490 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.020 \
  Outputs/file_base=out_ICreep_603K_490MPa \
  > log490.run

echo "Running sigmaV=490 MPa, stage 2 (recover)..."
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=490 \
  Executioner/end_time=3e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_603K_490MPa \
  --recover \
  > log490r1.run

