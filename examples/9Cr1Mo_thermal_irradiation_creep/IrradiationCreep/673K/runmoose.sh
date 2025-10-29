#!/bin/bash

# Export path to rhocp-opt
# provide appropriate path here
export RHOCPOPT=~/projects/rhocp/rhocp-opt
# change mpirun to appropriate mpi command

# --- sigmaV = 200 MPa ---
echo "Running sigmaV=200 MPa, stage 1..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=200 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_673K_200MPa \
  > log200.run

echo "Running sigmaV=200 MPa, stage 2 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=200 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_673K_200MPa \
  --recover \
  > log200r1.run

echo "Running sigmaV=200 MPa, stage 3 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=200 \
  Executioner/end_time=1e8 \
  Executioner/TimeStepper/log_dt=0.002 \
  Outputs/file_base=out_673K_200MPa \
  --recover \
  > log200r2.run


echo "Running sigmaV=30 MPa, stage 1..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=30 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_673K_30MPa \
  > log30.run

echo "Running sigmaV=30 MPa, stage 2 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=30 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.01 \
  Outputs/file_base=out_673K_30MPa \
  --recover \
  > log30r1.run

echo "Running sigmaV=30 MPa, stage 3 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=30 \
  Executioner/end_time=1e8 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_673K_30MPa \
  --recover \
  > log30r2.run


echo "Running sigmaV=60 MPa, stage 1..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=60 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_673K_60MPa \
  > log60.run

echo "Running sigmaV=60 MPa, stage 2 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=60 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.01 \
  Outputs/file_base=out_673K_60MPa \
  --recover \
  > log60r1.run

echo "Running sigmaV=60 MPa, stage 3 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=60 \
  Executioner/end_time=1e8 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_673K_60MPa \
  --recover \
  > log60r2.run


echo "Running sigmaV=100 MPa, stage 1..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=100 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_673K_100MPa \
  > log100.run

echo "Running sigmaV=100 MPa, stage 2 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=100 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.01 \
  Outputs/file_base=out_673K_100MPa \
  --recover \
  > log100r1.run

echo "Running sigmaV=100 MPa, stage 3 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=100 \
  Executioner/end_time=1e8 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_673K_100MPa \
  --recover \
  > log100r2.run


echo "Running sigmaV=140 MPa, stage 1..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=140 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_673K_140MPa \
  > log140.run

echo "Running sigmaV=140 MPa, stage 2 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=140 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_673K_140MPa \
  --recover \
  > log140r1.run

echo "Running sigmaV=140 MPa, stage 3 (recover)..."
mpirun -n 20 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=140 \
  Executioner/end_time=1e8 \
  Executioner/TimeStepper/log_dt=0.002 \
  Outputs/file_base=out_673K_140MPa \
  --recover \
  > log140r2.run


