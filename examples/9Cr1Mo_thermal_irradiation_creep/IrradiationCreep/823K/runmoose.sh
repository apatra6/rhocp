# Export path to rhocp-opt
# provide appropriate path here
export RHOCPOPT=~/projects/rhocp/rhocp-opt
# change mpirun to appropriate mpi command

# --- sigmaV = 13 MPa ---
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=13 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.1 \
  Outputs/file_base=out_ICreep_823K_13MPa \
  > log13.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=13 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_ICreep_823K_13MPa \
  --recover \
  > log13r1.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=13 \
  Executioner/end_time=5e7 \
  Executioner/TimeStepper/log_dt=0.01 \
  Outputs/file_base=out_ICreep_823K_13MPa \
  --recover \
  > log13r2.run


# --- sigmaV = 26 MPa ---
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=26 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_ICreep_823K_26MPa \
  > log26.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=26 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.01 \
  Outputs/file_base=out_ICreep_823K_26MPa \
  --recover \
  > log26r1.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=26 \
  Executioner/end_time=5e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_823K_26MPa \
  --recover \
  > log26r2.run


# --- sigmaV = 52 MPa ---
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=52 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.050 \
  Outputs/file_base=out_ICreep_823K_52MPa \
  > log52.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=52 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_823K_52MPa \
  --recover \
  > log52r1.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=52 \
  Executioner/end_time=5e7 \
  Executioner/TimeStepper/log_dt=0.001 \
  Outputs/file_base=out_ICreep_823K_52MPa \
  --recover \
  > log52r2.run


# --- sigmaV = 87 MPa ---
mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=87 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.050 \
  Outputs/file_base=out_ICreep_823K_87MPa \
  > log87.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=87 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.001 \
  Outputs/file_base=out_ICreep_823K_87MPa \
  --recover \
  > log87r1.run

mpirun -n 24 $RHOCPOPT \
  -i bcc_pxtal.i \
  sigmaV=87 \
  Executioner/end_time=5e7 \
  Executioner/TimeStepper/log_dt=0.001 \
  Outputs/file_base=out_ICreep_823K_87MPa \
  --recover \
  > log87r2.run
