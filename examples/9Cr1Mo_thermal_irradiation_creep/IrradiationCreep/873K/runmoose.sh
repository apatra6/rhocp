# Export path to rhocp-opt
# provide appropriate path here
export RHOCPOPT=~/projects/rhocp/rhocp-opt
# change mpirun to appropriate mpi command

# Simulation for sigmaV = 13 MPa
mpirun -n 20 $RHOCPOPT -i bcc_pxtal.i \
  sigmaV=13 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.1 \
  Outputs/file_base=out_ICreep_873K_13MPa \
  > log13.run

mpirun -n 20 $RHOCPOPT -i bcc_pxtal.i \
  sigmaV=13 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.008 \
  Outputs/file_base=out_ICreep_873K_13MPa \
  --recover \
  > log13r1.run

mpirun -n 20 $RHOCPOPT\
  -i bcc_pxtal.i \
  sigmaV=13 \
  Executioner/end_time=5e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_873K_13MPa \
  --recover \
  > log13r2.run


# # Simulation for sigmaV = 9 MPa
mpirun -n 20 $RHOCPOPT\
  -i bcc_pxtal.i \
  sigmaV=9 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.1 \
  Outputs/file_base=out_ICreep_873K_9MPa \
  > log9.run

mpirun -n 20 $RHOCPOPT\
  -i bcc_pxtal.i \
  sigmaV=9 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_ICreep_873K_9MPa \
  --recover \
  > log9r1.run

mpirun -n 20 $RHOCPOPT\
  -i bcc_pxtal.i \
  sigmaV=9 \
  Executioner/end_time=5e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_873K_9MPa \
  --recover \
  > log9r2.run


# # Simulation for sigmaV = 4 MPa
mpirun -n 20 $RHOCPOPT\
  -i bcc_pxtal.i \
  sigmaV=4 \
  Executioner/end_time=5e5 \
  Executioner/TimeStepper/log_dt=0.1 \
  Outputs/file_base=out_ICreep_873K_4MPa \
  > log4.run

mpirun -n 20 $RHOCPOPT\
  -i bcc_pxtal.i \
  sigmaV=4 \
  Executioner/end_time=1e7 \
  Executioner/TimeStepper/log_dt=0.05 \
  Outputs/file_base=out_ICreep_873K_4MPa \
  --recover \
  > log4r1.run

mpirun -n 20 $RHOCPOPT\
  -i bcc_pxtal.i \
  sigmaV=4 \
  Executioner/end_time=5e7 \
  Executioner/TimeStepper/log_dt=0.005 \
  Outputs/file_base=out_ICreep_873K_4MPa \
  --recover \
  > log4r2.run
