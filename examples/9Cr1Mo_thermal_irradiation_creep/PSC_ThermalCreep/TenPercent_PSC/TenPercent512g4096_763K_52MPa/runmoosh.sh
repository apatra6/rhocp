LoadRemovalTime=900
LoadRemovalTimeMinus=$(echo "$LoadRemovalTime - 0.01" | bc)
LoadRemovalTimePlus=$(echo "$LoadRemovalTime + 0.5" | bc) 
HeatTime=3500



mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i   Executioner/end_time=900 --recover

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
    Executioner/end_time=$LoadRemovalTimePlus \
    Materials/CPStressUpdate/tol=1e0 --recover

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
    Executioner/end_time=$HeatTime \
    Materials/CPStressUpdate/tol=5e-7 --recover    

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
    Executioner/end_time=1e8 \
    Materials/CPStressUpdate/tol=5e-7 --recover    
