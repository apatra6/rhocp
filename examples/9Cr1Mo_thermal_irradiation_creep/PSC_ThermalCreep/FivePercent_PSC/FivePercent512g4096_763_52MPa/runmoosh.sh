LoadRemovalTime=446
LoadRemovalTimeMinus=$(echo "$LoadRemovalTime - 0.01" | bc)
LoadRemovalTimePlus=$(echo "$LoadRemovalTime + 0.5" | bc) 
HeatupTime=3500


mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i   Executioner/end_time=445

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
    Executioner/end_time=$LoadRemovalTimePlus \
    Materials/CPStressUpdate/tol=1e-3 --recover

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
    Executioner/end_time=$HeatupTime \
    Materials/CPStressUpdate/tol=5e-7 --recover    

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
    Executioner/end_time=1e8 \
    Materials/CPStressUpdate/tol=5e-7 --recover        
