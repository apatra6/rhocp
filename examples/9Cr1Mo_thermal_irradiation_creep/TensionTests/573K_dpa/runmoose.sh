mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props_0p6.in \
Outputs/file_base=out_573K_0p6dpa

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props_1p5.in \
Outputs/file_base=out_573K_1p5dpa 

mpirun -n 24 ~/projects/rhocp/rhocp-opt -i bcc_pxtal.i \
Materials/CPStressUpdate/propsFile=bcc_props_0p06.in \
Outputs/file_base=out_573K_0p06dpa
