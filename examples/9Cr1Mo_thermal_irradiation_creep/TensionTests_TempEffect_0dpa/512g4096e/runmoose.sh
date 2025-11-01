NPROC=20
RHOCPOPT=~/projects/rhocp-vr-251014/rhocp-opt

# Strain rate: 3e-3 /s
for temp in {723..923..50}; do
    mpirun -n "$NPROC" $RHOCPOPT -i bcc_pxtal.i \
        Functions/top_pull/expression='0.8*3e-3' \
        Executioner/end_time=17 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=false\
        Materials/CPStressUpdate/propsFile=bcc_props.in \
        Outputs/file_base=out_${temp}K_3e-3
done

for temp in {473..673..100}; do
    mpirun -n "$NPROC" $RHOCPOPT -i bcc_pxtal.i \
        Functions/top_pull/expression='0.8*3e-3' \
        Executioner/end_time=17 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=true\
        Materials/CPStressUpdate/propsFile=bcc_props_473.in \
        Functions/dts/y='0.0001    0.02' \
        Outputs/file_base=out_${temp}K_1.26e-3
done

# Strain rate: 1.26e-3 /s
for temp in {723..923..50}; do
    mpirun -n "$NPROC" ~/projects/rhocp-vr-250902/rhocp-opt -i bcc_pxtal.i \
        Functions/top_pull/expression='0.8*1.26e-3' \
        Executioner/end_time=17 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=false\
        Materials/CPStressUpdate/propsFile=bcc_props.in \
        Outputs/file_base=out_${temp}K_3e-3
done

for temp in {473..673..100}; do
    mpirun -n "$NPROC" $RHOCPOPT -i bcc_pxtal.i \
        Functions/top_pull/expression='0.8*3e-3' \
        Executioner/end_time=17 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=true\
        Materials/CPStressUpdate/propsFile=bcc_props_473.in \
        Functions/dts/y='0.0001    0.04' \
        Outputs/file_base=out_${temp}K_1.26e-3
done

# Strain rate: 1e-4 /s
for temp in {473..573..100}; do
    mpirun -n "$NPROC" $RHOCPOPT -i bcc_pxtal.i \
        Functions/top_pull/expression='0.8*1e-4' \
        Executioner/end_time=420 \
        Materials/CPStressUpdate/temp=$temp \
        Materials/CPStressUpdate/deltaH_eV=true\
        Materials/CPStressUpdate/propsFile=bcc_props_473.in \
        Functions/dts/y='0.0001    0.2' \
        Outputs/file_base=out_${temp}K_1.26e-3
done


