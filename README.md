# &rho;-CP: Open Source Dislocation Density Based Crystal Plasticity Framework for Simulating Temperature- and Strain Rate-Dependent Deformation
### Anirban Patra<sup>1*</sup>, Suketa Chaudhary<sup>1</sup>, Namit Pai<sup>1</sup>, Tarakram Ramgopal<sup>1</sup>, Sarthak Khandelwal<sup>1</sup>, Adwitiya Rao<sup>1</sup>, David L. McDowell<sup>2,3**</sup>
#### <sup>1</sup>Department of Metallurgical Engineering and Materials Science, Indian Institute of Technology Bombay, Mumbai, India
#### <sup>2</sup>School of Materials Science and Engineering, Georgia Institute of Technology, Atlanta, USA
#### <sup>3</sup>GWW School of Mechanical Engineering, Georgia Institute of Technology, Atlanta, USA
#### <sup>*</sup>anirbanpatra@iitb.ac.in; anirbanpatra@gmail.com
#### <sup>**</sup>david.mcdowell@me.gatech.edu

&rho;-CP is a crystal plasticity solver that interfaces with the open source finite element solver, MOOSE (https://github.com/idaholab/moose), for crystal plasticity finite element modeling of anisotropic, heterogeneous deformation in polycrystalline ensembles. Source codes for the dislocation density-based crystal plasticity solver are provided in this repository. 

There are several constitutive models implemented for the different examples provided: 

(a) mobile and immobile dislocation density based crystal plasticity model, with threshold lattice resistance (DDCPStressUpdate, DDCPHCPStressUpdate) (Ref. [1])

(b) mobile and immobile dislocation density based crystal plasticity model, without threshold lattice resistance (DDCPTSTStressUpdate) (Ref. [2,9])

(c) statistically stored dislocation (SSD) density based Kocks-Mecking crystal plasticity model (DDCP_SSD_StressUpdate) (Ref. [3]) (code developed by Namit Pai)

(d) crystal plasticity model for hardening and creep under thermal and irradiation environments (ThermalIrradiationCPUpdate) (Ref. [5]) (code developed by Vikram Roy)

(e) mobile and immobile dislocation density based J<sub>2</sub> plasticity model (DDJ2StressUpdate) (Refs. [6,7])

Details of the framework and numerical implementation are available at:
https://doi.org/10.1016/j.commatsci.2023.112182  
https://arxiv.org/abs/2303.02441

Details of the material properties/model parameters and their input to the model are given in: <a href="rhoCP_model_parameters.pdf" target="_blank">Model Parameters</a>

Details of pre- and post-processing are given in: <a href="rho-CP_Pre_and_Post_Processing.pdf" target="_blank">Pre- and Post-Processing</a>

![Screenshot](rhoCP_graphic.png)

All input files tested with MOOSE version: 35c8c4a9a0 (2025-09-17), PETSc version: 3.23.4, SLEPc version: 3.23.2

Note: If you are using an older version of MOOSE, you may need to uncomment some lines in src/base/RhocpApp.C

## Installation
The user needs to install MOOSE first (https://mooseframework.inl.gov/getting_started/installation/index.html), then clone and compile &rho;-CP alongside MOOSE in the `projects` directory:
- Following installation of MOOSE and the required `conda` environment, the source files can be obtained either using the following commands from the `home` directory:   
`cd projects`  
`git clone https://github.com/apatra6/rhocp.git`  
or directly downloading the repository from github in the `projects` directory.
- The executable can be compiled using:  
`cd rhocp`  
`make -j 4`  
to get the executable `rhocp-opt` (here 4 represents the number of processors used for compiling and can be modified appropriately).
- If the user wishes to perform code developement and debug the application using `gdb`, the executable should be compiled in `debug` mode using the following coomand:  
`METHOD=dbg make -j 4`  
to get the executable `rhocp-dbg` (more details can be found at: https://mooseframework.inl.gov/application_development/debugging.html).

## Running Simulations
- The user is suggested to first go through the basics of running MOOSE simulations (https://mooseframework.inl.gov/getting_started/examples_and_tutorials/index.html).
- Example simulation files for magnesium, copper, tantalum, 304L stainless steel, DX54 ferritic steel, 316 and 317L austenitic stainless steel are located in the `examples` directory.
- The following input files are required to run a &rho;-CP simulation: (a) MOOSE input file, with `.i` extension, (b) slip system information file (`bcc_slip_sys.in`, for example), (c) material properties file (`bcc_props.in`, for example), (d) grain orientations in the form of Bunge Euler angles (`orientations.in`, for example). Additionally, the mesh may be: (i) created in the MOOSE input file itself, (ii) imported from an Exodus file (`64grains_512elements.e`, for example), or (iii) imported from an EBSD mesh file (`tantalum_input_original_euler.txt` in `examples/tantalum/EBSD_simulation`, for example). For the last case, Euler angles need not be imported separately.
- The EBSD mesh file can be created using DREAM3D. See: https://mooseframework.inl.gov/source/userobjects/EBSDReader.html and http://www.dream3d.io/2_Tutorials/EBSDReconstruction/ for additional details.
- Simulations can be run using the following example command:  
`mpiexec -n 4 ~/projects/rhocp/rhocp-opt -i Cu_compression_sim.i`  
for running the example given in  `rhocp/examples/copper/strain_rate_effect/compression_sr_1e-1ps/`.
- Output files in the form of `.csv` files can be used for plotting averaged values of various quantities and Exodus `.e` files can be visualized using Paraview (https://www.paraview.org/) for the deformation contours (the user is advised to use Paraview version 5.9 or lower).
- Spatio-temporal data can also be extracted from the `.e` output files using the Python SEACAS (https://github.com/sandialabs/seacas) libraries (an example script `extract_data.py` is provided in `examples/tantalum/temperature_effect/compression_512/298K_sr_5000_512grains`) or using the GUI-based data extraction tools in Paraview.

## References
For general details of the &rho;-CP framework and numerical implementation: Ref. [1]

For mobile and immobile dislocation density based crystal plasticity model, without threshold lattice resistance: Ref. [2]

For SSD-based Kocks-Mecking crystal plasticity model: Ref. [3]

For thermal Eigen strains and prediction of residual/internal strains: Refs. [3,4]

For crystal plasticity model for irradiation hardening and creep: Ref. [5]

For the dislocation density-based J<sub>2</sub> plasticity model: Refs. [6,7]

For numerical integration of the J<sub>2</sub> plasticity model: Ref. [8]

[1] Patra, A., Chaudhary, S., Pai, N., Ramgopal, T., Khandelwal, S., Rao, A., McDowell, D.L., “&rho;-CP: Open source dislocation density based crystal plasticity framework for simulating temperature- and strain rate-dependent deformation”, Computational Materials Science, Vol. 224, 2023, 112182.

[2] Patra, A., Tomé, C.N., “A dislocation density-based crystal plasticity constitutive model: Comparison of VPSC effective medium predictions with ρ-CP finite element predictions”, Modelling and Simulation in Materials Science and Engineering, Vol. 32, 2024, 045014.

[3] Pai, N., Samajdar, I., Patra, A., “Study of orientation-dependent residual strains during tensile and cyclic deformation of an austenitic stainless steel”, International Journal of Plasticity, Vol. 185, 2025, 104228.

[4] Pokharel, R., Patra, A., Brown, D.W., Clausen, B., Vogel, S.C., Gray, G.T., “An analysis of phase stresses in additively manufactured 304L stainless steel using neutron diffraction measurements and crystal plasticity finite element simulations”, International Journal of Plasticity, Vol. 121, 2019, pp. 201-217.

[5] Roy, V., Khan, I.A., Patra, A., "Crystal plasticity modeling of hardening and creep in ferritic - martensitic alloys under thermal and irradiation environments", International Journal of Plasticity, 2025, 104513.

[6] Khandelwal, S., Basu, S., Patra, A., “A machine learning-based surrogate modeling framework for predicting the history-dependent deformation of dual phase microstructures”, Materials Today Communications, Vol. 29, 2021, 102914.

[7] Basu, S., Patra, A., Jaya, B.N., Ganguly, S., Dutta, M., Samajdar, I., “Study of microstructure - property correlations in dual phase steels for achieving enhanced strength and reduced strain partitioning”, Materialia, Vol. 25, 2022, 101522.

[8] Patra, A., Pai, N., Sharma, P., “Modeling intrinsic size effects using dislocation density-based strain gradient plasticity”, Mechanics Research Communications, Vol. 127, 2023, 104038.

[9] Patra, A., Tomé, C.N., “Multiscale crystal plasticity modeling of deformation in an austenitic stainless steel”, Mechanics Research Communications, Vol. 148, 2025, 104490.
