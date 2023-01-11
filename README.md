# &rho;-CP: Open Source Dislocation Density Based Crystal Plasticity Framework for Simulating Thermo-Mechanical Deformation
## Anirban Patra<sup>1</sup>, Suketa Chaudhary<sup>1</sup>, Namit Pai<sup>1</sup>, Tarakram Ramgopal<sup>1</sup>, Sarthak Khandelwal<sup>1</sup>, Adwitiya Rao<sup>1</sup>, David L. McDowell<sup>2,3</sup>
## <sup>1</sup>Department of Metallurgical Engineering and Materials Science, Indian Institute of Technology Bombay, Mumbai, India
## <sup>2</sup>School of Materials Science and Engineering, Georgia Institute of Technology, Atlanta, USA
## <sup>3</sup>GWW School of Mechanical Engineering, Georgia Institute of Technology, Atlanta, USA

&rho;-CP is a crystal plasticity solver that interfaces with the open source finite element solver, MOOSE (https://github.com/idaholab/moose), for crystal plasticity finite element modeling of anisotropic, heterogeneous deformation in polycrystalline ensembles. Source codes for the dislocation density-based crystal plasticity solver are provided in this repository, along with example applications for the thermo-mechanical deformation of Mg single and polycrystals, polycrystalline OFHC Cu and polycrystalline Ta.

## Installation
The user needs to install MOOSE first (https://mooseframework.inl.gov/getting_started/installation/index.html), then clone and compile &rho;-CP alongside MOOSE in the projects directory:
- The source files can be obtained either using the following command: `git clone https://github.com/apatra6/rhocp.git` or directly downloading the repository from github.
- The executable can be compiled using: `make -j 4` to get the executable `rhocp-opt` (here 4 represents the number of processors used for compiling and can be modified appropriately).
- If the user wishes to perform code developement and debug the application using `gdb`, the executable should be compiled in `debug` mode using the following coomand: `METHOD=dbg make -j 4` to get the executable `rhocp-dbg` (more details can be found at: https://mooseframework.inl.gov/application_development/debugging.html).

## Running Simulations
- The user is advised to first go through the basics of running MOOSE simulations (https://mooseframework.inl.gov/getting_started/examples_and_tutorials/index.html).
- Example simulation files for Mg, Cu and Ta are located in the `examples` directory.
- In order to run a &rho;-CP simulation, the following input files are required: (a) MOOSE input file, with `.i` extension, (b) slip system information file (`bcc_slip_sys.in`, for example), (c) material properties file (`bcc_props.in`, for example), (d) grain orientations in the form of Bunge Euler angles (`orientations.in`, for example). Additionally, the mesh may be (i) created in the MOOSE input file itself, (ii) imported from an Exodus file (`64grains_512elements.e`, for example), or (iii) imported from an EBSD mesh file (`tantalum_input_ms.txt`, for example). For the last case, Euler angles need not be imported separately.
- The EBSD mesh file can be created using DREAM3D. See: https://mooseframework.inl.gov/source/userobjects/EBSDReader.html and http://www.dream3d.io/2_Tutorials/EBSDReconstruction/ for additional details.
- Simulations can be run using the following example command: `mpiexec -n 4 ~/projects/rhocp/rhocp-opt -i Cu_compression_sim.i` for running the example given in  `rhocp/examples/copper/compression_sr_1e-1ps/`.
- Output files in the form of `.csv` files can be used for plotting averaged values of various quantities and `.e` files can be visualized using Paraview (https://www.paraview.org/) for the deformation contours (the user is advised to use Paraview version 5.9 or lower).
- Spatio-temporal data can also be extracted from the `.e` output files using the Python SEACAS (https://github.com/sandialabs/seacas) libraries (an example script `extract_data.py` is provided) or using the GUI-based data extraction tools in Paraview.
