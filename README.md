# &rho;-CP: Open Source Dislocation Density Based Crystal Plasticity Framework for Simulating Thermo-Mechanical Deformation
## Anirban Patra<sup>1</sup>, Suketa Chaudhary<sup>1</sup>, Namit Pai<sup>1</sup>, Tarakram Ramgopal<sup>1</sup>, Sarthak Khandelwal<sup>1</sup>, Adwitiya Rao<sup>1</sup>, David L. McDowell<sup>2,3</sup>
## <sup>1</sup>Department of Metallurgical Engineering and Materials Science, Indian Institute of Technology Bombay, Mumbai, India
## <sup>2</sup>School of Materials Science and Engineering, Georgia Institute of Technology, Atlanta, USA
## <sup>3</sup>GWW School of Mechanical Engineering, Georgia Institute of Technology, Atlanta, USA

&rho;-CP is a crystal plasticity solver that interfaces with the open source finite element solver, MOOSE (https://github.com/idaholab/moose), for crystal plasticity finite element modeling of anisotropic, heterogeneous deformation in polycrystalline ensembles. Source codes for the dislocation density-based crystal plasticity solver are provided in this repository, along with example applications for the thermo-mechanical deformation of Mg single and polycrystals, polycrystalline OFHC Cu and polycrystalline Ta.

## Installation
The user needs to install MOOSE first (https://mooseframework.inl.gov/getting_started/installation/index.html), then clone and compile &rho;-CP alongside MOOSE in the projects directory:
- git clone https://github.com/apatra6/rhocp.git
- make
