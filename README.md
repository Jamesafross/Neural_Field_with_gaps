# Neural_Field_with_gaps
JULIA :
The code for simulations were written in the Julia language.  A download link for Julia can be found at https://julialang.org/downloads/.
The simulations done for this paper were performed on the stable release at the time: version 1.5.1.
On top of julia, the following julia packages were used:

Plots,
SparseArrays,
ProgressMeter,
DifferentialEquations,
NPZ,
JLD.

PDE FOLDER:
The PDE folder contains two subfolder for simulations
in one spatial dimension and two spatial dimensions.
You will need to install Julia if you do not have it.
Julia download can be found here https://julialang.org/downloads/.  The simulations done for this paper were performed on the stable release at the time: version 1.5.1.
On top of julia, the following julia
packages were used in the PDE folder:

Plots,
SparseArrays,
ProgressMeter,
DifferentialEquations,
NPZ,
JLD.

These can be installed first by going into julia and
typing 'using Pkg' in the command line.
Then 'Pkg.add("PACKAGE_NAME")'.  For example, to add Plots and DifferentialEquations we would do the following:

julia> using Pkg
julia> Pkg.add("Plots")
julia> Pkg.add("DifferentialEquations")

MassModel FOLDER:
This folder includes code for simulating a one population mass model and a two population mass model in the 1_population and 2_population subfolders respectively.
Julia Packages used in these folders are:

Plots,
NPZ,
SparseArrays,
DifferentialEquations.
