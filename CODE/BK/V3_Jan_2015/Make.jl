#### MODULES
using HDF5, JLD


####!! NUTS AND BOLTS
include("Parameters.jl")
include("sub_init.jl")
include("sub_functions.jl")
include("Experiments.jl")


####!! SWITCHES
#make_local_switch = true


####!! EXPERIMENTS
#! intialize
PRM_PI,PRM_PL,PRM_DE = make_parameters()
BIOMASS = sub_init_local(PRM_PI,PRM_PL,PRM_DE)

#! run model
make_local()




