#### MODULES
using HDF5, JLD, Devectorize


####!! NUTS AND BOLTS
include("Parameters.jl")
include("sub_init.jl")
include("sub_functions.jl")
include("sub_routines.jl")
include("Experiments.jl")


####!! SWITCHES
#make_local_switch = true

####!! EXPERIMENTS
@time make_local()




