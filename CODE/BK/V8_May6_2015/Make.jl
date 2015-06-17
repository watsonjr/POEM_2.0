#### MODULES
using HDF5, JLD, Devectorize

####!! NUTS AND BOLTS
include("Parameters.jl")
include("sub_init.jl")
include("sub_functions.jl")
include("sub_model.jl")
include("Experiments.jl")

include("sub_routines_2.jl")
####!! EXPERIMENTS
@time make_local()




