
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF
include("Parameters.jl")
include("sub_init.jl")
include("sub_routines.jl")
include("sub_functions.jl")
include("Experiments.jl")


####!! EXPERIMENTS
testoneloc = false
testalllocs_map = true
testalllocs_time = false

if testoneloc
	@time run_testoneloc()
end
if testalllocs_map
	@time run_testalllocs_map()
end
if testalllocs_time
	@time run_testalllocs_time()
end




