
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NetCDF
include("Parameters.jl")
include("sub_init.jl")
include("sub_routines.jl")
include("sub_functions.jl")
include("Experiments.jl")

####!! EXPERIMENTS
testoneloc = false
spinup_pristine = false
forecast_pristine = true
spinup_fished = false
forecast_fished = false

if testoneloc
	@time Testoneloc()
end
if spinup_pristine
	@time Spinup_pristine()
end
if spinup_fished
	@time Spinup_fished()
end
if forecast_pristine
	@time Forecast_pristine()
end
if forecast_fished
	@time Forecast_fished()
end




