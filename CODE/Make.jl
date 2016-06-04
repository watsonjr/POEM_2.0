
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT
include("Parameters.jl")
include("sub_init.jl")
include("sub_routines.jl")
include("sub_functions.jl")
include("Experiments.jl")


##! Test
#make_parameters(0) # make core parameters/constants
#COBALT = load("./Data/Data_000001.jld"); # if on laptop
#const global YEARS = 20; # integration period in years
#const global NX = 48111
#const global ID = collect(1:NX)
#Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, BENT = sub_init_fish(ID);
#ENVR = sub_init_env(ID);
#DY = 1
#get_COBALT!(COBALT,ID,DY,ENVR)
#JD = 1
#



####!! EXPERIMENTS
testoneloc = true
oneloc_hind_pristine = false
oneloc_fore_pristine = false
spinup_pristine = false
forecast_pristine = false
spinup_fished = false
forecast_fished = false

if testoneloc
	@time Testoneloc()
end
if oneloc_hind_pristine
	@time Oneloc_hindcast_pristine()
end
if oneloc_fore_pristine
	@time Oneloc_forecast_pristine()
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
