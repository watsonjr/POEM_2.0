# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

#include("Advect_diff_upwind_2D_u200ms_merge.jl")
include("Advect_diff_upwind_2D_u200ms_loop.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_hindcast_grid_cp2D.jld");
COBALT = load("/Volumes/GFDL/POEM_JLD/Data_hindcast_vel200_000001.jld"); # yr3=1990; m/d

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
bio[ID] = 1.0e6*ones(Float64,size(ID));
#bio[:,84:109] = 1.0e6; #seed equator
#bio[220:240,:] = 1.0e6; #seed Atl
#bio[59:79,:] = 1.0e6; #seed Pac
#bio[5:25,:] = 1.0e6; #seed Indian W
#bio[340:360,:] = 1.0e6; #seed Indian E
#bio[:,181:200] = 1.0e6; #seed Arctic
#bio[:,12:32] = 1.0e6; #seed Antarctic
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
dep = GRD["Z"];
ni=GRD["Nlon"];
nj=GRD["Nlat"];

K = 600.0;

const global DAYS = 365; # number of days

bio2D = open("/Volumes/GFDL/CSV/advect_tests/bio_2Dadvec_diff_test_global_dt1hr_k600_nt1.csv","w")

tstart = now()
for DAY = 1:DAYS
	println(DAY)
	U[ID] = COBALT["U"][:,DAY] ./ (60 * 60 * 24); #m/d -> m/s
	V[ID] = COBALT["V"][:,DAY] ./ (60 * 60 * 24);
	bio = sub_advec_diff(GRD,bio,K,U,V,ni,nj)
	biov=collect(bio[ID])
	#! Save
	writecsv(bio2D,biov')
end
#! Close
close(bio2D)
tend = now()
etime = Int(tend-tstart) / (1000 * 60 * 60) #elapsed time in hours
println(etime)
