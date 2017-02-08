# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_diff_upwind_2D_u200ms_loop.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_hindcast_grid_cp2D.jld")

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
#bio[ID] = 1.0e2*ones(Float64,size(ID));
#bio[:,84:109] = 1.0e2; #seed equator
bio[220:240,:] = 1.0e2; #seed Atl
#bio[59:79,:] = 1.0e2; #seed Pac
#bio[5:25,:] = 1.0e2; #seed Indian W
#bio[340:360,:] = 1.0e2; #seed Indian E
#bio[:,181:200] = 1.0e2; #seed Arctic
#bio[:,12:32] = 1.0e2; #seed Antarctic
ni=GRD["Nlon"];
nj=GRD["Nlat"];

K = 600.0;
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);

const global DAYS = 365; # number of days

bio2D = open("/Volumes/GFDL/CSV/advect_tests/bio_2Ddiff_test_Atl_dt15min_k600_b100.csv","w")

tstart = now()
for DAY = 1:DAYS
	println(DAY)

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
