# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_diff_upwind_2D_u200ms.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_hindcast_grid_cp2D.jld")

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
#bio[ID] = 1.0e6*ones(Float64,size(ID));
#bio[:,84:109] = 1.0e6; #seed equator
bio[220:240,:] = 1.0e6; #seed Atl
#bio[59:79,:] = 1.0e6; #seed Pac
#bio[5:25,:] = 1.0e6; #seed Indian W
#bio[340:360,:] = 1.0e6; #seed Indian E
#bio[:,181:200] = 1.0e6; #seed Arctic
#bio[:,12:32] = 1.0e6; #seed Antarctic
ni=GRD["Nlon"];
nj=GRD["Nlat"];

K = 100.0;

const global DAYS = 365; # number of days

bio2D = open("/Volumes/GFDL/CSV/advect_tests/bio_2Ddiff_test_Atl_dt1hr_v1.csv","w")

tstart = now()
for DAY = 1:DAYS
	println(DAY)

	bio = sub_diffuse(GRD,bio,K,ni,nj)
	biov=collect(bio[ID])
	#! Save
		writecsv(bio2D,biov')
end
#! Close
close(bio2D)
tend = now()
etime = Int(tend-tstart) / (1000 * 60 * 60) #elapsed time in hours
println(etime)
