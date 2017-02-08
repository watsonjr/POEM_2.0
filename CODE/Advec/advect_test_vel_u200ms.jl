# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_upwind_2D_u200ms.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_hindcast_grid_cp2D.jld")
COBALT = load("/Volumes/GFDL/POEM_JLD/Data_hindcast_vel200_000001.jld"); # yr3=1990; m/d

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
dep = GRD["Z"];
#bio[ID] = 1.0e6*ones(Float64,size(ID));
#bio[:,84:109] = 1.0e6; #seed equator
#bio[220:240,:] = 1.0e6; #seed Atl
#bio[59:79,:] = 1.0e6; #seed Pac
bio[5:25,:] = 1.0e6; #seed Indian W
#bio[340:360,:] = 1.0e6; #seed Indian E
#bio[:,181:200] = 1.0e6; #seed Arctic
#bio[:,12:32] = 1.0e6; #seed Antarctic
ni, nj = size(U);

const global DAYS = 365; # number of days

bio2D = open("/Volumes/GFDL/CSV/advect_tests/bio_2Dadvect_test_WInd_vel200_dt1hr.csv","w")

tstart = now()
for YR = 1#:2
	for DAY = 1:DAYS
		#println(DAY)
		println(YR," , ", DAY)
		ti = string(1000000+YR)
		COBALT = load(string("/Volumes/GFDL/POEM_JLD/Data_hindcast_vel200_",ti[2:end],".jld"));
		U[ID] = COBALT["U"][:,DAY] ./ (60 * 60 * 24); #m/d -> m/s
		V[ID] = COBALT["V"][:,DAY] ./ (60 * 60 * 24);

		bio = sub_advection(GRD,bio,U,V,ni,nj,dep)
		biov=collect(bio[ID])
		#! Save
		writecsv(bio2D,biov')
	end
end
#! Close
close(bio2D)
tend = now()
etime = Int(tend-tstart) / (1000 * 60 * 60) #elapsed time in hours
println(etime)