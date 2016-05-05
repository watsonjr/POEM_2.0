# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_upwind_2D.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_grid_cp_2D.jld")
COBALT = load("./Data/JLD/Data_hindcast_000120.jld"); # 1980

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
#bio[ID] = 1.0e3*ones(Float64,size(ID));
bio[:,84:109] = 1.0e6; #seed equator
#bio[220:240,:] = 1.0e6; #seed Atl
#bio[59:79,:] = 1.0e6; #seed Pac
ni, nj = size(U);

const global DAYS = 365; # number of days

bio2D = open("./Data/CSV/bio_2Dadvect_test_eq.csv","w")

tstart = now()
for DAY = 1:DAYS
	println(DAY)
	U[ID] = 100.0.*COBALT["U"][:,DAY];
	V[ID] = 100.0.*COBALT["V"][:,DAY];

	bio = sub_advection(GRD,bio,U,V,ni,nj)
	biov=collect(bio[ID])
	#! Save
	writecsv(bio2D,biov')
end
#! Close
close(bio2D)
tend = now()
etime = Int(tend-tstart) / (1000 * 60 * 60) #elapsed time in hours
println(etime)