# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_upwind_2D.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_grid_cp_2D.jld")
#COBALT = load("./Data/JLD/Data_hindcast_000130.jld"); # 1990
#COBALT = load("./Data/JLD/Data_hindcast_surfvel_000120.jld"); # 1980
COBALT = load("./Data/JLD/Data_hindcast_vel200_000003.jld"); # 1990 m/d
#COBALT = load("./Data/JLD/Data_hindcast_velH200_000001.jld"); # yr3=1990 m/s

#bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
bio = 1.0e3*ones(Float64,GRD["Nlon"],GRD["Nlat"]);
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
#bio[ID] = 1.0e3*ones(Float64,size(ID));
#bio[:,84:109] = 1.0e6; #seed equator
#bio[220:240,:] = 1.0e6; #seed Atl
#bio[59:79,:] = 1.0e6; #seed Pac
#bio[5:25,:] = 1.0e6; #seed Indian W
#bio[340:360,:] = 1.0e6; #seed Indian E
#bio[:,181:200] = 1.0e6; #seed Arctic
#bio[:,12:32] = 1.0e6; #seed Antarctic
ni, nj = size(U);

const global DAYS = 365; # number of days

bio2D = open("./Data/CSV/advect_tests/bio_2Dadvect_test_global_vel200_dt1hr_j2.csv","w")

tstart = now()
for YR = 1 #:2
		for DAY = 1:DAYS
			println(YR," , ",DAY)
			U[ID] = COBALT["U"][:,DAY]; #m/d
			V[ID] = COBALT["V"][:,DAY];

			bio = sub_advection(GRD,bio,U,V,ni,nj)
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
