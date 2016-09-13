# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_upwind_2D_uh200ms_swim.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_grid_cp_2D.jld")
#COBALT = load("./Data/JLD/Data_hindcast_000130.jld"); # 1990
#COBALT = load("./Data/JLD/Data_hindcast_surfvel_000120.jld"); # 1980
COBALT = load("./Data/JLD/Data_hindcast_velH200_000001.jld"); # yr3=1990; m2/s

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
bio[ID] = 1.0e6*ones(Float64,size(ID));
#bio[:,84:109] = 1.0e6; #seed equator
#bio[220:240,:] = 1.0e6; #seed Atl
#bio[59:79,:] = 1.0e6; #seed Pac
#bio[5:25,:] = 1.0e6; #seed Indian W
#bio[340:360,:] = 1.0e6; #seed Indian E
#bio[:,181:200] = 1.0e6; #seed Arctic
#bio[:,12:32] = 1.0e6; #seed Antarctic
ni, nj = size(U);

wgt = 2.5; 	#M=2.5; L=2500.0
#nu = -1.0 * GRD["Z"] #swim towards shallowest area
#nu = GRD["Z"] #swim towards deepest area
#nu = randn(ni,nj);

const global DAYS = 365; # number of days

bio2D = open("/Volumes/GFDL/CSV/advect_tests/bio_2Dadvect_swim_rand2_test_Global_velH200_dt1hr_j2_nodiv_divdepth.csv","w")

tstart = now()
for DAY = 1:DAYS
	println(DAY)
	U[ID] = COBALT["Uh"][:,DAY]; #m2/s
	V[ID] = COBALT["Vh"][:,DAY];
	U[ID] = U[ID]./GRD["Z"][ID] #m/s
	V[ID] = V[ID]./GRD["Z"][ID]

	nu = randn(ni,nj);

	bio = sub_advection_swim(GRD,bio,U,V,ni,nj,wgt,nu)
	biov=collect(bio[ID])
	#! Save
	writecsv(bio2D,biov')
end
#! Close
close(bio2D)
tend = now()
etime = Int(tend-tstart) / (1000 * 60 * 60) #elapsed time in hours
println(etime)
