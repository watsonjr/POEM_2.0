# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_upwind_2D.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_grid_cp_2D.jld")

#! Ocean currents
u200 = ncread("./Data/GCM/Hindcast/feb152013_run25_ocean.198801-200712_u200_v200.nc","Ut_200");
v200 = ncread("./Data/GCM/Hindcast/feb152013_run25_ocean.198801-200712_u200_v200.nc","Vt_200");
u200 = u200[:,:,1:12];
v200 = v200[:,:,1:12];

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
#bio[:,84:109] = 1.0e6; #seed equator
#bio[220:240,:] = 1.0e6; #seed Atl
bio[59:79,:] = 1.0e6; #seed Pac
#bio[5:25,:] = 1.0e6; #seed Indian W
#bio[340:360,:] = 1.0e6; #seed Indian E
#bio[:,181:200] = 1.0e6; #seed Arctic
#bio[:,12:32] = 1.0e6; #seed Antarctic
ni, nj, nt = size(u200);

const global MOS = 12; # number of mos
const global MNTH = collect([31,28,31,30,31,30,31,31,30,31,30,31])

bio2D = open("./Data/CSV/bio_2Dadvect_test_Pac_vel200_mo.csv","w")

tstart = now()
for MO = 1:MOS
	#vel200 m/s to m/d
	U = u200[:,:,MO] * 60 *60 *24;
	V = v200[:,:,MO] * 60 *60 *24;
	for DY = 1:MNTH[MO]
		println(MO," , ", DY)
		bio = sub_advection(GRD,bio,U,V,ni,nj)
		biov=collect(bio[ID])
		if (isnan(sum(biov)))
			break
		end
		#! Save
		writecsv(bio2D,biov')
	end
end
#! Close
close(bio2D)
tend = now()
etime = Int(tend-tstart) / (1000 * 60 * 60) #elapsed time in hours
println(etime)
