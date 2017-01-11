# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF, MAT

include("Advect_upwind_2D_uh200ms_swim.jl")

ID = load("./Data/Data_grid_hindcast_NOTflipped.jld","ID");
GRD = load("./Data/Data_grid_cp_2D.jld")
#COBALT = load("./Data/JLD/Data_hindcast_000130.jld"); # 1990
#COBALT = load("./Data/JLD/Data_hindcast_surfvel_000120.jld"); # 1980
COBALT = load("/Volumes/GFDL/POEM_JLD/Data_hindcast_velH200_1990.jld"); # yr3=1990; m2/s
COB2 = load("/Volumes/GFDL/POEM_JLD/esm2m_hist/Data_ESM2Mhist_1990.jld");

bio = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
prey = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
U = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
V = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
dep = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
dep[ID] = GRD["Z"][ID];
bio[ID] = 1.0e6*ones(Float64,size(ID));
#bio[:,84:109] = 1.0e6; #seed equator
#bio[220:240,:] = 1.0e6; #seed Atl
#bio[59:79,:] = 1.0e6; #seed Pac
#bio[5:25,:] = 1.0e6; #seed Indian W
#bio[340:360,:] = 1.0e6; #seed Indian E
#bio[:,181:200] = 1.0e6; #seed Arctic
#bio[:,12:32] = 1.0e6; #seed Antarctic
# PROBABLY NEED TO MULT BY LMASK TO GET ONLY OCEAN
bio = bio .* GRD["lmask"][:,:,1];

ni, nj = size(U);

#! Calc swimming speed (m/s)
#T = (Tp.*tdif) + (Tb.*(1.0-tdif))
# wgt = 2.5; 	#M=2.5; L=2500.0
# T = 15.0;
# w = exp(0.063*(T100-15.0)) * 0.5*L_m*1e-3;
# Q = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
# Q[ID] = w;
#L=10^((log10(20)+log10(200))/2); #medium
L=10^((log10(200)+log10(2000))/2); #large
Q = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);
T = zeros(Float64,GRD["Nlon"],GRD["Nlat"]);

#nu = -1.0 * dep #swim towards shallowest area
nu = dep #swim towards deepest area
#nu = randn(ni,nj);

const global DAYS = 365; # number of days

bio2D = open("/Volumes/GFDL/CSV/advect_tests/bio_2Dadvect_swim_deep_global_vel0_dt1hr_sep_lBLtemp.csv","w")
#prey2D = open("/Volumes/GFDL/CSV/advect_tests/prey_2Dadvect_swim_Zl_global_vel0_dt1hr_sep_mBLtemp.csv","w")

tstart = now()
for DAY = 1:DAYS
	println(DAY)
	# U[ID] = COBALT["Uh"][:,DAY]; #m2/s
	# V[ID] = COBALT["Vh"][:,DAY];
	T[ID] = COB2["Tp"][:,DAY];
	Q = exp(0.063*(T-15.0)) .* 0.5.*L.*1e-3 .* GRD["lmask"][:,:,1];
	prey[ID] = COB2["Zl"][:,DAY];
	#nu = prey;

	#nu = randn(ni,nj);
	bio = sub_advection_swim(GRD,bio,U,V,ni,nj,Q,nu,dep)
	biov=collect(bio[ID])
	nuv=collect(nu[ID])
	#! Save
	if (length(biov) == 48111)
		writecsv(bio2D,biov')
		#writecsv(prey2D,nuv')
	else
		println("biov != 48111")
		break
	end
end
#! Close
close(bio2D)
tend = now()
etime = Int(tend-tstart) / (1000 * 60 * 60) #elapsed time in hours
println(etime)
