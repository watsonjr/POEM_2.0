# TEST ADVECTION
####!! NUTS AND BOLTS
using HDF5, JLD, Devectorize, NPZ, NetCDF
include("Parameters.jl")
include("sub_init.jl")
include("sub_routines.jl")
include("sub_functions.jl")
include("Experiments.jl")
include("JW_adv_diff_update.jl")


make_parameters(0) # make core parameters/constants
#! setup spinup (loop first year of COBALT)
COBALT = load("./Data/JLD/Data_hindcast_PC_000120.jld"); # 1980
#! Add phenology params from csv file with ID as row
Tref = readdlm("./Data/grid_phenol_T0raw_NOflip.csv",','); #min temp for each yr at each location
#global TrefP = readdlm("./Data/grid_phenol_T0p_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for upper 100m
#global TrefB = readdlm("./Data/grid_phenol_T0b_clim_min_NOflip.csv",','); #1901-1950 climatological min temp at each location for bottom
global TrefP = Tref;
global TrefB = Tref;
global Dthresh = readdlm("./Data/grid_phenol_DTraw_NOflip.csv",',');
global Sp = readdlm("./Data/Gaussian_spawn_2mo.csv",',');
global GRD = load("./Data/Data_grid_hindcast_NOTflipped.jld");
YEARS = 50
DAYS = 365
#! choose where and when to run the model
const global NX = 48111
const global ID = collect(1:NX);

#! Initialize
phen=0;
Sml_f, Sml_p, Sml_d, Med_f, Med_p, Med_d, Lrg_p, Lrg_d, BENT = sub_init_fish(ID,phen);
Med_d.td[1:NX] = 0.0;
Lrg_d.td[1:NX] = 0.0;
ENVR = sub_init_env(ID);

bio = zeros(Float64,NX);
#seed = readdlm("./Data/Eq_ids.csv",','); #seed equator
seed = readdlm("./Data/Atl_ids.csv",','); #seed Atl
#seed = readdlm("./Data/Pac_ids.csv",','); #seed Pac
#seed = readdlm("./Data/EInd_ids.csv",','); #seed Indian W
#seed = readdlm("./Data/WInd_ids.csv",','); #seed Indian E
#seed = readdlm("./Data/Arc_ids.csv",','); #seed Arctic
#seed = readdlm("./Data/AntArc_ids.csv",','); #seed Antarctic
seed = round(Int64,seed);
bio[seed] = 1.0;

#Horizontal diffusivity m/s -> m/d
A = 1.0e-4 * 60 * 60 * 24;

#! Internal time step for advection (in days)
global dtime = (1/24.0)
bio2D = open("/Volumes/GFDL/NC/AdvectTests/bio_JWadvect_test_Atl_noadvec_1hr_noreflect_negs.csv","w")

tstart = now()
writecsv(bio2D,bio')
###################### Run the Model
#! Run model with no fishing
for YR = 1#:YEARS # years
	for DAY = 1:DT:DAYS # days
		###! Future time step
		DY  = Int(ceil(DAY))
		println(YR," , ", mod(DY,365))
		# Run biology to get nu values
		sub_futbio!(ID,DY,COBALT,ENVR,Sml_f,Sml_p,Sml_d,Med_f,Med_p,Med_d,Lrg_p,Lrg_d,BENT);

		#! Internal time step for advection
		for time = dtime:dtime:1
			# Use adult forage fish nus and swimming speed
			bio = sub_advection(bio,Med_f.nu,ENVR.U,ENVR.V,GRD["dxtn"],GRD["dyte"],ENVR.Tp,ENVR.Tb,Med_f.td,M_m)
			#bio = sub_diffuse(bio,A,GRD["dxtn"],GRD["dyte"])
		end
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
