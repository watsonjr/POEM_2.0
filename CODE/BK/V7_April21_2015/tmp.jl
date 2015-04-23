using HDF5, JLD, Devectorize

####!! NUTS AND BOLTS
include("Parameters.jl")
include("sub_init.jl")
include("sub_functions.jl")
include("sub_model.jl")
include("sub_routines.jl")

####!! RUN MODEL
DT = 1.; # time step
YEARS = 1; # integration period
DAY = 1
DY  = int(ceil(DAY))

####!! SETUP
PRM_PI,PRM_PL,PRM_DE = make_parameters()
COBALT = load("./Data/Data_000001.jld");    
GRD = load("./Data/Data_grid.jld")
ID = [1:48111]
PISC,PLAN,DETR,W = sub_init(PRM_PI,PRM_PL,PRM_DE,ID);


###! Get COBALT data
function get_COBALT(COBALT::Dict,ID::Int,X::Int,DY::Int)
    TEMP_p   = COBALT["Tp"][ID,DY]
    TEMP_b   = COBALT["Tb"][ID,DY]
    ZOO      = float64([COBALT["Zm"][ID,DY],COBALT["Zl"][ID,DY]])
    DZc      = float64([COBALT["dZm"][ID,DY],COBALT["dZl"][ID,DY]])
    WI       = float64(COBALT["det"][ID,DY])
    return TEMP_p::Float32, TEMP_b::Float32,
           ZOO::Array{Float64,1}, DZc::Array{Float64,1},WI::Float64
end

###! Run over all grid cells
function sub_allspace(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,COBALT,ID)
	for X = 1:length(ID) # space

		#! COBALT information
		TEMP_p, TEMP_b, ZOO, DZc, W.I[1] = get_COBALT(COBALT,ID[X],X,DY);

		#! calculate a function
		#sub_metabolism!(PLAN,X,PRM_PL,COBALT["Tp"]);
		PISC = sub_enc_pi(PISC,PLAN,DETR,ZOO,PRM_PI,X);

	end
end

sub_enc_pi(PISC,PLAN,DETR,ZOO,PRM_PI,X)

####!! Test speed
@time sub_allspace(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,COBALT,ID)

