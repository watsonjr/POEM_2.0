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
make_parameters()
COBALT = load("./Data/Data_000001.jld");    
GRD = load("./Data/Data_grid.jld")
ID = [1:48111]

PISC,PLAN,DETR,W = sub_init(ID);





###! Get COBALT data
function get_COBALT!(COBALT,ID,DY,TEMP_p,TEMP_b,ZOO,DZc,WI)
	## Get data
    TEMP_p[1]   = COBALT["Tp"][ID,DY]
    TEMP_b[1]   = COBALT["Tb"][ID,DY]
    ZOO[1,1]    = COBALT["Zm"][ID,DY]
    ZOO[2,1]    = COBALT["Zl"][ID,DY]
    DZc[1,1]    = COBALT["dZm"][ID,DY]
    DZc[2,1] 	= COBALT["dZl"][ID,DY]
    WI[1]       = COBALT["det"][ID,DY]
end

function sub_metabolism!(FISH,s,N,TEMP_p,U,activity,temp,basal)
	U = 3.9*s.^0.13
	activity = exp(0.03*U)
	basal = 0.0033*s.^-0.13
	temp = exp(0.0548*TEMP_p[1])
    for i = 1:N
        FISH.met[i] = basal[i] * temp * activity[i] * 5.258
    end
end

###! SERIAL: Run over all grid cells
function sub_allspace(PLAN,PRM_PL,COBALT,ID,TEMP_p,TEMP_b,ZOO,DZc,WI,U,activity,temp,basal)
	# run
	for X in ID # space

		#! COBALT information
		get_COBALT!(COBALT,X,DY,TEMP_p,TEMP_b,ZOO,DZc,WI);

		#! calculate a function
		sub_metabolism!(PLAN,PRM_PL.s,PRM_PL.N,TEMP_p,U,activity,temp,basal);

	end
end

####!! Test speed
## Preallocate
TEMP_p  = Array(Float32,1)
TEMP_b  = Array(Float32,1)
ZOO 	= Array(Float64,2)
DZc 	= Array(Float64,2)
WI  	= Array(Float64,1)
U 		= Array(Float64,1)
activity = Array(Float64,1)
temp 	= Array(Float64,1)
basal 	= Array(Float64,1)
@time sub_allspace(PLAN,PRM_PL,COBALT,ID,TEMP_p,TEMP_b,ZOO,DZc,WI,U,activity,temp,basal)


###! Learning about map and pmap
function easy(IN,DATA)
	OUT = Array(Float64,10)
	for i = 1:length(IN)
		OUT[i] = IN[i] * i * DATA
	end
	return OUT
end
a = cell(2)
DATA = [100, 200]# some data
a[1] = ones(10,1) # Could be the metabolism at each grid cell
a[2] = ones(10,1)
b = map(easy,a,DATA) # broadcasts function over all cells




###! Apply map to metabolism
ID = [1:48111] # number of grid cells
MET = cell(length(ID))
for i = 1:length(ID)
	MET[i] = Array(Float64,PL_N)
end

DY = 1
TEMP_p   = float64(COBALT["Tp"][:,DY])
temp = exp(0.0548*TEMP_p)
U = 3.9*PL_s.^0.13

const N = 15
const act = exp(0.03*U)
const bas = 0.0033*PL_s.^-0.13

function met(met,temp)
	for i = 1:N
		met[i] = bas[i] * temp * act[i] * 5.258
	end
	return met
end
@time map(met,MET,temp);


###! Apply to MODEL
TEMP_p   = float64(COBALT["Tp"][:,DY])
@devec temp = exp(0.0548.*TEMP_p)
function sub_metabolism!(met::Array{Float64},temp::Array{Float64})
	for i = 1:PL_N
		met[i] = PL_bas[i] * temp * PL_act[i] * 5.258
	end
	nothing
end
@time map(sub_metabolism!,PLAN.met,PLAN.tmet);



