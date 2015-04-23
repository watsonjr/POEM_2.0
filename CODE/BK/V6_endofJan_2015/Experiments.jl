

#========= LOCAL EXPERIMENT =========#
#! size-based model run at one location
function make_local()

	##! Make parameters
	PRM_PI,PRM_PL,PRM_DE = make_parameters()

	##! Load COBALT
	#COBALT = load("./Data/Data_daily_averages.jld");
	println("Loading data")
	COBALT = load("./Data/Data_000001.jld");
	GRD = load("./Data/Data_grid.jld")

	##! Extract a particular location
	XY = zeros(360,200);
	XY[GRD["ID"]] =[1:GRD["N"]]
	#ID = XY[272,156] # Iberian location
	#ID = [XY[272,156],XY[272,156],XY[272,156]] # Iberian location
	ID = [1:48111]

	####! Initialize variables
	PISC,PLAN,DETR,W = sub_init(PRM_PI,PRM_PL,PRM_DE,ID)
	
	#! Forward Euler integration
	sub_spinup(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,COBALT,ID)

end
