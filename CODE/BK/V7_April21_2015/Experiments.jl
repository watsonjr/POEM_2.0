

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
	PISC,PLAN,DETR,W = sub_init(PRM_PI,PRM_PL,PRM_DE,ID);
	
	#! Forward Euler integration
	#sub_spinup(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,COBALT,ID)
	#DT = 1.
	#for DAY = 1:10;
	#	DY = int(ceil(DAY))

	function test(PISC,PRM_PI,COBALT,ID)
		for X = 1:length(ID)
				# Metabolism
				#PISC.met[X,:] = sub_metabolism2(PRM_PI.N,PRM_PI.s,COBALT["Tp"][X]);
				sub_metabolism!(PLAN,X,PRM_PL,COBALT["Tp"]);
				#sub_metabolism!(DETR,X,PRM_DE,COBALT["Tp"]);

				# Handling time
				#sub_tau!(PISC,X,PRM_PI);
				#sub_tau!(PLAN,X,PRM_PL);
				#sub_tau!(DETR,X,PRM_DE);

				# Encounter rates
				#sub_enc_pipi!(PISC,X,PRM_PI)
				#sub_enc_pipl!(PISC,PLAN,X,PRM_PI,PRM_PL)
		end
	end

	@time test(PISC,PRM_PI,COBALT,ID);

	end

end
