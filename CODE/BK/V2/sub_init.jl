
#============== LOCAL EXP INITIAL CONDITIONS =============#
function sub_init_local(PRM_PI,PRM_PL,PRM_DE)

	#===== GCM DATA =====#
	#GRD = load("./Data/Data_grid_local.jld")
	GRD = load("./Data/Data_grid.jld") # Change this
	GRD["N"] = length(GRD["GRD_LON"]);

	#===== VARIABLES =====#
	#! Piscivore
	bio = ones(GRD["N"],PRM_PI.N) .* eps()
	met = ones(GRD["N"],PRM_PI.N) .* eps()
	tau = ones(GRD["N"],PRM_PI.N) .* eps()
	PISC = fish(bio,met,tau);

	#! Planktivore
	bio = ones(GRD["N"],PRM_PL.N) .* eps()
	met = ones(GRD["N"],PRM_PL.N) .* eps()
	tau = ones(GRD["N"],PRM_PL.N) .* eps()
	PLAN = fish(bio,met,tau);

	#! Detritivore
	bio = ones(GRD["N"],PRM_DE.N) .* eps()
	met = ones(GRD["N"],PRM_DE.N) .* eps()
	tau = ones(GRD["N"],PRM_DE.N) .* eps()
	DETR = fish(bio,met,tau);

	return PISC, PLAN, DETR
end


##============== GLOBAL EXP INITIAL CONDITIONS =============#
#function sub_init_local(PRM_PI,PRM_PL,PRM_DE)
#
#	#===== GCM DATA =====#
#	GRD = load("./Data/Data_grid_global.jld")
#
#	#===== VARIABLES =====#
#	#! Piscivore
#	bio = ones(GRD["N"],PRM_PI.N) .* eps()
#	met = ones(GRD["N"],PRM_PI.N) .* eps()
#	tau = ones(GRD["N"],PRM_PI.N) .* eps()
#	PISC = fish(bio,met,tau);
#
#	#! Planktivore
#	bio = ones(GRD["N"],PRM_PL.N) .* eps()
#	met = ones(GRD["N"],PRM_PL.N) .* eps()
#	tau = ones(GRD["N"],PRM_PL.N) .* eps()
#	PLAN = fish(bio,met,tau);
#
#	#! Detritivore
#	bio = ones(GRD["N"],PRM_DE.N) .* eps()
#	met = ones(GRD["N"],PRM_DE.N) .* eps()
#	tau = ones(GRD["N"],PRM_DE.N) .* eps()
#	DETR = fish(bio,met,tau);
#
#	return PISC, PLAN, DETR
#end
#


