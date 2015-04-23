
#============== LOCAL EXP INITIAL CONDITIONS =============#
function sub_init_local(PRM_PI,PRM_PL,PRM_DE)

	#===== VARIABLES =====#
	X = 0.17981663628808964; # mean Zm off Spain
	#! Piscivore
	bio_pisc = ones(PRM_PI.N) .* X
	bio_plan = ones(PRM_PL.N) .* X
	bio_detr = ones(PRM_DE.N) .* X
	bio_W = ones(1) .* X;
	BIOM = biom(bio_pisc,bio_plan,bio_detr,bio_W)
	return BIOM
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


