

#========= LOCAL EXPERIMENT =========#
#! size-based model run at one location
function make_local()

	##! Make parameters
	PRM_PI,PRM_PL,PRM_DE = make_parameters()

	####! Initialize variables
	PISC,PLAN,DETR,W = sub_init_local(PRM_PI,PRM_PL,PRM_DE)

	##! Load COBALT
	COBALT = load("./Data/Data_daily_averages.jld");
	GRD = load("./Data/Data_grid.jld")

	##! Extract a particular location
	XY = zeros(360,200);
	XY[GRD["ID"]] =[1:GRD["N"]]
	id = XY[272,156] # Iberian location
	
	#! get data
	Tp  = squeeze(COBALT["Tp"][id,:],1);
	Tb  = squeeze(COBALT["Tb"][id,:],1);
	Zm  = squeeze(COBALT["Zm"][id,:],1);
	Zl  = squeeze(COBALT["Zl"][id,:],1);
	dZm = squeeze(COBALT["dZm"][id,:],1);
	dZl = squeeze(COBALT["dZl"][id,:],1);
	det = squeeze(COBALT["det"][id,:],1)

	#! Setup data files to save
    Data_PISC = open("./Data/CSV/Data_PISC.csv","w")
    Data_PLAN = open("./Data/CSV/Data_PLAN.csv","w")
    Data_DETR = open("./Data/CSV/Data_DETR.csv","w")
    Data_W    = open("./Data/CSV/Data_W.csv","w")

	YEARS = 35;
	for YR = 1:YEARS # years
		for DY = 1:365# days

			#! ticker
			println(YR/YEARS," , ",DY/365)

			#! Write to file
			writecsv(Data_PISC,float32(PISC.bio'))
			writecsv(Data_PLAN,float32(PLAN.bio'))
			writecsv(Data_DETR,float32(DETR.bio'))
			writecsv(Data_W,float32(W.bio'))

			#! COBALT information
			TEMP_p = Tp[DY]
			TEMP_b = Tb[DY]
			ZOO    = [Zm[DY],Zl[DY]] ## zooplankton biomasses
			####### Get rid of time thing when rerun make_daily_data
			DZc    = [dZm[DY],dZl[DY]] .* 60 .* 60 .* 24## zooplankton mortality rates
			W.I[1] = det[DY] ## account for change in detrital pool here

			#! metabolism 
			PISC,PLAN,DETR = sub_metabolism(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE,TEMP_p)

			###! Encounter rates
			PISC,PLAN,DETR = sub_enc(PISC,PLAN,DETR,W,ZOO,PRM_PI,PRM_PL,PRM_DE)

			###! Total prey biomass encountered
			PISC,PLAN,DETR = sub_ENC(PISC,PLAN,DETR)

			###! Handling times
			PISC,PLAN,DETR = sub_tau(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)

			###! Consumption/Predation
			PISC,PLAN,DETR,W = sub_consume(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,ZOO,DZc)

			####! Energy available for growth
			PISC,PLAN,DETR = sub_nu(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)

			#! Energy available for somatic growth
			PISC,PLAN,DETR = sub_gamma(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)

			####! TOTAL RATES OF CHANGE (factoring in biomass densities)
			# eggs produced
			PISC,PLAN,DETR = sub_rep(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)

			#! total biomass to somatic growth
			PISC,PLAN,DETR = sub_grw(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)

            #! total biomass that matures to next size class
            PISC,PLAN,DETR = sub_mat(PISC,PLAN,DETR,PRM_PI,PRM_PL,PRM_DE)

            ###! MASS BALANCE
            PISC,PLAN,DETR,W = sub_update_bio(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE)

            ###! Forward Euler checks
            PISC,PLAN,DETR,W = sub_euler_checks(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE)

        end
    end
    close(Data_PISC)
    close(Data_PLAN)
    close(Data_DETR)
    close(Data_W)
end
