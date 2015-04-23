
function sub_spinup(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,COBALT,ID)

	##! Setup data files to save
    #Spinup_PISC = open("./Data/CSV/Spinup_PISC.csv","w")
    #Spinup_PLAN = open("./Data/CSV/Spinup_PLAN.csv","w")
    #Spinup_DETR = open("./Data/CSV/Spinup_DETR.csv","w")
    #Spinup_W    = open("./Data/CSV/Spinup_W.csv","w")

	####! SPINUP 
    DT = 1.; # time step
    YEARS = 100; # integration period

    for YR = 1:YEARS # years

        for DAY = 1:1;#365/DT # days
            
            #! ticker
			DY 	= int(ceil(DAY))
            println(YR," , ", mod(DY,365))
			
			function sub_allspace(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,COBALT,ID)
				for X = 1:length(ID) # space
					
					##! Write to file
					#writecsv(Spinup_PISC,float32(PISC.bio'))
					#writecsv(Spinup_PLAN,float32(PLAN.bio'))
					#writecsv(Spinup_DETR,float32(DETR.bio'))
					#writecsv(Spinup_W,float32(W.bio'))

					#! COBALT information
					TEMP_p, TEMP_b, ZOO, DZc, W.I[1] = get_COBALT(COBALT,ID[X],X,DY);

					#! calculate demographics
					sub_demog(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,
							TEMP_p,TEMP_b,ZOO,DZc,X)
				end
			end
			@time sub_allspace(PISC,PLAN,DETR,W,PRM_PI,PRM_PL,PRM_DE,COBALT,ID);

			@time sub_enc_pi!(PISC,PLAN,DETR,ZOO,PRM_PI,X)

        end
    end
	
	### close save
	#close(Spinup_PISC)
    #close(Spinup_PLAN)
    #close(Spinup_DETR)
    #close(Spinup_W)

end

