
function sub_spinup(PISC,PLAN,DETR,W,COBALT,ID)

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
            
            ###! ticker
			DY 	= int(ceil(DAY))
            println(YR," , ", mod(DY,365))
			
	
			###! COBALT information
			get_COBALT!(COBALT,ID[X],X,DY);

			###! DEMOGRAPHIC CALCULATIONS
			#! metabolism
			map(sub_metabolism_pl!,PLAN.met,TEMP_p)

        end
    end
	
	### close save
	#close(Spinup_PISC)
    #close(Spinup_PLAN)
    #close(Spinup_DETR)
    #close(Spinup_W)

end

