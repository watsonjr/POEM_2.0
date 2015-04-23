# sub_model.jl
# running the model in time


################# POEM 2.0 ################
###### Packages to use
using PyPlot, NPZ
using Parameters



####### include sub routines
include("sub_feeding.jl")
include("sub_tau.jl")
include("size_based.jl")
include("size_structured.jl")
include("bio_change.jl")

PropEaten = sub_diet()
tauAll = sub_tau(PropEaten)

growth2 = fill(0.0, int(noGroups), int(P_Tmax))
Beaten2 = fill(0.0, int(noGroups), int(P_Tmax))
bio = fill(0.0, int(noGroups), int(P_Tmax))
bio[:, 1] = Binit
sumPrey = fill(0.0, int(noGroups), int(P_Tmax))
growth = fill(0.0, int(noGroups), int(P_Tmax))
Beaten = fill(0.0, int(noGroups), int(P_Tmax))
reprodLoss = fill(0.0, int(noGroups), int(P_Tmax))
reprodGain = fill(0.0, int(noGroups), int(P_Tmax))
matLoss = fill(0.0, int(noGroups), int(P_Tmax))
matGain = fill(0.0, int(noGroups), int(P_Tmax))

dBeaten2 = fill(0.0, int(noGroups), int(P_Tmax))
dgrowth2 = fill(0.0, int(noGroups), int(P_Tmax))
dsumPrey = fill(0.0, int(noGroups), int(P_Tmax))
dgrowth = fill(0.0, int(noGroups), int(P_Tmax))
dBeaten = fill(0.0, int(noGroups), int(P_Tmax))
dreprodLoss = fill(0.0, int(noGroups), int(P_Tmax))
dreprodGain = fill(0.0, int(noGroups), int(P_Tmax))
dmatLoss = fill(0.0, int(noGroups), int(P_Tmax))
dmatGain = fill(0.0, int(noGroups), int(P_Tmax))
dbio = fill(0.0, int(noGroups), int(P_Tmax))

model_time = [int(current_time):int(P_Tmax)]

function simulate()

for t = (model_time[1] + 1):P_Tmax
for i = 1:noGroups

	dgrowth2, dBeaten2 = size_based(growth2, Beaten2, i, t, PropEaten, bio)
	dgrowth, dBeaten, dreprodLoss, dreprodGain, dmatLoss, dmatGain = size_structured(tauAll, PropEaten, sumPrey, growth, Beaten, reprodLoss, reprodGain, matLoss, matGain, i, t, bio)
	dbio = Bchange(dgrowth, dgrowth2, dBeaten, dBeaten2, dreprodLoss, dreprodGain, dmatLoss, dmatGain, i, t) 

	bio[i, t] = bio[i, t-1] + (dbio[i, t] * dt)
end
end
end

simulate()
writecsv("testRes1.csv", hcat(model_time, bio'))
