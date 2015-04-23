######## Functions for POEM 2.0
### size_structured.jl 

function size_structured(tauAll, PropEaten, sumPrey, growth, Beaten, reprodLoss, reprodGain, matLoss, matGain, i, t, bio)

# sum predated by i prey
for j = 1:noGroups
	if orgTypes[i] < 4 
		dsumPrey[i, t] = sumPrey[i, t] + tauAll[i, j]*alphaAll[i]*bio[j, t-1] # was Binit and sub_tau(sizes)
		sumPrey[i, t] = dsumPrey[i, t]
end
end

for j = 1:noGroups
	if orgTypes[i] < 4
		dgrowth[i, t] = growth[i, t] + PropEaten[i, j] * lambdaFish * alphaAll[i] * bio[i, t-1] * bio[j, t-1]/(1 + sumPrey[i, t-1])
		growth[i, t] = dgrowth[i, t]	
end
end

# Beaten (loss by predation)

for j = 1:noGroups
	if orgTypes[j] < 4 # also zpl and benthos do not eat fish (pred. on eggs could be an exception)
	dBeaten[i, t] = Beaten[i, t] + PropEaten[j, i] * alphaAll[j] * bio[j, t-1] * bio[i, t-1]/(1 + sumPrey[i, t-1]) 
	Beaten[i, t] = dBeaten[i, t]
	end
end

# reproduction loss
dreprodLoss[i, t] = reprod[i] * bio[i, t-1] # zero for zpl(5) and benthos(4) ; - SN: bio OR growth???


#reproduction gain
for j = 1:noGroups
	if i < noGroups && orgTypes[i] != orgTypes[i+1] && orgTypes[j] == orgTypes[i]
		dreprodGain[i, t] = reprodGain[i, t] + dreprodLoss[j, t-1]
		reprodGain[i, t] = dreprodGain[i, t]
	end
end

# maturation loss
dmatLoss[i, t] = matAll[i] * bio[i, t-1] # zero for zpl(5) and benthos(4)
# SN: bio or growth??

# maturation gain
	if i < noGroups && orgTypes[i] == orgTypes[i+1]
		dmatGain[i, t] = matLoss[i+1, t]
	end

return dgrowth, dBeaten, dreprodLoss, dreprodGain, dmatLoss, dmatGain

end
