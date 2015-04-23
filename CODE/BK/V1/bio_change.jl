# bio_change.jl
# incl the functions for size-structured and size-based

function Bchange(growth, growth2, Beaten, Beaten2, reprodLoss, reprodGain, matLoss, matGain, i, t)

	if orgTypes[i] < 4 # size-structured
		dbio[i, t] = dgrowth[i, t] - dBeaten[i, t] - dreprodLoss[i, t] + dreprodGain[i, t] - dmatLoss[i, t] + dmatGain[i, t] - mAll[i] - m0All[i]
	else # size-based
		dbio[i, t] = dgrowth2[i, t] - dBeaten2[i, t] - mAll[i] - m0All[i]
	end
return dbio
end
