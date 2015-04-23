# size_based.jl
# function for groups that are described by size-based functions
# i.e., benthos, zpl

function size_based(growth2, Beaten2, i, t, PropEaten, bio)

for j = 1:noGroups
	if orgTypes[i] > 3
		dgrowth2[i, t] = growth2[i, t] + PropEaten[i, j] * imax[i] * ((lambdaDet*bio[i, t-1]* bio[j, t-1])/(bio[j, t-1]+KAll[i]))
		growth2[i, t] = dgrowth2[i, t]
	end
end

## benthos and zpl predated by themselves, i.e., in the final need to add pred. by fish

for j = 1:noGroups
	if orgTypes[j] > 3
		dBeaten2[i, t] = Beaten2[i, t] + PropEaten[j, i] * imax[j] * ((bio[j, t-1]*bio[i, t-1])/(bio[i, t-1]+KAll[j]))
		Beaten2[i, t] = dBeaten[i, t]
	end
end

return dgrowth2, dBeaten2

end
