# sub_tau.jl
# calculating the tau for entire pred-prey matrix

function sub_tau(PropEaten)
	
tauAll = zeros(int(noGroups), int(noGroups))

for i = 1:noGroups
for j = 1:noGroups

	if PropEaten[i, j] != 0 # excluding the non-dietary relationships
		tauAll[i, j] = t0*sizes[i]^83*sizes[j]^0.5 # large numbers --> check units!(0.83?)
	else
		tauAll[i, j] = 0
end
end
end
return tauAll
end
