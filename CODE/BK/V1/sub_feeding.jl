# sub_feeding.jl
# routine for size based feeding pref.

function sub_diet() # who eats whom? -takes habitat difr. into account

PropEaten = zeros(int(noGroups),int(noGroups))

for i = 1:3 # piscivores
for j = 1:noGroups
	if sizes[i] > sizes[j] 
		PropEaten[i, j] = 1/(stdev*(2*pi)^(1/2))*exp(-(log(sizes[i]) - log(sizes[j]) - PPMR)/(2*stdev^2))
	end
end
end

for i in [4, 5] #, 10, 11] # planktivores and zpl
for j in [1, 2, 3, 4, 5, 10, 11] # pelagic prey
	if sizes[i] > sizes[j]
		PropEaten[i, j] = 1/(stdev*(2*pi)^(1/2))*exp(-(log(sizes[i]) - log(sizes[j]) - PPMR)/(2*stdev^2))
	end
end
end

for i = 6:9 # benthos
for j = 6:9
	if sizes[i] > sizes[j]
		PropEaten[i, j] = 1/(stdev*(2*pi)^(1/2))*exp(-(log(sizes[i]) - log(sizes[j]) - PPMR)/(2*stdev^2))
	end
end
end

#for i = 10:11
#for j = 12
#	PropEaten[i, j]
#end
#end

return PropEaten
end

