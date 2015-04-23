module Types

export Fish, Fishers, Output

#### Define Fish type
type Fish{T}
size::T
maturation::T # need to split
reproduction::T # need to split
production::T # need to split
metabolism::T # need to split
otherM::T # need to split
end

#### Define plankton type

#### Define benthos type

#### Define Fisher type
type Fishers
end

#### Output variable for plotting
type Output
end

end
