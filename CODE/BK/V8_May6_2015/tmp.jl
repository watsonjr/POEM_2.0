A = reshape(1:120, 15, 8)
h5write("./test2.h5", "mygroup2/A", A)

h5open("./test2.h5", "w") do file
	g = g_open(file, "mygroup2/A")
	g["A"] = A
end

h5open("test.h5", "w") do file
    g = g_create(file, "mygroup") # create a group
    g["dset1"] = 3.2              # create a scalar dataset inside the group
    attrs(g)["Description"] = "This group contains only a single dataset" # an attribute
end

fid = h5open("test.h5", "r+")
dset = fid["mygroup/dset1"]



data = h5read("./test2.h5", "mygroup2/A")


