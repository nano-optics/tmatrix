using HDF5

include("_fun.jl")




# simple objects
a = collect(400:50:800)
Nl = length(a)

m = reshape(repeat(collect(1:900.0), Nl), (30, 30, Nl))


g = (shape="spheroid",
    a=20.0,
    c=40.0)


# this one is nested
e = (
    medium=(epsilon=1.0, mu=1.0),
    particle=(epsilon=1.0, mu=1.0))


allthestuff = (a=a, m=m, g=g, e=e)

## ideal situation



# for (k,v) in zip(keys(e), e)
#     println(k, " => ", v)
# end


h5open("allthestuff.h5", "w") do h5f
    write_namedtuples_to_hdf5_groups(h5f, allthestuff)
end








f = "debug.h5"
fid = h5open(f, "w")
# save_dict_toh5(fid, allthestuff)
fid["a"] = e
close(fid)

## manual process

f = "debug.h5"
fid = h5open(f, "w")


fid["a"] = a
fid["m"] = m

group1 = create_group(fid, "g")
for (key, value) in g
    g[key] = value
end

group2 = create_group(fid, "e")
for (key, value) in e
    # ridiculous way to handle nested group
    if (typeof(value) <: Dict)
        subgroup = create_group(fid, "e" * "/" * key)
        for (subkey, subvalue) in value
            subgroup[subkey] = subvalue
        end
    else
        group2[key] = value
    end
end


attributes(fid)["what's this"] = "some attributes"
attributes(fid["e"])["unit"] = "kg"


close(fid)



function assign_dict(groupname, dict::Dict)

    g = create_group(fid, groupname)

    for (key, value) in dict

        keyname = groupname * "/" * key

        if typeof(value) <: Dict
            assign_dict(keyname, value)
        else
            write(fid, keyname, value)
        end
    end
end


f = "debug_auto.h5"
fid = h5open(f, "w")


fid["a"] = a
fid["m"] = m

# group1 = create_group(fid, "g") 
assign_dict("g", g)

# group2 = create_group(fid, "e") 
assign_dict("e", e)

close(fid)

fid

# assign_dict(g, geometry)

# geometry["shape"] = "spheroid"
# geometry["radiusxy"] = 20.0
# geometry["radiusz"] = 40.0




# modes = create_group(fid, "modes") 
# modes["l"] = collect(1:30)
# modes["m"] = collect(1:30)
# modes["polarization"] = repeat(["electric","magnetic"],  15)


# computation = create_group(fid, "computation") 
# computation["method"] = "EBCM"
# computation["software"] = "SMARTIES"
# computation["version"] = "1.1"
# computation["Ntheta"] = 40
# computation["accuracy"] = 1e-10
# computation["analytical_zeros/p"] = [z[1] for z in analytical_zeros[:] ]
# computation["analytical_zeros/pp"] = [z[2] for z in analytical_zeros[:] ]

# computation["script"] = read("test_dummy.jl", String)
# fid["uuid"] = string(UUIDs.uuid4())

# using Pkg

# pkg_name = "HDF5"
# ver = Pkg.Operations.Context().env.manifest


# attributes(fid)["name"] = "Au prolate spheroid in water"
# attributes(fid)["created_with"] = "HDF5.jl"
# attributes(fid)["storage_format_version"] = string(ver[findfirst(v->v.name == pkg_name, ver)].version)
# attributes(fid)["description"] = "Computation using SMARTIES, a numerically robust EBCM implementation for spheroids"
# attributes(fid)["keywords"] = "gold, spheroid, ebcm"
# attributes(fid["vacuum_wavelength"])["unit"] = "nm"
# attributes(fid["uuid"])["version"] = "4"
# attributes(fid["geometry"])["name"] = "prolate spheroid"


# close(fid)


# using FileIO
# e = Dict(
#     "medium" => Dict("epsilon" => 1.0, "mu" => 1.0),
#     "particle" => Dict("epsilon" => 1.0, "mu" => 1.0))


# g = Dict("shape" => "spheroid", 
# "a" => 20.0,
# "c" => 40.0)


# save("testgroup.h5", g)

