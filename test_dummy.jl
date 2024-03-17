

## mockup data

using Pkg, UUIDs

# possibly multiple wavelengths
wavelength = collect(400:50:800)
Nl = length(wavelength)

# dummy 30x30 matrix values for each wavelength
tmatrix = reshape(repeat(collect(1:900.0) + collect(1:900) * 1im, Nl),
    (30, 30, Nl))

modes = (l=collect(1:30),
    m=collect(1:30),
    polarisation=repeat(["electric", "magnetic"], 15))

# dummy 'analytical zeros' for e.g. EBCM methods
azeros = collect(Iterators.product(1:2:30, 1:2:30))
analytical_zeros = (q=[z[1] for z in azeros[:]],
    qp=[z[2] for z in azeros[:]])

embedding = (relative_permeability=1.0,
    relative_permittivity=1.33^2)

materials = (embedding=embedding,
    Au=(relative_permeability=1.0,
        relative_permittivity=repeat([-11.4 + 1.181im], Nl)))

geometry = (shape="spheroid",
    radiusxy=20.0,
    radiusz=40.0)

# details about computation, including full script
computation = (method="EBCM",
    software="SMARTIES",
    version="1.1",
    Ntheta=40,
    accuracy=1e-10,
    analytical_zeros=analytical_zeros,
    script=read("test_dummy.jl", String))


pkgs = Pkg.Operations.Context().env.manifest
hdf5version = string(pkgs[findfirst(v -> v.name == "HDF5", pkgs)].version)

all = (vacuum_wavelength=wavelength,
    tmatrix=tmatrix,
    modes=modes,
    embedding=embedding,
    materials=materials,
    geometry=geometry,
    computation=computation,
    uuid=string(UUIDs.uuid4()))


## saving
using HDF5

# custom write_dicts_to_hdf5_groups()
include("_fun.jl")

h5open("aj.tmat.h5", "w") do fid

    # write all data in groups
    write_namedtuples_to_hdf5_groups(fid, all)

    # write custom attributes
    attributes(fid)["name"] = "Au prolate spheroid in water"
    attributes(fid)["created_with"] = "HDF5.jl"
    attributes(fid)["storage_format_version"] = hdf5version
    attributes(fid)["description"] = "Computation using SMARTIES, a numerically robust EBCM implementation for spheroids"
    attributes(fid)["keywords"] = "gold, spheroid, ebcm"
    attributes(fid["vacuum_wavelength"])["unit"] = "nm"
    attributes(fid["uuid"])["version"] = "4"
    attributes(fid["geometry"])["name"] = "prolate spheroid"

end


## manual step by step version

# embedding = create_group(fid, "embedding") 
# embedding["relative_permeability"] = 1.0
# embedding["relative_permittivity"] = 1.33^2

# materials = create_group(fid, "materials") 
# materials["embedding/relative_permeability"] = 1.0
# materials["embedding/relative_permittivity"] = 1.33^2
# materials["Au/relative_permeability"] = 1.0
# materials["Au/relative_permittivity"] = -11.4+1.181im

# geometry = create_group(fid, "geometry") 
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
# computation["analytical_zeros/q"] = [z[1] for z in analytical_zeros[:] ]
# computation["analytical_zeros/qp"] = [z[2] for z in analytical_zeros[:] ]
# computation["script"] = read("test_dummy.jl", String)

