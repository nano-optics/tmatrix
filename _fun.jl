using HDF5

function write_namedtuples_to_hdf5_groups(
    parent::Union{HDF5.File,HDF5.Group},
    tuple::NamedTuple
)
    for (key, value) in zip(keys(tuple), tuple)
        write_namedtuples_to_hdf5_groups(parent, key, value)
    end
end


function write_namedtuples_to_hdf5_groups(
    parent::Union{HDF5.File,HDF5.Group},
    key::Symbol,
    tuple::NamedTuple
)
    group = create_group(parent, string(key))
    write_namedtuples_to_hdf5_groups(group, tuple)
end


function write_namedtuples_to_hdf5_groups(
    parent::Union{HDF5.File,HDF5.Group},
    key::Symbol,
    value
)
    write(parent, string(key), value)
end


# earlier version with Dicts, but tuples are better memory-wise
# function write_dict(groupname, dict::Dict)
#     g = create_group(fid, groupname) 
#     for (key, value) in dict
#         keyname = groupname * "/" * key
#         if typeof(value) <: Dict # nested, recurse
#             write_dict(keyname, value)
#         else
#             write(fid, keyname, value)
#         end
#     end
# end

# # write top-level datasets
# fid["vacuum_wavelength"] = wavelength
# fid["tmatrix"] = tmatrix
# fid["uuid"] = uuid

# # write groups (some nested)
# write_dict("embedding", embedding)
# write_dict("materials", materials)
# write_dict("geometry", geometry)
# write_dict("modes", modes)
# write_dict("computation", computation)