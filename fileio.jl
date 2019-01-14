using HDF5
ComplexOrFloat = Union{Complex{Float64},Float64}
function copy_to_sharedarray!(arr::Array{T,N}) where {T<:ComplexOrFloat,N}
    return convert(SharedArray,arr)
end

function file2sharedarray(filename,variable)
    var = h5open(filename,"r") do file
        read(file,variable)
    end
    var_s = copy_to_sharedarray!(var)
    var = 0
    GC.gc()
    return var_s
end

function h5todict(filename)
    #read variables under /
    @assert filename[end-2:end] == ".h5"
    data = Dict()
    fid = h5open(filename,"r")
    varnames = names(fid)
    for name in varnames
        data[name] = read(fid,name)
    end
    close(fid)
    return data
end

function dicttoh5(filename, data)
    @assert filename[end-2:end] == ".h5"
    varnames = collect(keys(data))
    fid = h5open(filename,"w")
    for name in varnames
        write(fid,name,data[name])
    end
    close(fid)
end
# function savemat(filename,variable,var_name)
#     file = matopen(filename, "w")
#     write(file, var_name, variable)
#     close(file)
# end
