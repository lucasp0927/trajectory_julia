using HDF5
using SharedArrays
ComplexOrFloat = Union{Complex{Float64},Float64}
function copy_to_sharedarray!(arr::Array{T,N}) where {T<:ComplexOrFloat,N}
    return convert(SharedArray,arr)
end

function finalize_shared_array!(shared_array::SharedArray)
    foreach(shared_array.refs) do r
        @spawnat r.where finalize(fetch(r))
    end
    finalize(shared_array.s)
    finalize(shared_array)    
end

function file2sharedarray(filename,variable)
    # var = h5open(filename,"r") do file
    #     read(file,variable)
    # end
    fid = h5open(filename,"r")
    if isa(fid[variable],HDF5Group)
        var = read(fid,variable*"/real") .+ (read(fid,variable*"/imag").*1im)
    elseif isa(fid[variable],HDF5Dataset)
        var = read(fid,variable)
    else
        @error "unexpected HDF5 file structure"
    end
    close(fid)

    var_s = copy_to_sharedarray!(var)
    var = 0
    GC.gc()
    return var_s
end

function h5todict(filename)
    #read variables under /
    @assert filename[end-2:end] == ".h5"
    data = Dict()
    @info "Opening file "*filename
    fid = h5open(filename,"r")
    varnames = names(fid)
    for name in varnames
        if isa(fid[name],HDF5Group)
            data[name] = read(fid,name*"/real") .+ (read(fid,name*"/imag").*1im)
        elseif isa(fid[name],HDF5Dataset)
            data[name] = read(fid,name)
        else
            @error "unexpected HDF5 file structure"
        end
    end
    close(fid)
    return data
end

function dicttoh5(filename, data)
    @assert filename[end-2:end] == ".h5"
    varnames = collect(keys(data))
    fid = h5open(filename,"w")
    for name in varnames
        write_dataset(fid,name,data[name])
    end
    close(fid)
end

function write_dataset(fid,name,data::T) where T<:Number
    write(fid,name,data)
end

function write_dataset(fid,name,data::Array{T}) where T<:Real
    write(fid,name,data)
end

function write_dataset(fid,name,data::Array{T}) where T<:Complex
    write(fid,name*"/real",real(data))
    write(fid,name*"/imag",imag(data))
end
