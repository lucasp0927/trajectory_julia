using MAT
ComplexOrFloat = Union{Complex{Float64},Float64}
function copy_to_sharedarray!{T<:ComplexOrFloat,N}(arr::Array{T,N})
    return convert(SharedArray,arr)
end

function mat2sharedarray(filename,variable)
    file = matopen(filename)
    var = read(file, variable) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    var_s = copy_to_sharedarray!(var)
    var = 0
    gc()
    return var_s
end

function savemat(filename,variable,var_name)
    file = matopen(filename, "w")
    write(file, var_name, variable)
    close(file)
end
