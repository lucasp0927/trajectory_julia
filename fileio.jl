using MAT
function copy_to_sharedarray!{T<:ComplexOrFloat,N}(arr::Array{T,N})
    arr_s = SharedArray(T,size(arr))
    arr_s[:] = arr
    return arr_s
end

function mat2sharedarray(filename,variable)
    file = matopen(filename)
    var = read(file, variable) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    var_s = copy_to_sharedarray!(var)
    return var_s
end

function savemat(filename,variable,var_name)
    file = matopen(filename, "w")
    write(file, var_name, variable)
    close(file)
end

