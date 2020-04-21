function copyfield(f_arr::Vector{ScalarFieldNode{N}}) where N
    return map(copyfield, f_arr)
end

function copyfield(f::ScalarField{T,N}) where {T <: ComplexOrFloat, N}
    @debug "copy ScalarField "*f.name
    return ScalarField{T,N}(f.field,f.position,f.size,scaling_expr=f.scaling_expr,name=f.name)
end

function copyfield(f::ScalarFieldFunc{T,N}) where {T <: ComplexOrFloat, N}
    @debug "copy ScalarFieldFunc "*f.name
    return ScalarFieldFunc{T,N}(f.field,f.position,f.res,f.size,f.func_expr,f.gradx_expr,f.grady_expr,scaling_expr=f.scaling_expr,name=f.name)
end


function copyfield(f::VectorField{T,N}) where {T <: ComplexOrFloat,N}
    @debug "copy VectorField "*f.name
    return VectorField{T,N}(f.field,f.position,f.size,scaling_expr=f.scaling_expr,name=f.name)
end

function copyfield(f::ScalarFieldNode{N}) where N
    @debug "copy ScalarFieldNode "*f.name
    f_arr = map(copyfield,f.fields)
    fn = ScalarFieldNode{N}(f_arr,scaling_expr = f.scaling_expr,name=f.name)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end

function copyfield(f::VectorFieldNode{N}) where N
    @debug "copy VectorFieldNode "*f.name
    f_arr = map(copyfield,f.fields)
    fn = VectorFieldNode{N}(f_arr,scaling_expr = f.scaling_expr, name=f.name)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end
