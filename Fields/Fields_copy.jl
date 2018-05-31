function copyfield{T<:ComplexOrFloat,N}(f::ScalarField{T,N})
    @debug "copy ScalarField "*f.name
    return ScalarField{T,N}(f.field,f.position,f.size,scaling_expr=f.scaling_expr,name=f.name)
end

function copyfield{T<:ComplexOrFloat,N}(f::VectorField{T,N})
    @debug "copy VectorField "*f.name
    return VectorField{T,N}(f.field,f.position,f.size,scaling_expr=f.scaling_expr,name=f.name)
end

function copyfield{N}(f::ScalarFieldNode{N})
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

function copyfield{N}(f::VectorFieldNode{N})
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
