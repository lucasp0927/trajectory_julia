
function copyfield{T<:ComplexOrFloat,N}(f::ScalarField{T,N})
    return ScalarField{T,N}(f.field,f.position,f.size,scaling=f.scaling,name=f.name)
end

function copyfield{T<:ComplexOrFloat,N}(f::VectorField{T,N})
    return VectorField{T,N}(f.field,f.position,f.size,scaling=f.scaling,name=f.name)
end

function copyfield{N}(f::ScalarFieldNode{N})
    f_arr = map(copyfield,f.fields)
    fn = ScalarFieldNode{N}(f_arr,scaling = f.scaling,name=f.name)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end

function copyfield{N}(f::VectorFieldNode{N})
    f_arr = map(copyfield,f.fields)
    fn = VectorFieldNode{N}(f_arr,scaling = f.scaling, name=f.name)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end
