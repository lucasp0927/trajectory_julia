#####composition: for FieldNode object return composite field at time t, replace scaling with x->1.0
function composite{T<:Union{VectorField,ScalarField}}(f::T,t::Real) 
    T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end

function composite{N}(f::ScalarFieldNode{N},t::Real)
    #remember scaling
    output_type = typeoffield(f)
    geo = geometry(f)
    res = geo["res"]
    pos = geo["pos"]
    sz = geo["size"]
    arr_sz = floor(Integer,collect(sz)./collect(res))
    new_sz = arr_sz.*collect(res)
    output = Array(output_type,(arr_sz...))
    ####handle output
    ff = map(x->composite(x,t),f.fields)
    ####
    #ScalarField{Float64,N}(output,pos,new_sz;scaling = t->1.0)
end
