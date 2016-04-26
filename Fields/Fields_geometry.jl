function geometry{T<:FieldNode}(f::T)
    #find the new optimal position and size
    #minimum(cat(2,[1,2,3],[2,3,4]),2)
    minpos::Vector{Float64} = vec(minimum(cat(2,map(x->collect(geometry(x)[1]),f.fields)...),2))
    maxpos::Vector{Float64} = vec(maximum(cat(2,map(x->(collect(geometry(x)[1]).+collect(geometry(x)[2])),f.fields)...),2))
    return tuple(minpos...), tuple((maxpos-minpos)...)
end

geometry{T<:Union{VectorField,ScalarField}}(f::T) = f.position,f.size
