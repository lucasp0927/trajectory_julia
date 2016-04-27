####type of field function
typeoffield{T<:ComplexOrFloat,N}(f::VectorField{T,N}) = T

typeoffield{T<:ComplexOrFloat,N}(f::ScalarField{T,N}) = T

function typeoffield(f::VectorFieldNode)
    f_member_type = map(typeoffield,f.fields)
    @assert all(x->x<:ComplexOrFloat,f_member_type) "Fields have to be either Complex{Float64} or Float64"
    typeintersect(reduce((x,y)->Union{x,y},f_member_type),ComplexOrFloat)==Float64?Float64:Complex{Float64}
end

function typeoffield(f::ScalarFieldNode)
    f_member_type = map(typeoffield,f.fields)
    @assert all(x->x<:ComplexOrFloat,f_member_type) "Fields have to be either Complex{Float64} or Float64"
    if all(x->typeof(x)<:ScalarField,f.fields)
        return typeintersect(reduce((x,y)->Union{x,y},f_member_type),ComplexOrFloat)==Float64?Float64:Complex{Float64}            
    else
        return Float64
    end
end
