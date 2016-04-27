####type of field function
typeoffield{T<:ComplexOrFloat,N}(f::VectorField{T,N}) = T
typeoffield{T<:ComplexOrFloat,N}(f::ScalarField{T,N}) = T
function typeoffield{T<:FieldNode}(f::T)
    f_member_type = map(typeoffield,f.fields)
    if all(x->x<:ComplexOrFloat,f_member_type)
        typeintersect(reduce((x,y)->Union{x,y},f_member_type),ComplexOrFloat)==Float64?Float64:Complex{Float64}
    else
        error("Fields have to be either Complex{Float64} or Float64")
    end
end
