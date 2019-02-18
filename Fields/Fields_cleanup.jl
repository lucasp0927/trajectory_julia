function cleanupfield!(f_arr::Vector{ScalarFieldNode{N}}) where N
    map(cleanupfield!,f_arr)
end

function cleanupfield!(f::ScalarField{T,N}) where {T <: ComplexOrFloat, N}
    finalize_shared_array!(f.field)
end

function cleanupfield!(f::VectorField{T,N}) where {T <: ComplexOrFloat,N}
    finalize_shared_array!(f.field)
end

function cleanupfield!(f::ScalarFieldNode{N}) where N
    map(cleanupfield!,f.fields)
end

function cleanupfield!(f::VectorFieldNode{N}) where N
    map(cleanupfield!,f.fields)
end
