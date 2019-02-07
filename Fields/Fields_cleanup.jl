function cleanupfield(f::ScalarField{T,N}) where {T <: ComplexOrFloat, N}
    @debug "cleanup ScalarField "*f.name
    finalize_shared_array!(f.field)
end

function cleanupfield(f::VectorField{T,N}) where {T <: ComplexOrFloat,N}
    @debug "cleanup VectorField "*f.name
    finalize_shared_array!(f.field)
end

function cleanupfield(f::ScalarFieldNode{N}) where N
    @debug "cleanup ScalarFieldNode "*f.name
    map(cleanupfield,f.fields)
end

function cleanupfield(f::VectorFieldNode{N}) where N
    @debug "cleanup VectorFieldNode "*f.name
    map(cleanupfield,f.fields)
end
