# type definitions
export ComplexOrFloat
export Field, VectorField, ScalarField
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract Field
type VectorField{T <: ComplexOrFloat, N} <: Field
    field::Array{T}
    position::Tuple{Vararg{Float64}}
    size::Tuple{Vararg{Float64}}
    scaling::Float64
    dim::Integer
    function VectorField(f::Array{T},pos::Tuple{Vararg{Real}},sz::Tuple{Vararg{Real}};scaling::Real = 1.0)
        length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,scaling,N):error("dimension error!")
    end
end

type ScalarField{T <: ComplexOrFloat,N} <: Field
    field::Array{T,N}
    position::Tuple{Vararg{Float64}}
    size::Tuple{Vararg{Float64}}
    scaling::Float64
    dim::Integer
    function ScalarField(f::Array{T,N},pos::Tuple{Vararg{Real}},sz::Tuple{Vararg{Real}};scaling::Real = 1.0)
        length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,scaling,N):error("dimension error!")
    end
end

abstract FieldNode <: Field
type VectorFieldNode <: FieldNode
    fields::Vector{Union{VectorFieldNode, VectorField}}
end

type ScalarFieldNode <: FieldNode
    fields::Vector{Field}
end
