# type definitions
export ComplexOrFloat
export Field, AbstractVectorField, AbstractScalarField, VectorField, ScalarField, VectorFieldNode, ScalarFieldNode, FieldNode
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract Field
abstract AbstractVectorField <: Field
abstract AbstractScalarField <: Field

type VectorField{T <: ComplexOrFloat, N} <: AbstractVectorField
    field::Array{T}
    position::Tuple{Vararg{Float64}}
    size::Tuple{Vararg{Float64}}
    scaling::Float64
    dim::Integer
    function VectorField(f::Array{T},pos::Tuple{Vararg{Real}},sz::Tuple{Vararg{Real}};scaling::Real = 1.0)
        length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,scaling,N):error("dimension error!")
    end
end

type ScalarField{T <: ComplexOrFloat,N} <: AbstractScalarField
    field::Array{T,N}
    position::Tuple{Vararg{Float64}}
    size::Tuple{Vararg{Float64}}
    scaling::Float64
    dim::Integer
    function ScalarField(f::Array{T,N},pos::Tuple{Vararg{Real}},sz::Tuple{Vararg{Real}};scaling::Real = 1.0)
        length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,scaling,N):error("dimension error!")
    end
end

type VectorFieldNode{N} <: AbstractVectorField
    fields::Vector{AbstractVectorField}
    dim::Integer
    function VectorFieldNode{T<:AbstractVectorField}(f::Vector{T}) 
        all(x->x.dim==N,f)?new(f,N):error("dimension error!")
    end
end

type ScalarFieldNode{N} <: AbstractScalarField
    fields::Vector{Field}
    dim::Integer
    function ScalarFieldNode{T<:Field}(f::Vector{T}) 
        all(x->x.dim==N,f)?new(f,N):error("dimension error!")
    end
end
FieldNode = Union{VectorFieldNode,ScalarFieldNode}

