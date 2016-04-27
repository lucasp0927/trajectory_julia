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
    res::Tuple{Vararg{Float64}}
    scaling::Function
    dim::Integer
    function VectorField(f::Array{T},pos::Tuple{Vararg{Real}},sz::Tuple{Vararg{Real}};scaling::Function = t->1.0)
        res = tuple((collect(sz)./(collect(size(f))[1:N]-1))...)
        if all(x->x!=0,res)
            length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,res,scaling,N):error("dimension error!")
        else
            error("zero resolution!")
        end
    end
end

type ScalarField{T <: ComplexOrFloat,N} <: AbstractScalarField
    field::Array{T,N}
    position::Tuple{Vararg{Float64}}
    size::Tuple{Vararg{Float64}}
    res::Tuple{Vararg{Float64}}
    scaling::Function
    dim::Integer
    function ScalarField(f::Array{T,N},pos::Tuple{Vararg{Real}},sz::Tuple{Vararg{Real}};scaling::Function = t->1.0)
        res = tuple((collect(sz)./(collect(size(f))[1:N]-1))...)
        if all(x->x!=0,res)
            length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,res,scaling,N):error("dimension error!")
        else
            error("zero resolution!")
        end
    end
end

type VectorFieldNode{N} <: AbstractVectorField
    fields::Vector{AbstractVectorField}
    scaling::Function
    dim::Integer
    function VectorFieldNode{T<:AbstractVectorField}(f::Vector{T};scaling::Function = t->1.0) 
        all(x->x.dim==N,f)?new(f,scaling,N):error("dimension error!")
    end
end

type ScalarFieldNode{N} <: AbstractScalarField
    fields::Vector{Field}
    scaling::Function
    dim::Integer
    function ScalarFieldNode{T<:Field}(f::Vector{T},scaling::Function = t->1.0) 
        all(x->x.dim==N,f)?new(f,scaling,N):error("dimension error!")
    end
end
FieldNode = Union{VectorFieldNode,ScalarFieldNode}

