# type definitions
export ComplexOrFloat
export Field, AbstractVectorField, AbstractScalarField, VectorField, ScalarField, VectorFieldNode, ScalarFieldNode, FieldNode
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract Field
abstract AbstractVectorField <: Field
abstract AbstractScalarField <: Field

type VectorField{T <: ComplexOrFloat, N} <: AbstractVectorField
    field::Array{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
    function VectorField(f::Array{T},pos::Vector{Float64},sz::Vector{Float64};scaling::Function = t->1.0)
        res = sz./(collect(size(f))[2:N+1]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,res,scaling,N):error("dimension error!")
    end
end

type ScalarField{T <: ComplexOrFloat,N} <: AbstractScalarField
    field::Array{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
    function ScalarField(f::Array{T,N},pos::Vector{Float64},sz::Vector{Float64};scaling::Function = t->1.0)
        res = sz./(collect(size(f))[1:N]-1)
        @assert all(x->x!=0,res) "zero resolution!"        
        length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,res,scaling,N):error("dimension error!")
    end
end

function setfield!{T<:ComplexOrFloat,N}(f::ScalarField{T,N},A::Array{T},pos::Vector{Float64},sz::Array{Float64};scaling::Function = t->1.0)
    res = sz./(collect(size(A))[1:N]-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A) "dimension error!"
    f.field = Array(T,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

function setfield!{T<:ComplexOrFloat,N}(f::VectorField{T,N},A::Array{T},pos::Vector{Float64},sz::Vector{Float64};scaling::Function = t->1.0)
    res = sz./(collect(size(A))[2:N+1]-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A)-1 "dimension error!"
    f.field = Array(T,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

type VectorFieldNode{N} <: AbstractVectorField
    fields::Vector{AbstractVectorField}
    scaling::Function
    dim::Integer
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    typeof::DataType
    function VectorFieldNode{T<:AbstractVectorField}(f::Vector{T};scaling::Function = t->1.0)
        @assert all(x->x.dim==N,f) "dimension error!"
        new(f,scaling,N,[],[],[],Complex{Float64})
    end
end

type ScalarFieldNode{N} <: AbstractScalarField
    fields::Vector{Field}
    scaling::Function
    dim::Integer
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    typeof::DataType    
    function ScalarFieldNode{T<:Field}(f::Vector{T};scaling::Function = t->1.0)
        @assert all(x->x.dim==N,f) "dimension error!"        
        new(f,scaling,N,[],[],[],Complex{Float64})
    end
end
FieldNode = Union{VectorFieldNode,ScalarFieldNode}

