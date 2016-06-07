# type definitions
export ComplexOrFloat
export Field, AbstractVectorField, AbstractScalarField, VectorField, ScalarField, VectorFieldNode, ScalarFieldNode, FieldNode
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract Field
abstract AbstractVectorField <: Field
abstract AbstractScalarField <: Field
#TODO: right now sample is only for 2d, also the ScalarFieldNode.sample datatype is Float64
type VectorField{T <: ComplexOrFloat, N} <: AbstractVectorField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
    sample::Array{T,3}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Complex{Float64}
    name::ASCIIString
    function VectorField(f::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling::Function =  t->1.0, name::ASCIIString = "VectorField")
        res = sz./(collect(size(f))[2:N+1]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,res,scaling,N,zeros(T,(3,4,4)),[0.0,0.0],[0,0,0,0],zero(Complex{Float64}),name):Lumberjack.error("dimension error!")
    end
end

type ScalarField{T <: ComplexOrFloat,N} <: AbstractScalarField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
    sample::Array{T,2}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Float64
    name::ASCIIString
    function ScalarField(f::SharedArray{T,N},pos::Vector{Float64},sz::Vector{Float64};scaling::Function =  t->1.0, name::ASCIIString = "ScalarField")
        res = sz./(collect(size(f))[1:N]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,res,scaling,N,zeros(T,(4,4)),[0.0,0.0],[0,0,0,0],zero(Float64),name):Lumberhack.error("dimension error!")
    end
end

type VectorFieldNode{N} <: AbstractVectorField
    fields::Vector{AbstractVectorField}
    scaling::Function
    dim::Integer
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    typeof::DataType
    sample::Array{Complex{Float64},3}
    s::Complex{Float64}
    name::ASCIIString
    function VectorFieldNode{T<:AbstractVectorField}(f::Vector{T};scaling::Function  =  t->1.0+0.0im,name::ASCIIString="VectorFieldNode")
        @assert all(x->x.dim==N,f) "dimension error!"
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Complex{Float64},(3,4,4)),zero(Complex{Float64}),name)
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
    sample::Array{Float64,2}
    vf_sample::Array{Complex{Float64},3}
    s::Float64
    one_vf_flag::Bool
    name::ASCIIString
    function ScalarFieldNode{T<:Field}(f::Vector{T};scaling::Function =  t->1.0, name::ASCIIString="ScalarFieldNode")
        @assert all(x->x.dim==N,f) "dimension error!"
        flag = length(f) == 1 && typeof(f[1])<:AbstractVectorField
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Float64,(4,4)),zeros(Complex{Float64},(3,4,4)),zero(Float64),flag, name)
    end
end

FieldNode = Union{VectorFieldNode,ScalarFieldNode}
