# type definitions
export ComplexOrFloat
export Field, AbstractVectorField, AbstractScalarField, VectorField, ScalarField, VectorFieldNode, ScalarFieldNode, FieldNode
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract Field
abstract AbstractVectorField <: Field
abstract AbstractScalarField <: Field
#TODO: right now sample is only for 2d, also the ScalarFieldNode.sample datatype is Float64

type ScalarField{T <: ComplexOrFloat,N} <: AbstractScalarField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
    sample::Array{T,N}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Float64
    name::ASCIIString
    function ScalarField(f::SharedArray{T,N},pos::Vector{Float64},sz::Vector{Float64};scaling::Function =  t->1.0, name::ASCIIString = "ScalarField")
        res = sz./(collect(size(f))[1:N]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,res,scaling,N,zeros(T,(4,4)),repmat([0.0],N),repmat([0],N*2),zero(Float64),ascii(name)):Lumberhack.error("dimension error!")
    end
end

type VectorField{T <: ComplexOrFloat, N} <: AbstractVectorField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
#    sample::Array{T,N+1}
    sample::Array{T}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Complex{Float64}
    name::ASCIIString
    function VectorField(f::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling::Function =  t->1.0, name::ASCIIString = "VectorField")
        res = sz./(collect(size(f))[2:N+1]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,res,scaling,N,zeros(T,(3,4,4)),repmat([0.0],N),repmat([0],N*2),zero(Complex{Float64}),ascii(name)):Lumberjack.error("dimension error!")
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
#    sample::Array{Complex{Float64},N+1}
    sample::Array{Complex{Float64}}
    s::Complex{Float64}
    name::ASCIIString
    function VectorFieldNode{T<:AbstractVectorField}(f::Vector{T};scaling::Function  =  t->1.0+0.0im,name::ASCIIString="VectorFieldNode")
        @assert all(x->x.dim==N,f) "dimension error!"
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Complex{Float64},(3,4,4)),zero(Complex{Float64}),ascii(name))
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
    sample::Array{Float64,N}
#    vf_sample::Array{Complex{Float64},N+1}
    vf_sample::Array{Complex{Float64}}
    s::Float64
    one_vf_flag::Bool
    name::ASCIIString
    function ScalarFieldNode{T<:Field}(f::Vector{T};scaling::Function =  t->1.0, name::ASCIIString="ScalarFieldNode")
        @assert all(x->x.dim==N,f) "dimension error!"
        flag = (length(f) == 1) && typeof(f[1])<:AbstractVectorField
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Float64,(4,4)),zeros(Complex{Float64},(3,4,4)),zero(Float64),flag, ascii(name))
    end
end

FieldNode = Union{VectorFieldNode,ScalarFieldNode}
FieldNode2D = Union{VectorFieldNode{2},ScalarFieldNode{2}}
FieldNode3D = Union{VectorFieldNode{3},ScalarFieldNode{3}}
ComplexField2D = Union{VectorField{Complex{Float64},2},ScalarField{Complex{Float64},2}}
ComplexField3D = Union{VectorField{Complex{Float64},3},ScalarField{Complex{Float64},3}}
FloatField2D = Union{VectorField{Float64,2},ScalarField{Float64,2}}
FloatField3D = Union{VectorField{Float64,3},ScalarField{Float64,3}}
