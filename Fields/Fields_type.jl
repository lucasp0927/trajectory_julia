# type definitions
using SharedArrays
export ComplexOrFloat
export Field, AbstractVectorField, AbstractScalarField, VectorField, ScalarField, VectorFieldNode, ScalarFieldNode, FieldNode
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract type Field end
abstract type AbstractVectorField <: Field end
abstract type AbstractScalarField <: Field end
#TODO: right now sample is only for 2d, also the ScalarFieldNode.sample datatype is Float64

mutable struct ScalarField{T <: ComplexOrFloat,N} <: AbstractScalarField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    scaling_expr::Expr
    dim::Integer
    sample::Array{T,N}
#    sample_views::SubArray{T,N}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Float64
    name::String
    function ScalarField{T,N}(f::SharedArray{T,N},pos::Vector{Float64},sz::Vector{Float64};scaling_expr::Expr = parse("t->1.0"), name::String = "ScalarField")  where {T <: ComplexOrFloat,N}
        res = sz./(collect(size(f))[1:N].-1)
        @assert all(x->x!=0,res) "zero resolution!"
        @assert length(pos)==length(sz)==N==ndims(f)
#        new(f,pos,sz,res,eval(scaling_expr),scaling_expr,N,zeros(T,repmat([4],N)...),repmat([0.0],N),repmat([0],N*2),zero(Float64),ascii(name))
        new(f,pos,sz,res,eval(scaling_expr),scaling_expr,N,zeros(T,repeat([4],N)...),repeat([0.0],N),repeat([0],N*2),zero(Float64),ascii(name))
    end
end

mutable struct VectorField{T <: ComplexOrFloat, N} <: AbstractVectorField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    scaling_expr::Expr
    dim::Integer
    sample::Array{T}
    sample_views::SubArray{T}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Complex{Float64}
    name::String
    function VectorField{T,N}(f::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling_expr::Expr =  parse("t->1.0"), name::String = "VectorField") where {T <: ComplexOrFloat, N}
        res = sz./(collect(size(f))[2:N+1].-1)
        @assert all(x->x!=0,res) "zero resolution!"
        @assert length(pos)==length(sz)==N==ndims(f)-1
        s = zeros(T,[3,repeat([4],N)...]...)
        sv = @views s[:]
        new(f,pos,sz,res,eval(scaling_expr),scaling_expr,N,s,sv,repeat([0.0],N),repeat([0],N*2),zero(Complex{Float64}),ascii(name))
    end
end

mutable struct VectorFieldNode{N} <: AbstractVectorField
    fields::Vector{AbstractVectorField}
    scaling::Function
    scaling_expr::Expr
    dim::Integer
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    typeof::DataType
#    sample::Array{Complex{Float64},N+1}
    sample::Array{Complex{Float64}}
    s::Complex{Float64}
    name::String
    function VectorFieldNode{N}(f::Vector{T};scaling_expr::Expr  = parse("t->1.0+0.0im"),name::String="VectorFieldNode") where {T<:AbstractVectorField, N}
        @assert all(x->x.dim==N,f) "dimension error!"
        new(f,eval(scaling_expr),scaling_expr,N,[],[],[],Complex{Float64},zeros(Complex{Float64},[3,repeat([4],N)...]...),zero(Complex{Float64}),ascii(name))
    end
end

mutable struct ScalarFieldNode{N} <: AbstractScalarField
    fields::Vector{Field}
    scaling::Function
    scaling_expr::Expr
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
    name::String
    #TODO: restrict type
    function ScalarFieldNode{N}(f::Vector{T};scaling_expr::Expr = parse("t->1.0"), name::String="ScalarFieldNode") where {T<:Field, N}
        @assert all(x->x.dim==N,f) "dimension error!"
        flag = (length(f) == 1) && typeof(f[1])<:AbstractVectorField
        new(f,eval(scaling_expr),scaling_expr,N,[],[],[],Complex{Float64},zeros(Float64,repeat([4],N)...),zeros(Complex{Float64},[3,repeat([4],N)...]...),zero(Float64),flag, ascii(name))
    end
end

FieldNode = Union{VectorFieldNode,ScalarFieldNode}
FieldNode2D = Union{VectorFieldNode{2},ScalarFieldNode{2}}
FieldNode3D = Union{VectorFieldNode{3},ScalarFieldNode{3}}
ComplexField2D = Union{VectorField{Complex{Float64},2},ScalarField{Complex{Float64},2}}
ComplexField3D = Union{VectorField{Complex{Float64},3},ScalarField{Complex{Float64},3}}
FloatField2D = Union{VectorField{Float64,2},ScalarField{Float64,2}}
FloatField3D = Union{VectorField{Float64,3},ScalarField{Float64,3}}
