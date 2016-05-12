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
    function VectorField(f::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling =  t->1.0)
        res = sz./(collect(size(f))[2:N+1]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,res,scaling,N,zeros(T,(3,4,4)),[0.0,0.0],[0,0,0,0],zero(Complex{Float64})):error("dimension error!")
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
    function ScalarField(f::SharedArray{T,N},pos::Vector{Float64},sz::Vector{Float64};scaling =  t->1.0)
        res = sz./(collect(size(f))[1:N]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,res,scaling,N,zeros(T,(4,4)),[0.0,0.0],[0,0,0,0],zero(Float64)):error("dimension error!")
    end
end

function setfield!{T<:ComplexOrFloat,N}(f::ScalarField{T,N},A::SharedArray{T},pos::Vector{Float64},sz::Array{Float64};scaling = t->1.0)
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

function setfield!{T<:ComplexOrFloat,N}(f::VectorField{T,N},A::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling = t->1.0)
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
    sample::Array{Complex{Float64},3}
    s::Complex{Float64}
    function VectorFieldNode{T<:AbstractVectorField}(f::Vector{T};scaling  =  t->1.0)
        @assert all(x->x.dim==N,f) "dimension error!"
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Complex{Float64},(3,4,4)),zero(Complex{Float64}))
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
    function ScalarFieldNode{T<:Field}(f::Vector{T};scaling =  t->1.0)
        @assert all(x->x.dim==N,f) "dimension error!"
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Float64,(4,4)),zeros(Complex{Float64},(3,4,4)),zero(Float64))
    end
end

FieldNode = Union{VectorFieldNode,ScalarFieldNode}

function copyfield{T<:ComplexOrFloat,N}(f::ScalarField{T,N})
    return ScalarField{T,N}(f.field,f.position,f.size,scaling=f.scaling)
end

function copyfield{T<:ComplexOrFloat,N}(f::VectorField{T,N})
    return VectorField{T,N}(f.field,f.position,f.size,scaling=f.scaling)
end

function copyfield{N}(f::ScalarFieldNode{N})
    fields = map(copyfield,f.fields)
    fn = ScalarFieldNode{N}(fields,scaling = f.scaling)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end

function copyfield{N}(f::VectorFieldNode{N})
    fields = map(copyfield,f.fields)
    fn = VectorFieldNode{N}(fields,scaling = f.scaling)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end
