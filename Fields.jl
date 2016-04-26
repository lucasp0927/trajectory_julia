module Fields
using MAT
#TODO: resample field, for faster interpolation
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

# variables
export fields, U_total, resolution, size
global fields, U_total
global resolution, size
# functions
function initialize(res::Tuple{Integer,Integer}, sz::Tuple{Real,Real})
    global resolution, size, fields
    resolution = res
    size = sz
    U_total = zero(ScalarField{Float64,2},res,(0,0),size)
    fields = ScalarFieldNode(Vector{Field}())
end

function composite(f::ScalarFieldNode)
end

function composite(f::VectorFieldNode)
    #find the optimal new position and size, and return the composite vectorfield
    new_min_x = minimum(map((x)->getfield(x,:position)[1],f.fields))
    new_min_y = minimum(map((x)->getfield(x,:position)[2],f.fields))
    new_max_x = maximum(map((x)->(getfield(x,:position)[1]+getfield(x,:size)[1]),f.fields))
    new_max_y = maximum(map((x)->(getfield(x,:position)[2]+getfield(x,:size)[2]),f.fields))
    return new_min_x, new_min_y, new_max_x, new_max_y
end

composite(f::ScalarField) = f
composite(f::VectorField) = f

# utility functions, for simple fields
function zero{T<:ComplexOrFloat,N}(::Type{ScalarField{T,N}},res::Tuple{Vararg{Integer}},pos::Tuple{Vararg{Real}},size::Tuple{Vararg{Real}};scaling::Real = 1.0)
    return ScalarField{T,N}(zeros(T,res),pos,size,scaling=scaling)
end

function zero{T<:ComplexOrFloat,N}(::Type{VectorField{T,N}},res::Tuple{Vararg{Integer}},pos::Tuple{Vararg{Real}},size::Tuple{Vararg{Real}};scaling::Real = 1.0)
    return VectorField{T,N}(zeros(T,(res...,3)),pos,size,scaling = scaling)
end

function func2field{T<:ComplexOrFloat}(::Type{ScalarField{T}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Real = 1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    f = [func(x,y)::T for x in xx, y in yy]
    return ScalarField{T}(f::Array{T,2},pos,size,scaling)    
end

function func2field{T<:ComplexOrFloat}(::Type{VectorField{T}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Real = 1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    f = zeros(T,(res...,3))
    for x in enumerate(xx), y in enumerate(yy)
        v = func(x[2],y[2])::Vector{T}
        f[x[1],y[1],:] = v
    end
    return VectorField{T}(f::Array{T,3},pos,size,scaling)    
end
# file2field?
end
