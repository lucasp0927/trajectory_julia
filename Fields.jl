module Fields
using MAT
# type definitions
export ComplexOrFloat
export Field, VectorField, ScalarField
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract Field
type VectorField{T <: ComplexOrFloat} <: Field
    field::Array{T,3}
    position::Tuple{Float64,Float64}
    size::Tuple{Float64,Float64}
    scaling::Float64
end

type ScalarField{T <: ComplexOrFloat} <: Field
    field::Array{T,2}
    position::Tuple{Float64,Float64}
    size::Tuple{Float64,Float64}
    scaling::Float64
end

abstract FieldNode <: Field
type VectorFieldNode <: FieldNode
    fields::Vector{Union{VectorFieldNode, VectorField}}
end

type ScalarFieldNode <: FieldNode
    fields::Vector{Union{ScalarFieldNode, ScalarField}}
end

# variables
export fields, U_total, resolution, size
global fields, U_total
global resolution, size
# functions
function initialize(res::Tuple{Integer,Integer}, sz::Tuple{Real,Real})
    global resolution, size
    resolution = res
    size = sz
    U_total = zero(ScalarField{Float64},res,(0,0),size)
end

# utility functions, for simple fields
function zero{T<:ComplexOrFloat}(::Type{ScalarField{T}},res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Real = 1.0)
    return ScalarField{T}(zeros(T,res),pos,size,scaling)
end

function zero{T<:ComplexOrFloat}(::Type{VectorField{T}},res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Real = 1.0)
    return VectorField{T}(zeros(T,(res...,3)),pos,size,scaling)
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
