# functions
include("Fields_geometry.jl")
function initialize(res::Tuple{Integer,Integer}, sz::Tuple{Real,Real};dim = 2)
    global fields
    fields = ScalarFieldNode{dim}(Vector{Field}())
end

#####composition: for FieldNode object return composite field at time t, replace scaling with x->1.0
function composite{T<:Union{VectorField,ScalarField}}(f::T,t::Real) 
    T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end

function composite{N}(f::ScalarFieldNode{N})
    #TODO: implement ability to output both float64 and complex field
    #output scalar field, for now only output float64 field
    #remember scaling
    geo = geometry(f)
    res = geo["res"]
    pos = geo["pos"]
    sz = geo["size"]
    arr_sz = floor(Integer,collect(sz)./collect(res))
    new_sz = arr_sz.*collect(res)
    output = Array{Float64,N}(arr_sz...)
    ####handle output
    ff = map(x->composite(x),f.fields)
    ####
    ScalarField{Float64,N}(output,pos,new_sz;scaling = t->1.0)
end

# utility functions, for simple fields
function zero{T<:ComplexOrFloat,N}(::Type{ScalarField{T,N}},res::Tuple{Vararg{Integer}},pos::Tuple{Vararg{Real}},size::Tuple{Vararg{Real}};scaling::Function = t->1.0)
    return ScalarField{T,N}(zeros(T,res),pos,size,scaling=scaling)
end

function zero{T<:ComplexOrFloat,N}(::Type{VectorField{T,N}},res::Tuple{Vararg{Integer}},pos::Tuple{Vararg{Real}},size::Tuple{Vararg{Real}};scaling::Function = t->1.0)
    return VectorField{T,N}(zeros(T,(res...,3)),pos,size,scaling = scaling)
end

function func2field{T<:ComplexOrFloat}(::Type{ScalarField{T,2}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Function = t->1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    f = [func(x,y)::T for x in xx, y in yy]
    return ScalarField{T,2}(f::Array{T,2},pos,size,scaling=scaling)    
end

function func2field{T<:ComplexOrFloat}(::Type{ScalarField{T,3}},func::Function,res::Tuple{Integer,Integer,Integer},pos::Tuple{Real,Real,Real},size::Tuple{Real,Real,Real};scaling::Function = t->1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    zz = linspace(pos[3],pos[3]+size[3],res[3])
    f = [func(x,y,z)::T for x in xx, y in yy, z in zz]
    return ScalarField{T,3}(f::Array{T,3},pos,size,scaling=scaling)    
end

function func2field{T<:ComplexOrFloat}(::Type{VectorField{T,2}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Function = t->1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    f = zeros(T,(res...,3))
    for x in enumerate(xx), y in enumerate(yy)
        v = func(x[2],y[2])::Vector{T}
        f[x[1],y[1],:] = v
    end
    return VectorField{T,2}(f::Array{T,3},pos,size,scaling=scaling)    
end

function func2field{T<:ComplexOrFloat}(::Type{VectorField{T,3}},func::Function,res::Tuple{Integer,Integer,Integer},pos::Tuple{Real,Real,Real},size::Tuple{Real,Real,Real};scaling::Real = 1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    zz = linspace(pos[3],pos[3]+size[3],res[3])
    f = zeros(T,(res...,3))
    for x in enumerate(xx), y in enumerate(yy), z in enumerate(zz)
        v = func(x[2],y[2],z[2])::Vector{T}
        f[x[1],y[1],z[1],:] = v
    end
    return VectorField{T,3}(f::Array{T,4},pos,size,scaling=scaling)    
end
# file2field?