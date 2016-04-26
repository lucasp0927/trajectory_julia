# functions
function initialize(res::Tuple{Integer,Integer}, sz::Tuple{Real,Real};dim = 2)
    global resolution, size, fields
    resolution = res
    size = sz
    U_total = zero(ScalarField{Float64,2},res,(0,0),size)
    fields = ScalarFieldNode{dim}(Vector{Field}())
end

function geometry{T<:FieldNode}(f::T)
    #find the new optimal position and size
    #minimum(cat(2,[1,2,3],[2,3,4]),2)
    minpos::Vector{Float64} = vec(minimum(cat(2,map(x->collect(geometry(x)[1]),f.fields)...),2))
    maxpos::Vector{Float64} = vec(maximum(cat(2,map(x->(collect(geometry(x)[1]).+collect(geometry(x)[2])),f.fields)...),2))
    return tuple(minpos...), tuple((maxpos-minpos)...)
end

geometry{T<:Union{VectorField,ScalarField}}(f::T) = f.position,f.size

# utility functions, for simple fields
function zero{T<:ComplexOrFloat,N}(::Type{ScalarField{T,N}},res::Tuple{Vararg{Integer}},pos::Tuple{Vararg{Real}},size::Tuple{Vararg{Real}};scaling::Real = 1.0)
    return ScalarField{T,N}(zeros(T,res),pos,size,scaling=scaling)
end

function zero{T<:ComplexOrFloat,N}(::Type{VectorField{T,N}},res::Tuple{Vararg{Integer}},pos::Tuple{Vararg{Real}},size::Tuple{Vararg{Real}};scaling::Real = 1.0)
    return VectorField{T,N}(zeros(T,(res...,3)),pos,size,scaling = scaling)
end

function func2field{T<:ComplexOrFloat}(::Type{ScalarField{T,2}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Real = 1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    f = [func(x,y)::T for x in xx, y in yy]
    return ScalarField{T,2}(f::Array{T,2},pos,size,scaling=scaling)    
end

function func2field{T<:ComplexOrFloat}(::Type{ScalarField{T,3}},func::Function,res::Tuple{Integer,Integer,Integer},pos::Tuple{Real,Real,Real},size::Tuple{Real,Real,Real};scaling::Real = 1.0)
    xx = linspace(pos[1],pos[1]+size[1],res[1])
    yy = linspace(pos[2],pos[2]+size[2],res[2])
    zz = linspace(pos[3],pos[3]+size[3],res[3])
    f = [func(x,y,z)::T for x in xx, y in yy, z in zz]
    return ScalarField{T,3}(f::Array{T,3},pos,size,scaling=scaling)    
end


function func2field{T<:ComplexOrFloat}(::Type{VectorField{T,2}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Real = 1.0)
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