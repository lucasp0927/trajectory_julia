# functions
include("Fields_geometry.jl")
include("Fields_typeof.jl")
include("Fields_composite.jl")
include("Fields_utility.jl")
# function initialize(res::Tuple{Integer,Integer}, sz::Tuple{Real,Real};dim = 2)
#     global fields
#     fields = ScalarFieldNode{dim}(Vector{Field}())
# end

function align_field_tree{T<:FieldNode}(f::T)
    unalign_geo = geometry(f)
    println("unalign_geo:",unalign_geo)
    arr_sz = floor(Integer,collect(unalign_geo["size"])./collect(unalign_geo["res"]))
    new_sz = arr_sz.*collect(unalign_geo["res"])
    align_geo = unalign_geo
    align_geo["size"] = tuple(new_sz...)
    println("align_geo",align_geo)
end

function align_field{T<:FieldNode}(f::T,res::Tuple{Vararg{Float64}},pos::Tuple{Vararg{Float64}})
    
end

function align_field{T<:ComplexOrFloat,N}(f::ScalarField{T,N})
end
