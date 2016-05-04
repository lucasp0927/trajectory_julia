function set_geometry!{T<:FieldNode}(f::T)
    #memoize geo parameters
    geo = geometry(f)
    f.position = geo["pos"]
    f.size = geo["size"]
    f.res = geo["res"]    
    for ff in f.fields
        if typeof(ff) <: FieldNode
            set_geometry!(ff)
        end
    end
end

function geometry{T<:FieldNode}(f::T;lowest_res=true)
    #find the new optimal position and size
    #minimum(cat(2,[1,2,3],[2,3,4]),2)
    minpos::Vector{Float64} = vec(minimum(cat(2,map(x->geometry(x)["pos"],f.fields)...),2))
    maxpos::Vector{Float64} = vec(maximum(cat(2,map(x->geometry(x)["pos"].+geometry(x)["size"],f.fields)...),2))
    myminmax = (lowest_res?maximum:minimum)::Function
    res::Vector{Float64} = vec(myminmax(cat(2,map(x->geometry(x)["res"],f.fields)...),2))
    @assert all(x->x!=0,res) "zero resolution!"
    return Dict("pos"=>minpos, "size"=>(maxpos-minpos), "res"=>res)
end

geometry{T<:Union{VectorField,ScalarField}}(f::T) = Dict("pos"=>f.position,"size"=>f.size,"res"=>f.res)

