function geometry{T<:FieldNode}(f::T;lowest_res=true)
    #find the new optimal position and size
    #minimum(cat(2,[1,2,3],[2,3,4]),2)
    minpos::Vector{Float64} = vec(minimum(cat(2,map(x->collect(geometry(x)["pos"]),f.fields)...),2))
    maxpos::Vector{Float64} = vec(maximum(cat(2,map(x->(collect(geometry(x)["pos"]).+collect(geometry(x)["size"])),f.fields)...),2))
    myminmax = (lowest_res?minimum:maximum)::Function
    res::Vector{Float64} = vec(myminmax(cat(2,map(x->collect(geometry(x)["res"]),f.fields)...),2))
    if all(x->x!=0,res)
        return Dict("pos"=>tuple(minpos...), "size"=>tuple((maxpos-minpos)...), "res"=>tuple(res...))
    else
        error("zero resolution!")
    end
end

geometry{T<:Union{VectorField,ScalarField}}(f::T) = Dict("pos"=>f.position,"size"=>f.size,"res"=>f.res)
