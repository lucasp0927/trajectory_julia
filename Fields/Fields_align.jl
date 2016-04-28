function align_field_tree!{T<:FieldNode}(f::T)
    unalign_geo = geometry(f)
    arr_sz = floor(Integer,collect(unalign_geo["size"])./collect(unalign_geo["res"]))
    new_sz = arr_sz.*collect(unalign_geo["res"])
    align_geo = unalign_geo
    align_geo["size"] = tuple(new_sz...)
    align_field!(f,align_geo["res"],align_geo["pos"])
    gc()
end

function align_field!{T<:FieldNode}(f::T,res::Tuple{Vararg{Float64}},pos::Tuple{Vararg{Float64}})
    @assert length(res) == length(pos) "dimension mismatch!"
    for ff in f.fields
        align_field!(ff,res,pos)
    end
end

function align_field!{T<:ComplexOrFloat,N}(f::ScalarField{T,N},res::Tuple{Vararg{Float64}},pos::Tuple{Vararg{Float64}})
    @assert length(res) == length(pos) "dimension mismatch!"
    unalign_geo = geometry(f)
    new_pos = ceil(collect(unalign_geo["pos"])./collect(res)).*collect(res)
    old_end = collect(unalign_geo["pos"]).+collect(unalign_geo["size"])
    new_end = floor(collect(old_end)./collect(res)).*collect(res)
    new_size = new_end.-new_pos
    align_geo = Dict("pos"=>tuple(new_pos...),"size"=>tuple(new_size...),"res"=>res)
    ##### interpolate field
    new_arr_size = round(Int64,new_size ./ collect(res))+one(Int64)
    new_field = Array(T,new_arr_size...)
    old_field = f.field
    old_field_itp = interpolate(old_field, BSpline(Cubic(Flat())), OnGrid())
    itp_field!(new_field,old_field_itp,unalign_geo,align_geo)
    setfield!(f,new_field,align_geo["pos"],align_geo["size"],scaling=f.scaling)
    @assert f.res == align_geo["res"] "resolution check failed!"
end

function align_field!{T<:ComplexOrFloat,N}(f::VectorField{T,N},res::Tuple{Vararg{Float64}},pos::Tuple{Vararg{Float64}})
    @assert length(res) == length(pos) "dimension mismatch!"
    unalign_geo = geometry(f)
    new_pos = ceil(collect(unalign_geo["pos"])./collect(res)).*collect(res)
    old_end = collect(unalign_geo["pos"]).+collect(unalign_geo["size"])
    new_end = floor(collect(old_end)./collect(res)).*collect(res)
    new_size = new_end.-new_pos
    align_geo = Dict("pos"=>tuple(new_pos...),"size"=>tuple(new_size...),"res"=>res)
    ##### interpolate field
    new_arr_size = round(Int64,new_size ./ collect(res))+one(Int64)
    new_field = Array(T,new_arr_size...,3)
    old_field = f.field
    #loop over three components
    for i = 1:3
        old_field_itp = interpolate(myslice(old_field,i), BSpline(Cubic(Flat())), OnGrid())
        itp_field!(myslice(new_field,i),old_field_itp,unalign_geo,align_geo)
    end
    setfield!(f,new_field,align_geo["pos"],align_geo["size"],scaling=f.scaling)
    @assert f.res == align_geo["res"] "resolution check failed!"
end

#TODO clean this up with meta programming?
function myslice{T<:ComplexOrFloat,N}(A::Array{T,N},i::Integer)
    if N == 3
        return slice(A,:,:,i)
    elseif N == 4
        return slice(A,:,:,:,i)
    end
end

@inbounds function transform_coordinate(apos::Vector{Float64},ares::Vector{Float64},uapos::Vector{Float64},uares::Vector{Float64},index::Vector{Int64})
    #2D
    #position at current index
    @devec    idx_pos = apos.+ ares.*(index-1)
    #calculate the old field index
    @devec    old_idx = ((idx_pos.-uapos)./uares)+1
    #test
    @devec    un_idx_pos = uapos.+ uares.*(old_idx-1)
    @assert   un_idx_pos == idx_pos "coordinate transform check"
    return tuple(old_idx...)
end

@generated function itp_field!{T<:Union{Array,SubArray}}(new_field::T,old_field_itp::Interpolations.BSplineInterpolation,unalign_geo,align_geo)
    N::Int64 = ndims(new_field)
    quote
        apos::Vector{Float64} = collect(align_geo["pos"])
        ares::Vector{Float64} = collect(align_geo["res"])
        uapos::Vector{Float64} = collect(unalign_geo["pos"])
        uares::Vector{Float64} = collect(unalign_geo["res"])                                
        @nloops $N i new_field begin
            old_idx = transform_coordinate(apos,ares,uapos,uares,collect((@ntuple $N i)))
            (@nref $N new_field i) = old_field_itp[old_idx...]
        end
    end
end

# @generated function itp_field!{T<:ComplexOrFloat,N}(new_field::SubArray{T,N},old_field_itp::Interpolations.BSplineInterpolation,unalign_geo,align_geo)
#     quote
#         @nloops $N i new_field begin
#             old_idx = transform_coordinate(unalign_geo,align_geo,(@ntuple $N i))
#             (@nref $N new_field i) = old_field_itp[old_idx...]
#         end
#     end
# end
