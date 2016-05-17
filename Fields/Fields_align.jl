function align_field_tree!{T<:FieldNode}(f::T)
    unalign_geo = geometry(f)
    arr_sz = floor(Integer,unalign_geo["size"]./unalign_geo["res"])
    new_sz = arr_sz.*unalign_geo["res"]
    align_geo = unalign_geo
    align_geo["size"] = new_sz
    align_field!(f,align_geo["res"],align_geo["pos"])
    gc()
end

function align_field!{T<:FieldNode}(f::T,res::Vector{Float64},pos::Vector{Float64})
    @assert length(res) == length(pos) "dimension mismatch!"
    for ff in f.fields
        align_field!(ff,res,pos)
    end
end

function align_field!{T<:ComplexOrFloat,N}(f::ScalarField{T,N},res::Vector{Float64},pos::Vector{Float64})
    @assert length(res) == length(pos) "dimension mismatch!"
    unalign_geo = geometry(f)
    new_pos = ceil(unalign_geo["pos"]./res).*res
    old_end = unalign_geo["pos"].+unalign_geo["size"]
    new_end = floor(old_end./res).*res
    new_size = new_end.-new_pos
    align_geo = Dict("pos"=>new_pos,"size"=>new_size,"res"=>res)
    ##### interpolate field
    new_arr_size = round(Int64,new_size ./ res)+one(Int64)
    old_field_itp = interpolate(f.field, BSpline(Cubic(Flat())), OnGrid())
    new_field = SharedArray(T,new_arr_size...)
    itp_field!(new_field,old_field_itp,unalign_geo,align_geo)
    @assert collect(size(new_field)) == new_arr_size
    setfield!(f,new_field,align_geo["pos"],align_geo["size"],scaling=f.scaling)
    #check
    res_delta = mean((f.res.-align_geo["res"])./f.res)
    @assert res_delta<=1e-10 "resolution check failed!"
end

function align_field!{T<:ComplexOrFloat,N}(f::VectorField{T,N},res::Vector{Float64},pos::Vector{Float64})
    @assert length(res) == length(pos) "dimension mismatch!"
    unalign_geo = geometry(f)
    new_pos = ceil(unalign_geo["pos"]./res).*res
    old_end = unalign_geo["pos"].+unalign_geo["size"]
    new_end = floor(old_end./res).*res
    new_size = new_end.-new_pos
    align_geo = Dict("pos"=>new_pos,"size"=>new_size,"res"=>res)
    ##### interpolate field
    new_arr_size = round(Int64,new_size ./ res)+one(Int64)
    #loop over three components
    new_field = SharedArray(T,3,new_arr_size...)
    for i = 1:3
        old_field_itp = interpolate(myslice(f.field,i), BSpline(Cubic(Flat())), OnGrid())
        itp_field!(myslice(new_field,i),old_field_itp,unalign_geo,align_geo)
    end
    @assert collect(size(new_field)[2:end]) == new_arr_size
    setfield!(f,new_field,align_geo["pos"],align_geo["size"],scaling=f.scaling)
    #check
    res_delta = mean((f.res.-align_geo["res"])./f.res)
    @assert res_delta<=1e-10 "resolution check failed!"
end

#TODO clean this up with meta programming?
@generated function myslice{T<:ComplexOrFloat,N}(A::Array{T,N},i::Integer)
    ex_str = "slice(A,i"*repeat(",:",N-1)*")"
    parse(ex_str)
end

@generated function myslice{T<:ComplexOrFloat,N}(A::SharedArray{T,N},i::Integer)
    ex_str = "slice(A,i"*repeat(",:",N-1)*")"
    parse(ex_str)
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
    return old_idx
end

@generated function itp_field!{T<:Union{Array,SubArray,SharedArray}}(new_field::T,old_field_itp::Interpolations.BSplineInterpolation,unalign_geo,align_geo)
    N::Int64 = ndims(new_field)
    quote
        apos::Vector{Float64} = align_geo["pos"]
        ares::Vector{Float64} = align_geo["res"]
        uapos::Vector{Float64} = unalign_geo["pos"]
        uares::Vector{Float64} = unalign_geo["res"]
        start_idx = transform_coordinate(apos,ares,uapos,uares,ones(Int64,$N))
        end_idx = transform_coordinate(apos,ares,uapos,uares,collect(size(new_field)))
        @nexprs $N j-> x_j = linspace(start_idx[j],end_idx[j],size(new_field,j))
        @nloops $N i new_field begin
            (@nref $N new_field i) = (@nref $N old_field_itp j->x_j[i_j])
        end
    end
end
