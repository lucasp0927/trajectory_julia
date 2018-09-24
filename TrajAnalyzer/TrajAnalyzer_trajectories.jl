using SharedArrays
mutable struct Trajectories
    traj::SharedArray{Float64}  #traj[1,:,:] =>x coordinate, traj[2,:,:] => y coordinate
    atom_num::Int64
    tspan::Vector{Float64}
    t_div::Float64
    pos::Vector{Float64}
    siz::Vector{Float64}
    function Trajectories(result::Dict,traj_s::SharedArray{Float64})
        traj = traj_s
        atom_num = size(traj,3)
        tspan = copy(result["tspan"])
        t_div = mean(diff(tspan))
        pos = copy(result["pos"])
        siz = copy(result["siz"])
        new(traj,atom_num,vec(tspan),t_div,vec(pos),vec(siz))
    end
end

import Base.getindex
@inbounds function getindex(tr::Trajectories,t::Float64,traj_id::Int64)
    #Linear interpolation
    #TODO: higher order interpolation
#    @assert t>=tr.tspan[1] && t<=tr.tspan[end]
    t_idx = searchsortedlast(tr.tspan,t)
    @assert t_idx >=1 && t_idx <= length(tr.tspan)
    t_div = tr.t_div::Float64
    result = Array{Float64}(size(tr.traj,1))
    if t_idx != length(tr.tspan)
        r = (t-tr.tspan[t_idx])/t_div
        result[:] = tr.traj[:,t_idx,traj_id].*(1.0-r).+tr.traj[:,t_idx+1,traj_id].*r
    else
        result[:] = tr.traj[:,t_idx,traj_id]
    end
    return result
end

@inbounds function getindex(tr::Trajectories,t::Float64,traj_id::T) where {T<:Range}
    #Linear interpolation
    #TODO: higher order interpolation
#    @assert t>=tr.tspan[1] && t<=tr.tspan[end]
    t_idx = searchsortedlast(tr.tspan,t)
    @assert t_idx >=1 && t_idx <= length(tr.tspan)
    t_div = tr.t_div::Float64
    result = Array{Float64}(size(tr.traj,1),length(traj_id))
    if t_idx != length(tr.tspan)
        r = (t-tr.tspan[t_idx])/t_div
        result .= tr.traj[:,t_idx,traj_id].*(1.0-r)
        result .+= tr.traj[:,t_idx+1,traj_id].*r
#        result[:,:] = tr.traj[:,t_idx,traj_id].*(1.0-r)
#        result[:,:] += tr.traj[:,t_idx+1,traj_id].*r
        # @simd for i = 1:size(tr.traj,1)
        #     result[i,:] = tr.traj[i,t_idx,traj_id].*(1.0-r)
        #     result[i,:] += tr.traj[i,t_idx+1,traj_id].*r
        # end
    else
        #result[:,:] = tr.traj[:,t_idx,traj_id]
        result .= tr.traj[:,t_idx,traj_id]
        # @simd for i = 1:size(tr.traj,1)
        #     result[i,:] = tr.traj[i,t_idx,traj_id]
        # end
    end
    return result
end

function getindex(tr::Trajectories,t::Float64,traj_id::T) where {T <: Colon}
    getindex(tr,t,range(1,size(tr.traj,3)))
    #Linear interpolation
    #TODO: higher order interpolation
    #    @assert t>=tr.tspan[1] && t<=tr.tspan[end]
    #=
    t_idx = searchsortedlast(tr.tspan,t)
    t_div = tr.t_div::Float64
    r = (t-tr.tspan[t_idx])/t_div
    result = Array(Float64,size(tr.traj,1),size(tr.traj,3))
    @simd for i = 1:size(tr.traj,1)
        result[i,:] = tr.traj[i,t_idx,:].*(1.0-r)
        result[i,:] += tr.traj[i,t_idx+1,:].*r
    end
    return result
=#
end

function getindex(tr::Trajectories,t_range::T,traj_id) where {T <: Range}
    #Linear interpolation
    #TODO: higher order interpolation
    reduce((x,y)->cat(2,x,y),[tr[t,traj_id] for t in t_range])
end

function removenan(traj_snapshot::Array{Float64,2})
    x1 = traj_snapshot[1,:]
    x2 = traj_snapshot[2,:]
    x3 = traj_snapshot[3,:]
    x4 = traj_snapshot[4,:]
    filter!(x->~isnan(x),x1)
    filter!(x->~isnan(x),x2)
    filter!(x->~isnan(x),x3)
    filter!(x->~isnan(x),x4)
    traj_snapshot = cat(1,x1',x2',x3',x4')
    return traj_snapshot
end
