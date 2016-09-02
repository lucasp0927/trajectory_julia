type Trajectories
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
        new(traj,atom_num,tspan,t_div,pos,siz)
    end
end
import Base.getindex
@inbounds function getindex(tr::Trajectories,t::Float64,traj_id::Int64)
    #Linear interpolation
    #TODO: higher order interpolation
#    @assert t>=tr.tspan[1] && t<=tr.tspan[end]
    t_idx = searchsortedlast(tr.tspan,t)
    t_div = tr.t_div::Float64
    r = (t-tr.tspan[t_idx])/t_div
    result = Array(Float64,size(tr.traj,1))
    @simd for i = 1:size(tr.traj,1)
        result[i] = tr.traj[i,t_idx,traj_id]*(1.0-r)
        result[i] += tr.traj[i,t_idx+1,traj_id]*r
    end
    return result
end

@inbounds function getindex{T<:Range}(tr::Trajectories,t::Float64,traj_id::T)
    #Linear interpolation
    #TODO: higher order interpolation
#    @assert t>=tr.tspan[1] && t<=tr.tspan[end]
    t_idx = searchsortedlast(tr.tspan,t)
    t_div = tr.t_div::Float64
    r = (t-tr.tspan[t_idx])/t_div
    result = Array(Float64,size(tr.traj,1),length(traj_id))
    @simd for i = 1:size(tr.traj,1)
        result[i,:] = tr.traj[i,t_idx,traj_id].*(1.0-r)
        result[i,:] += tr.traj[i,t_idx+1,traj_id].*r
    end
    return result
end

@inbounds function getindex{T<:Colon}(tr::Trajectories,t::Float64,traj_id::T)
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


function getindex{T<:Range}(tr::Trajectories,t_range::T,traj_id)
    #Linear interpolation
    #TODO: higher order interpolation
    reduce((x,y)->cat(2,x,y),[tr[t,traj_id] for t in t_range])
end
