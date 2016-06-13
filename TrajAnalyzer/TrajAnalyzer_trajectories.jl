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
function getindex(tr::Trajectories,t::Float64,traj_id)
    #Linear interpolation
    #TODO: higher order interpolation
    t_idx = indmin(abs(tr.tspan-t))::Int64
    t_div = tr.t_div::Float64
    if t>=tr.tspan[t_idx]
        dt = t-tr.tspan[t_idx]
        return tr.traj[:,t_idx,traj_id]*(1.0-dt/t_div) .+ tr.traj[:,t_idx+1,traj_id]*dt/t_div
    else
        dt = tr.tspan[t_idx]-t
        return tr.traj[:,t_idx,traj_id]*(1.0-dt/t_div) .+ tr.traj[:,t_idx-1,traj_id]*dt/t_div
    end
end

function getindex{T<:Range}(tr::Trajectories,t_range::T,traj_id)
    #Linear interpolation
    #TODO: higher order interpolation
    reduce((x,y)->cat(2,x,y),[tr[t,traj_id] for t in t_range])
end
