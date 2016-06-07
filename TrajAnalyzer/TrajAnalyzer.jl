module TrajAnalyzer
using Lumberjack
include("../fileio.jl")
global Trajs
type Trajectories
    traj::SharedArray{Float64}
    tspan::Vector{Float64}
    pos::Vector{Float64}
    siz::Vector{Float64}
    function Trajectories(result::Dict)
        Lumberjack.debug("construct Trajectories object")
        traj = copy_to_sharedarray!(result["traj"])
        tspan = copy(result["tspan"])
        pos = copy(result["pos"])
        siz = copy(result["siz"])
    end
end

function set_trajs_parallel!(result::Dict)
        @sync begin
        for p = 1:nprocs()
            @async remotecall_wait(p,set_trajs!,result)
        end
    end
end

function set_trajs!(result::Dict)
    global Trajs
    Trajs = Trajectories(result)
end
end
