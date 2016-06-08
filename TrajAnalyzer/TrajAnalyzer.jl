module TrajAnalyzer
using Fields
using Lumberjack
using HDF5
include("../fileio.jl")
include("../TrajSolver/polygon.jl")
include("TrajAnalyzer_set.jl")
include("TrajAnalyzer_output.jl")

global Trajs
type Trajectories
    traj::SharedArray{Float64}
    tspan::Vector{Float64}
    pos::Vector{Float64}
    siz::Vector{Float64}
    function Trajectories(result::Dict,traj_s::SharedArray{Float64})
        traj = traj_s
        tspan = copy(result["tspan"])
        pos = copy(result["pos"])
        siz = copy(result["siz"])
        new(traj,tspan,pos,siz)
    end
end

function calc_score(area)
    Lumberjack.debug("In TrajAnalyzer.calc_score()")
    pp = Polygon([promote(area...)...])
    score = @parallel (+) for i = 1:size(Trajs.traj,3)
        sum(map(j->pointInPolygon(pp,Trajs.traj[1:2,j,i])?1:0,1:size(Trajs.traj,2)))
    end
    Lumberjack.debug("score: $score")
    return score
end

end
