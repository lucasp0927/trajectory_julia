module TrajAnalyzer
using Fields
using Lumberjack
using HDF5
include("../fileio.jl")
include("TrajAnalyzer_trajectories.jl")
include("../TrajSolver/polygon.jl")
include("TrajAnalyzer_set.jl")
include("TrajAnalyzer_output.jl")

#include("TrajAnalyzer_spectrum.jl")
global Trajs, Probe, ForceFields, TA_Config
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
