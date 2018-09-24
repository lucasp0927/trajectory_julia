module TrajAnalyzer
using Fields
using Logging
#using PyCall
using ProgressMeter
#using MAT
using HDF5
include("../fileio.jl")
include("TrajAnalyzer_trajectories.jl")
include("../TrajSolver/polygon.jl")
include("TrajAnalyzer_set.jl")
include("TrajAnalyzer_output.jl")
include("TrajAnalyzer_spectrum.jl")
include("TrajAnalyzer_ngamma1d.jl")
include("TrajAnalyzer_trajtype.jl")

global Trajs, Probe, ForceFields, TA_Config
global avg_atom_num,lattice_width,lattice_unit,k_ratio,gamma_1d,gamma_prime
global range_i, range_j

#for TrajAnalyzer_output to determin using ffmpeg or avconv.
@pyimport platform
function calc_score(area)
    pp = Polygon([promote(area...)...])
    score = @parallel (+) for i = 1:size(Trajs.traj,3)
        sum(map(j->pointInPolygon(pp,Trajs.traj[1:2,j,i])?1:0,1:size(Trajs.traj,2)))
    end
    @debug "score: $score"
    return score
end
end
