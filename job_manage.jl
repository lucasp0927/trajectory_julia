using PyCall
include("output.jl")
include("flux.jl")
@everywhere using HDF5
function single_scan_scaling(config::Dict,sfn::ScalarFieldNode,output_file,calc_traj_flag::Bool,movie_flag::Bool)
    range = config["range"]
    field_name = config["field"]
    scaling = config["scaling"]
    score = [ascii(k)=>zeros(Int64,range) for k in keys(config["score"])]
    for i = 1:range
        s = replace(scaling,"@i",float(i))
        Lumberjack.info("change scaling of field $field_name to ",s)
        s_exp = eval(parse(s))
        Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        Fields.init_parallel!(sfn)
        if calc_traj_flag
            result = calculate_traj()
            TrajAnalyzer.set_trajs_parallel!(result)
            Lumberjack.info("save results...")
            matwrite(output_file*string(i)*".mat",result)
            traj = result["traj"]
            tspan = result["tspan"]
        else
            vars = matread(output_file*string(i)*".mat")
            traj = vars["traj"]
            tspan = vars["tspan"]
        end
        for (k,v) in config["score"]
            Lumberjack.info("calculating score for area $k...")
            @time (score[ascii(k)])[i] = calc_score(traj,v)
        end
        @time flux = calc_flux(traj,tspan,config["flux"],output_file*string(i)*"_flux.mat")
        if movie_flag
            @time output_movie_traj(config["movie-output"],output_file*string(i)*"_traj.mp4",traj,tspan)
        end
    end
    matwrite(output_file*"score.mat",score)
end

function calc_score(traj,area)
    @everywhere include("./TrajSolver/polygon.jl")
    pp = Polygon([promote(area...)...])
    traj_s = copy_to_sharedarray!(traj)
    score = @parallel (+) for i = 1:size(traj_s,3)
        sum(map(j->pointInPolygon(pp,traj_s[1:2,j,i])?1:0,1:size(traj_s,2)))
    end
    Lumberjack.debug("score: $score")
    return score
end
