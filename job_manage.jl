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
        println("change scaling of field $field_name to ",s)
        s_exp = eval(parse(s))
        Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        Fields.init_parallel!(sfn)
        println("start calculation...")
        if calc_traj_flag
            traj = calculate_traj()
            println("save results...")
            tspan = TrajSolver.get_tspan()
            result = Dict(
                          "traj"=>traj,
                          "tspan" =>tspan,
                          "pos"=>sfn.position,
                          "siz"=>sfn.size
                          )
            matwrite(output_file*string(i)*".mat",result)
        else
            vars = matread(output_file*string(i)*".mat")
            traj = vars["traj"]
            tspan = vars["tspan"]
        end
        for (k,v) in config["score"]
            println("calculating score for area $k...")
            @time (score[ascii(k)])[i] = calc_score(traj,v)
        end
        @time flux = calc_flux(traj,tspan,config["flux"],output_file*string(i)*"_flux.mat")
        if movie_flag
            @time output_movie_traj(config["movie-output"],output_file*string(i)*"_traj.mp4",traj,tspan)
        end
    end
    #save score and plot score
    matwrite(output_file*"score.mat",score)
    # pyimport("matplotlib")[:use]("Agg")
    # @pyimport matplotlib.pyplot as plt
    # prange = config["plot-range"]
    # xx = collect(linspace(prange...,range))
    # plt.plot(xx,score)
    # plt.savefig(output_file*"score.png")
    # plt.clf()
end


function calculate_traj()
    @time @sync begin
        for p = 2:nprocs()
            @async remotecall_wait(p,TrajSolver.solve_traj)
        end
    end
    temp = cell(nworkers())
    for p = 2:nprocs()
        temp[p-1] = remotecall_fetch(p,TrajSolver.get_result)
    end
    traj = cat(3,temp...)
    return traj
end

function calc_score(traj,area)
    @everywhere include("./TrajSolver/polygon.jl")
    pp = Polygon([promote(area...)...])
    traj_s = copy_to_sharedarray!(traj)
    score = @parallel (+) for i = 1:size(traj_s,3)
        sum(map(j->pointInPolygon(pp,traj_s[1:2,j,i])?1:0,1:size(traj_s,2)))
        # for j = 1:size(traj_s,2)
        #     pointInPolygon(pp,traj_s[1:2,j,i])?1:0
        # end
    end    
#    score = sum([pointInPolygon(pp,traj[1:2,j,i]) for i = 1:size(traj,3),j = 1:size(traj,2)])
    # score = 0
    # for i = 1:size(traj,3),j = 1:size(traj,2)
    #     if pointInPolygon(pp,traj[1:2,j,i])
    #         score+=1
    #     end
    # end
    println("score:",score)
    return score
end
