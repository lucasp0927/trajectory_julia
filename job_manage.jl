include("flux.jl")
#include("benchmark.jl")
using Optim
using PyCall

@pyimport scipy.optimize as opt
function job_inner_loop(config,sfn,probe_sfn,input_prefix,output_prefix,flags,idx::Vector)
    if config["type"] == "single-scan-scaling"
        @assert length(idx) == 1
    elseif config["type"] == "double-scan-scaling"
        @assert length(idx) == 2
    end
    sim_type = config["movie-output"]["simulation-type"]
    function prepare_var()
        if flags["calc_traj_flag"]
            result = calculate_traj()
            @info "save results..."
            dicttoh5(output_prefix*".h5",result)
            traj = result["traj"]
            tspan = result["tspan"]
        elseif flags["need_traj_flag"]
            result = h5todict(input_prefix*".h5");
            traj = result["traj"]
            tspan = result["tspan"]
        end
        return result,traj,tspan
    end
    result,traj,tspan = prepare_var();

    if flags["need_traj_flag"]
        @info "Initialize TrajAnalyzer..."
        #probe_sfn = Fields.buildAndAlign(config["probe"]["field"],0,name=ascii([k for k in keys(config["probe"])][1]))
        TrajAnalyzer.init_parallel!(result,probe_sfn,config)
        #crashed_num = length(TrajAnalyzer.traj_iscrashed())
        #gap_num = length(TrajAnalyzer.traj_ingap(false))
        #@info "$crashed_num trajectories crashed."
        #@info "$gap_num trajectories in gap."
    end
    result = nothing
    if flags["movie_flag"]
        @info "Outputing Movie..."
        if sim_type == "2D"
            movie_range = [promote(config["movie-output"]["range"]...)...]
            @time TrajAnalyzer.output_movie_traj_2d(config["movie-output"],output_prefix*"_traj.mp4")
        elseif sim_type == "3D"
            movie_range = [promote(config["movie-output"]["range"]...)...]
            @time TrajAnalyzer.output_movie_traj_3d(config["movie-output"],output_prefix*"_traj.mp4")
        else
            @error "unrecognized simulation type: "* sim_type
        end
    end
    if flags["spectrum_flag"]
        if sim_type == "2D"
            for gm_name in config["spectrum"]["name"]
                @info "Initialize TrajAnalyzer..."
                TrajAnalyzer.init_probe_parallel!(probe_sfn)
                @info "Calculating Spectrum for 2D probe "*gm_name*"..."
                TrajAnalyzer.spectrum2d(output_prefix, gm_name)
            end
        elseif sim_type == "3D"
            for gm_name in config["spectrum"]["name"]
                @info "Initialize TrajAnalyzer..."
                TrajAnalyzer.init_probe_parallel!(probe_sfn)
                @info "Calculating Spectrum for 3D probe "*gm_name*"..."
                TrajAnalyzer.spectrum3d(output_prefix, gm_name)
            end            
        else
            @error "unrecognized simulation type: "* sim_type            
        end
    end
    if flags["movie_data_flag"]
        @info "Calculating movie potentials..."
        TrajAnalyzer.output_movie_data_2d(config["movie-output"],output_prefix*"_moviedata.h5")
    end
    #cleanup sharedarray
    TrajAnalyzer.cleanup()
    @everywhere GC.gc()
    # if flags["ngamma1d_flag"]
    #     info("Calculating N Gamma1D.")
    #     TrajAnalyzer.ngamma1d(idx,output_prefix*"_ngamma1d.mat")
    # end
end

function single_scan_scaling(trajsolver_config::Dict,config::Dict,sfn::ScalarFieldNode,probe_sfn::ScalarFieldNode,input_file,output_file,flags)
    @assert config["type"] == "single-scan-scaling"
    range_i_start = config["range_i_start"]::Int
    range_i_end = config["range_i_end"]::Int
    range_i_step = config["range_i_step"]::Int
    jobs = config["jobs"]
    for i = range_i_start:range_i_step:range_i_end
        @info "i = "*string(i)
        for job in values(jobs)
            field_name = job["field"]
            s = job["scaling"]
            s = replace(s,"@i"=>float(i))
            @info "change scaling of field $field_name to ",s
            s_exp = Meta.parse(s)
            Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        end
        Fields.init_parallel_potential!(sfn)
        input_prefix = input_file*string(i)
        output_prefix = output_file*string(i)
        job_inner_loop(config,sfn,probe_sfn,input_prefix,output_prefix,flags,[i])
        GC.gc()
    end
end

function double_scan_scaling(trajsolver_config::Dict,config::Dict,sfn::ScalarFieldNode,input_file,output_file,flags)
    @assert config["type"] == "double-scan-scaling"
    range_i_start = config["range_i_start"]::Int
    range_j_start = config["range_j_start"]::Int
    range_i_end = config["range_i_end"]::Int
    range_j_end = config["range_j_end"]::Int
    jobs = config["jobs"]
    for i = range_i_start:range_i_end, j = range_j_start:range_j_end
        @info "i = "*string(i)*", j = "*string(j)
        for job in values(jobs)
            field_name = job["field"]
            s = job["scaling"]
            s = replace(s,"@i",float(i))
            s = replace(s,"@j",float(j))
            #info("change scaling of field $field_name to ",s)
            @info "change scaling of field $field_name to $s"
            s_exp = parse(s)
            Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        end
#        Fields.init_parallel!(sfn)
        input_prefix = input_file*string(i)*"_"*string(j)
        output_prefix = output_file*string(i)*"_"*string(j)
        job_inner_loop(config,sfn,input_prefix,output_prefix,flags,[i,j])
    end
end

function optimize_ngamma1d_func(x,knobs,config::Dict,sfn::ScalarFieldNode)
    @debug "in optimize_ngamma1d_func()"
    @debug "update scaling factor..."
    for (i,k) in enumerate(knobs)
        s = "t->$(x[i])"
        @info "change scaling of field $k to $s"
        s_exp = parse(s)
        Fields.setscaling!(Fields.find_field(x->x.name==ascii(k),sfn),s_exp)
    end
    Fields.init_parallel!(sfn)
    @debug "calculate trajectories...."
    result = calculate_traj()
    probe_sfn = Fields.buildAndAlign(config["probe"]["field"],0,name=ascii([k for k in keys(config["probe"])][1]))
    TrajAnalyzer.init_parallel!(result,probe_sfn,sfn,config)
    ngamma1d = TrajAnalyzer.calc_avg_ngamma1d_parallel()*-1.0
    @info string(ngamma1d)
    return ngamma1d
end

function optimize_ngamma1d(config::Dict,sfn::ScalarFieldNode)
#    result = optimize(f,[0.5,0.5],store_trace = true,extended_trace = true)
    knobs = config["knobs"]
#    result = optimize(x->optimize_ngamma1d_func(x,knobs,config,sfn),zeros(length(knobs)), NelderMead(),OptimizationOptions(store_trace = true,extended_trace = true))
    result = opt.basinhopping(x->optimize_ngamma1d_func(x,knobs,config,sfn),zeros(length(knobs)),niter=5,stepsize=1,minimizer_kwargs=Dict("method"=>"Powell"));
    println(result)
    #result = optimize(x->optimize_ngamma1d_func(x,knobs,config,sfn),[1.0,-1.0,5.0], NelderMead(),OptimizationOptions(store_trace = true,extended_trace = true))
end

function get_large_init_range(init_range)
    #find the largest area that contain the whole init_range
    init_range = convert(Array{Float64},cat(2,init_range...))
    init_range_large = zeros(Float64,4)
    init_range_large[1] = minimum(init_range[1,:])
    init_range_large[2] = maximum(init_range[2,:])
    init_range_large[3] = minimum(init_range[3,:])
    init_range_large[4] = maximum(init_range[4,:])
    return init_range_large
end
