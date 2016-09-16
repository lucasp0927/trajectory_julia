include("flux.jl")
function job_inner_loop(config,sfn,input_prefix,output_prefix,flags)
    if flags["calc_traj_flag"]
        result = calculate_traj(i)
        Lumberjack.info("save results...")
        matwrite(output_prefix*".mat",result)
        traj = result["traj"]
        tspan = result["tspan"]
    elseif flags["spectrum_flag"] || flags["movie_flag"]
        Lumberjack.info("read results...")
        result = matread(input_prefix*".mat")
        traj = result["traj"]
        tspan = result["tspan"]
    end
    if flags["spectrum_flag"] || flags["movie_flag"]
        #these function use trajectories data, so TrajAnalyzer needs to be initialized.
        Lumberjack.info("Initialize TrajAnalyzer...")
        probe_sfn = Fields.buildAndAlign(config["probe"]["field"],0,name=ascii([k for k in keys(config["probe"])][1]))
        TrajAnalyzer.init_parallel!(result,probe_sfn,sfn,config)
    end
    if flags["movie_flag"]
        Lumberjack.info("Outputing Movie...")
        movie_range = [promote(config["movie-output"]["range"]...)...]
        TrajAnalyzer.output_movie_traj(config["movie-output"],output_prefix*"_traj.mp4")
    end
    if flags["spectrum_flag"]
        Lumberjack.info("Calculating Spectrum...")
        TrajAnalyzer.spectrum(output_prefix)
    end
    if flags["movie_data_flag"]
        Lumberjack.info("Calculating movie potentials...")
        TrajAnalyzer.output_movie_data(config["movie-output"],output_prefix*"_moviedata.mat")
    end
end

function single_scan_scaling(trajsolver_config::Dict,config::Dict,sfn::ScalarFieldNode,input_file,output_file,flags)
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
        job_inner_loop(config,sfn,input_file*string(i),output_file*string(i),flags)
    end
end

function double_scan_scaling(trajsolver_config::Dict,config::Dict,sfn::ScalarFieldNode,input_file,output_file,flags)
    range_i_start = config["range_i_start"]
    range_j_start = config["range_j_start"]
    range_i_end = config["range_i_end"]
    range_j_end = config["range_j_end"]
    jobs = config["jobs"]
    for i = range_i_start:range_i_end, j = range_j_start:range_j_end
        for job in values(jobs)
            field_name = job["field"]
            scaling = job["scaling"]
            s = replace(scaling,"@i",float(i))
            s = replace(s,"@j",float(j))
            Lumberjack.info("change scaling of field $field_name to ",s)
            s_exp = eval(parse(s))
            Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        end
        Fields.init_parallel!(sfn)
        input_prefix = input_file*string(i)*"_"*string(j)
        output_prefix = output_file*string(i)*"_"*string(j)
        job_inner_loop(config,sfn,input_prefix,output_prefix,flags)
    end
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
