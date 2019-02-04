using ArgParse
using YAML
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config", "-C"
        help = "Configuration file."
        required = true
        "--infile", "-I"
        help = "input file prefix"
        default = ""
        "--outfile", "-O"
        help = "output file prefix."
        required = true
        "--trajectory", "-T"
        help = "calculate trajectories."
        action = :store_true
        "--trajinit"
        help = "just initialize TrajAnalyzer"
        action = :store_true
        "--spectrum", "-S"
        help = "calculate spectrum."
        action = :store_true
        "--movie", "-M"
        help = "render movies."
        action = :store_true
        "--moviedata"
        help = "output movies potential data in a 3D array."
        action = :store_true
        # "--benchmark"
        # help = "benchmark Fields module."
        #action = :store_true
        "--ngamma1d"
        help = "Calculate N*Gamma1D, and output a range_i by range_j matrix."
        action = :store_true
        # "--optim"
        # help = "Optimize"
        # action = :store_true
        "--procs", "-P"
        arg_type = Int
        required = true
        help = "processes number."
        "--irange"
        help = "override i range in config file."
        nargs = 2
        arg_type = Int
        "--jrange"
        help = "override j range in config file."
        nargs = 2
        arg_type = Int
    end
    parsed_args = parse_args(s)
    need_traj_flag = parsed_args["spectrum"] || parsed_args["movie"] || parsed_args["ngamma1d"] || parsed_args["trajinit"]
    if need_traj_flag && (parsed_args["trajectory"] == false && parsed_args["infile"] == "")
        @error "please provide input file prefix."
    end
    @debug "Parsed args:"
    for (key,val) in parsed_args
        @debug "  $key  =>  $(repr(val))"
    end
    flags = Dict("calc_traj_flag" => parsed_args["trajectory"],
                 "spectrum_flag" => parsed_args["spectrum"],
                 "movie_flag" => parsed_args["movie"],
                 "movie_data_flag" => parsed_args["moviedata"],
#                 "benchmark_flag" => parsed_args["benchmark"],
                 "ngamma1d_flag" => parsed_args["ngamma1d"],
                 "trajinit_flag" => parsed_args["trajinit"],
                 "need_traj_flag" => need_traj_flag)
    return parsed_args, flags
end

function parse_config(filename,parsed_args)
    config = YAML.load(open(filename))
    sim_type = config["simulation-type"]
    fields_config = config["fields-config"]
    if sim_type == "3D"
        material_config = config["material-config"]
    elseif sim_type == "2D"
        material_config =Dict{Any,Any}("material"=>
            Dict{Any,Any}(
            "scaling"    => "t->0.0",
            "field-type" => "ScalarFieldNode",
            "fields"     => Dict{Any,Any}("sphere"=>Dict{Any,Any}("scaling"=>"t->0.0","field-type"=>"ScalarField","D-type"=>"Float","init-type"=>"zero","dim"=>2,"pos"=>[0,0],"size"=>[10,10],"res"=>[1,1])),
            "dim"        => 2
        ))
    end
    @assert length(keys(fields_config)) == 1 "more than 1 top level fieldnode!"
    trajsolver_config = config["trajsolver-config"]
    trajsolver_config["simulation-type"] = sim_type
    job_config = config["job-config"]
    job_config["movie-output"]["simulation-type"] = sim_type
    @debug "fields config:" convert(Dict{Any,Any},copy(fields_config))
    @debug "trajsolver config:" convert(Dict{Any,Any},copy(trajsolver_config))
    @debug "job config:" convert(Dict{Any,Any},copy(job_config))
    if length(parsed_args["irange"]) != 0
        job_config["range_i_start"] = parsed_args["irange"][1]
        job_config["range_i_end"] = parsed_args["irange"][2]
    end
    if length(parsed_args["jrange"]) != 0
        job_config["range_j_start"] = parsed_args["jrange"][1]
        job_config["range_j_end"] = parsed_args["jrange"][2]
    end
    ##TODO: check config format
    return fields_config, trajsolver_config, job_config, material_config
end
