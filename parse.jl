using ArgParse
using YAML
using Lumberjack
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
        "--spectrum", "-S"
        help = "calculate spectrum."
        action = :store_true
        "--movie", "-M"
        help = "render movies."
        action = :store_true
        "--moviedata"
        help = "output movies potential data in a 3D array."
        action = :store_true
        "--benchmark"
        help = "benchmark Fields module."
        action = :store_true
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
    if parsed_args["trajectory"] == false && parsed_args["infile"] == ""
        Base.error("please provide input file prefix.")
    end
    println("Parsed args:")
    for (key,val) in parsed_args
        println("  $key  =>  $(repr(val))")
    end
    return parsed_args
end

function parse_config(filename,parsed_args)
    config = YAML.load(open(filename))
    simu_type = config["simulation-type"]
    fields_config = config["fields-config"]
    @assert length(keys(fields_config)) == 1 "more than 1 top level fieldnode!"
    trajsolver_config = config["trajsolver-config"]
    job_config = config["job-config"]
    Lumberjack.debug("fields config:",convert(Dict{Any,Any},copy(fields_config)))
    Lumberjack.debug("trajsolver config:",convert(Dict{Any,Any},copy(trajsolver_config)))
    Lumberjack.debug("job config:",convert(Dict{Any,Any},copy(job_config)))
    if length(parsed_args["irange"]) != 0
        job_config["range_i_start"] = parsed_args["irange"][1]
        job_config["range_i_end"] = parsed_args["irange"][2]
    end
    if length(parsed_args["jrange"]) != 0
        job_config["range_j_start"] = parsed_args["jrange"][1]
        job_config["range_j_end"] = parsed_args["jrange"][2]
    end
    ##TODO: check config format
    return fields_config, trajsolver_config, job_config
end
