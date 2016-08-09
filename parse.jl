using ArgParse
using YAML
using Lumberjack
function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config", "-C"
        help = "Configuration file."
        required = true
        "--outfile", "-O"
        help = "output file."
        required = true
        "--trajectory", "-T"
        help = "calculate trajectories"
        action = :store_true
        "--spectrum", "-S"
        help = "calculate spectrum"
        action = :store_true
        "--movie", "-M"
        help = "render movies"
        action = :store_true
    end
    return parse_args(s)
end

function parse_config(filename)
    config = YAML.load(open(filename))
    simu_type = config["simulation-type"]
    fields_config = config["fields-config"]
    trajsolver_config = config["trajsolver-config"]
    job_config = config["job-config"]
    Lumberjack.debug("fields config:",convert(Dict{Any,Any},copy(fields_config)))
    Lumberjack.debug("trajsolver config:",convert(Dict{Any,Any},copy(trajsolver_config)))
    Lumberjack.debug("job config:",convert(Dict{Any,Any},copy(job_config)))
    ##TODO: check config format
    return fields_config, trajsolver_config, job_config
end
