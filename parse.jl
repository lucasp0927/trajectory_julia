using ArgParse
using YAML

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--config", "-C"
        help = "Configuration file."
        required = true
        "--outfile", "-O"
        help = "output file."
        required = true
    end
    return parse_args(s)
end

function parse_config(filename,verbose)
    config = YAML.load(open(filename))
    fields_config = config["fields-config"]
    trajsolver_config = config["trajsolver-config"]
    job_config = config["job-config"]
    if verbose == true
        display(fields_config);println("")
        display(trajsolver_config);println("")
    end
    ##TODO: check config format
    return fields_config, trajsolver_config, job_config
end
