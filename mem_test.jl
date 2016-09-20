using Logging
Logging.configure(level=Logging.INFO)
#include("trajectory.jl")
include("parse.jl")
push!(LOAD_PATH, "./Fields")
using Fields
function parse_config(filename)
    config = YAML.load(open(filename))
    simu_type = config["simulation-type"]
    fields_config = config["fields-config"]
    @assert length(keys(fields_config)) == 1 "more than 1 top level fieldnode!"
    trajsolver_config = config["trajsolver-config"]
    job_config = config["job-config"]
    debug("fields config:",convert(Dict{Any,Any},copy(fields_config)))
    debug("trajsolver config:",convert(Dict{Any,Any},copy(trajsolver_config)))
    debug("job config:",convert(Dict{Any,Any},copy(job_config)))
    ##TODO: check config format
    return fields_config, trajsolver_config, job_config
end
parse_config("/home/lucaspeng/Desktop/config_20160919.yml")
fields_config,trajsolver_config,job_config = parse_config("/home/lucaspeng/Desktop/config_20160919.yml")
sfn = Fields.buildAndAlign(fields_config["field"],0,name=ascii([k for k in keys(fields_config)][1]))
include("benchmark.jl")
