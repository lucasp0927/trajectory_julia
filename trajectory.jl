using MicroLogging
configure_logging(min_level=:debug)
include("parse.jl")
parsed_args, flags = parse_commandline()

@info "Starting $(parsed_args["procs"]) processes."
#addprocs(parsed_args["procs"], exeflags=`--depwarn=no --compilecache=no`)
#addprocs(parsed_args["procs"], exeflags=`--depwarn=no`)
addprocs(parsed_args["procs"])

push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
push!(LOAD_PATH, "./TrajAnalyzer")
using Fields
using TrajSolver
using TrajAnalyzer
include("fileio.jl")
include("job_manage.jl")

function prepare()
    config_file,input_file,output_file = parsed_args["config"],parsed_args["infile"],parsed_args["outfile"]
    fields_config,trajsolver_config,job_config = parse_config(config_file,parsed_args)
    #info("building field $([k for k in keys(fields_config)][1])...")
    name = ascii([k for k in keys(fields_config)][1])
    @info "building field $name"
    sfn = Fields.buildAndAlign(fields_config["field"],0,name=ascii([k for k in keys(fields_config)][1]))
    probe_sfn = Fields.buildAndAlign(job_config["probe"]["field"],0,name=ascii([k for k in keys(job_config["probe"])][1]))
    return sfn,probe_sfn,input_file,output_file,job_config,trajsolver_config    
end

function main()
    sfn,probe_sfn,input_file,output_file,job_config,trajsolver_config = prepare()
    TrajSolver.init_parallel(trajsolver_config,probe_sfn)    
    @info "initialize fields"
    Fields.init_parallel!(sfn)
    @info "Start calculating trajectories..."
    if job_config["type"] == "single-scan-scaling"
        @info "job type: single-scan-scaling"
        single_scan_scaling(trajsolver_config,job_config,sfn,probe_sfn,input_file,output_file,flags)
    elseif job_config["type"] == "double-scan-scaling"
        @info "job type: double-scan-scaling"
        double_scan_scaling(trajsolver_config,job_config,sfn,input_file,output_file,flags)
    # elseif job_config["type"] == "optimization"
    #     @info "job type: optimization"
    #     optimize_ngamma1d(job_config,sfn)
    end
end
main()
