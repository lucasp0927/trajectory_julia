include("parse.jl")
parsed_args, flags = parse_commandline()

@info "Starting $(parsed_args["procs"]) processes."
using Distributed
addprocs(parsed_args["procs"])

@everywhere push!(LOAD_PATH, "./Fields")
@everywhere push!(LOAD_PATH, "./TrajSolver")
@everywhere push!(LOAD_PATH, "./TrajAnalyzer")
using Fields
using TrajSolver
using TrajAnalyzer
include("fileio.jl")
include("job_manage.jl")

function prepare()
    config_file,input_file,output_file = parsed_args["config"],parsed_args["infile"],parsed_args["outfile"]
    fields_config,trajsolver_config,job_config,material_config = parse_config(config_file,parsed_args)
    #info("building field $([k for k in keys(fields_config)][1])...")
    name = ascii([k for k in keys(fields_config)][1])
    @info "building field $name"
    sfn = Fields.buildAndAlign(fields_config["field"],0,name=ascii([k for k in keys(fields_config)][1]))
#    Fields.save_field(sfn)
    probe_sfn = Fields.buildAndAlign(job_config["probe"]["field"],0,name=ascii([k for k in keys(job_config["probe"])][1]))
    mat_sfn = Fields.buildAndAlign(material_config["material"],0,name=ascii([k for k in keys(material_config)][1]))
    return sfn,probe_sfn,mat_sfn,input_file,output_file,job_config,trajsolver_config
end

function main()
    sfn,probe_sfn,mat_sfn,input_file,output_file,job_config,trajsolver_config = prepare()
    @info "initializing TrajSolver..."
    TrajSolver.init_parallel(trajsolver_config,probe_sfn)
    Fields.init_parallel!(sfn, mat_sfn)
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
