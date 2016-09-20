using Logging
Logging.configure(level=Logging.INFO)
include("parse.jl")
parsed_args = parse_commandline()
info("Starting ",parsed_args["procs"]," processes.")
addprocs(parsed_args["procs"], exeflags="--depwarn=no")

push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
push!(LOAD_PATH, "./TrajAnalyzer")
using Fields
using TrajSolver
using TrajAnalyzer

include("fileio.jl")
include("job_manage.jl")

function prepare()
#    parsed_args = parse_commandline()
    config_file,input_file,output_file = parsed_args["config"],parsed_args["infile"],parsed_args["outfile"]
    flags = Dict("calc_traj_flag" => parsed_args["trajectory"],
                 "spectrum_flag" => parsed_args["spectrum"],
                 "movie_flag" => parsed_args["movie"],
                 "movie_data_flag" => parsed_args["moviedata"],
                 "benchmark_flag" => parsed_args["benchmark"])
    fields_config,trajsolver_config,job_config = parse_config(config_file,parsed_args)
    println("building field ",[k for k in keys(fields_config)][1],"...")
    sfn = Fields.buildAndAlign(fields_config["field"],0,name=ascii([k for k in keys(fields_config)][1]))
    return sfn,input_file,output_file,job_config,trajsolver_config,flags
end

function main()
    #preparation
    sfn,input_file,output_file,job_config,trajsolver_config,flags = prepare()
    TrajSolver.init_parallel(trajsolver_config)
    info("initialize fields")
    Fields.init_parallel!(sfn)
    info("Start calculating trajectories...")
    if job_config["type"] == "single-scan-scaling"
        info("single-scan-scaling")
        single_scan_scaling(trajsolver_config,job_config,sfn,input_file,output_file,flags)
    elseif job_config["type"] == "double-scan-scaling"
        info("double-scan-scaling")
        double_scan_scaling(trajsolver_config,job_config,sfn,input_file,output_file,flags)
    end
end
main()
