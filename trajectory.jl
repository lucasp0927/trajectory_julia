include("parse.jl")
parsed_args = parse_commandline()
println("Starting ",parsed_args["procs"]," processes.")
addprocs(parsed_args["procs"], exeflags="--depwarn=no")

push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
push!(LOAD_PATH, "./TrajAnalyzer")
using Fields
using TrajSolver
using TrajAnalyzer
@everywhere using Lumberjack
using Lumberjack
include("fileio.jl")
include("job_manage.jl")

function prepare()
#    parsed_args = parse_commandline()
    config_file,input_file,output_file = parsed_args["config"],parsed_args["infile"],parsed_args["outfile"]
    fields_config,trajsolver_config,job_config = parse_config(config_file)
    TrajSolver.init_parallel(trajsolver_config)
    @assert length(keys(fields_config)) == 1 "more than 1 top level fieldnode!"
    println("building field ",[k for k in keys(fields_config)][1],"...")
    sfn = Fields.buildAndAlign(fields_config["field"],0,name=[k for k in keys(fields_config)][1])
    return sfn,input_file,output_file,job_config,parsed_args["trajectory"],parsed_args["spectrum"],parsed_args["movie"],trajsolver_config
end

function main()
    @everywhere Lumberjack.remove_truck("console")
    #Lumberjack.remove_truck("console")
    @everywhere Lumberjack.add_truck(LumberjackTruck(STDOUT, "info"))
#    Lumberjack.add_truck(LumberjackTruck(STDOUT, "debug"))
    @everywhere Lumberjack.add_truck(LumberjackTruck("trajectory_logfile.log","debug"))
    #preparation
    sfn,input_file,output_file,job_config,calc_traj_flag,spectrum_flag,movie_flag,trajsolver_config = prepare()
    println("initialize fields")
    Fields.init_parallel!(sfn)
    println("Start calculating trajectories...")
    single_scan_scaling(trajsolver_config,job_config,sfn,input_file,output_file,calc_traj_flag,spectrum_flag,movie_flag)
end
main()
