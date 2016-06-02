push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
using Fields
using TrajSolver
include("fileio.jl")
include("parse.jl")
include("job_manage.jl")

function prepare()
    parsed_args = parse_commandline()
    config_file,output_file = parsed_args["config"],parsed_args["outfile"]
    fields_config,trajsolver_config,job_config = parse_config(config_file,false)
    TrajSolver.init_parallel(trajsolver_config)
    #build field with name "field"
    @assert length(keys(fields_config)) == 1 "more than 1 top level fieldnode!"
    println("building field ",[k for k in keys(fields_config)][1],"...")
    sfn = Fields.build_field(fields_config["field"],0,true,name=ascii([k for k in keys(fields_config)][1]))
    println("aligning...")
    Fields.align_field_tree!(sfn)
    Fields.set_geometry!(sfn)
    Fields.set_typeof!(sfn)
    return sfn,output_file,job_config
end

function main()
    #preparation
    sfn,output_file,job_config = prepare()
    println("Start calculating trajectories...")
    println("initialize fields")
    Fields.init_parallel!(sfn)
    single_scan_scaling(job_config,sfn,output_file)
end
main()
