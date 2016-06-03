push!(LOAD_PATH, "./Fields")
using Fields
using MAT
include("flux.jl")
include("parse.jl")
include("fileio.jl")
function main()
    println(ARGS)
    result = matread(ARGS[2])
    fields_config,trajsolver_config,job_config = parse_config(ARGS[1],false)
    @time calc_flux(result["traj"],result["tspan"],job_config["flux"],ARGS[3])
    @time calc_flux(result["traj"],result["tspan"],job_config["flux"],ARGS[3])
end

main()
