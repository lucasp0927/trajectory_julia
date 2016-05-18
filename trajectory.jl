push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
using Fields
using TrajSolver
include("fileio.jl")
include("parse.jl")

function main()
    println(nprocs()," processes running.")
    parsed_args = parse_commandline()
    config_file = parsed_args["config"]
    println("config file: ", config_file)
    output_file = parsed_args["outfile"]
    println("output file: ", output_file)
    fields_config,trajsolver_config = parse_config(config_file,true)
    TrajSolver.init_parallel(trajsolver_config)
    sfn = Fields.build_field(fields_config["field"],"field",0,true)
    println("aligning...")
    Fields.align_field_tree!(sfn)
    Fields.set_geometry!(sfn)
    Fields.set_typeof!(sfn)
    println("Start calculating trajectories...")

    for scale = (1.0:1.0:50.0).^2
        println("scale: ",scale)
        Fields.setscaling!(sfn,t->scale)
        Fields.init_parallel!(sfn)
        @time    @sync begin
            for p = 2:nprocs()
                @async remotecall_wait(p,TrajSolver.solve_traj)
            end
        end
        temp = cell(nworkers())
        for p = 2:nprocs()
            temp[p-1] = remotecall_fetch(p,TrajSolver.get_result)
        end
        result = cat(3,temp...)
        savemat(output_file*string(sqrt(scale))*".mat",result,"result")
    end
end
main()
