push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
using Fields
using TrajSolver
include("fileio.jl")
include("parse.jl")
include("job_manage.jl")
function main()
    #preparation
    println(nprocs()," processes running.")
    parsed_args = parse_commandline()
    config_file = parsed_args["config"]
    println("config file: ", config_file)
    output_file = parsed_args["outfile"]
    println("output file: ", output_file)
    fields_config,trajsolver_config,job_config = parse_config(config_file,true)
    println(trajsolver_config)
    TrajSolver.init_parallel(trajsolver_config)
    #build field with name "field"
    sfn = Fields.build_field(fields_config["field"],0,true)
    println("aligning...")
    Fields.align_field_tree!(sfn)
    Fields.set_geometry!(sfn)
    Fields.set_typeof!(sfn)
    println("Start calculating trajectories...")
    Fields.init_parallel!(sfn)
#    single_scan_scaling(job_config,sfn)
    # println("find field")
    # println(Fields.find_field(x->x.name==ascii("field")),sfn)
    # println("find lattice-beam")
    # println(Fields.find_field(x->x.name==ascii("lattice-beam")),sfn)
    # println("find gm")
    # println(Fields.find_field(x->x.name==ascii("gm")),sfn)
    # println("find non")
    # println(Fields.find_field(x->x.name==ascii("non")),sfn)
    for scale = 100.0^2
        println("scale: ",scale)
        Fields.setscaling!(Fields.find_field(x->x.name==ascii("field"),sfn),t->scale)
        println(sfn.scaling(0))
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
    output = Fields.composite_slow([40000.0, 60000.0, 15000.0, 35000.0],0.0)
    savemat("out.mat",output,"output")
end
main()
