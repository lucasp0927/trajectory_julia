push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
using Fields
using TrajSolver
include("fileio.jl")
include("benchmark.jl")
include("parse.jl")

function test(sfn::ScalarFieldNode)
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
    savemat("result.mat",result,"result")
#    output = Fields.composite_slow([200.0, 600.0, 100.0, 49895.0],0.0)
#    savemat("out.mat",output,"output")
    # @time    begin
    #     println("value")
    #     r =@spawnat 2 itp_test()
    #     println("gradient")
    #     r_grad = @spawnat 3 itp_grad_test2()
    #     output = fetch(r)

    #     grad = fetch(r_grad)
    #     savemat("out.mat",output,"output")
    #     savemat("grad_out.mat",grad,"grad")
    # end
end

function main()
    println(nprocs()," processes running.")
    fields_config,trajsolver_config = parse_config("config_200khz_lattice.yml",true)
    TrajSolver.init_parallel(trajsolver_config)
    field_config = fields_config["field"]
    sfn = Fields.build_field(field_config,"field",0,true)
    println("aligning...")
    Fields.align_field_tree!(sfn)
    Fields.set_geometry!(sfn)
    Fields.set_typeof!(sfn)
#    Profile.init(delay=0.01)
    test(sfn)
#    Profile.clear()
#    open("profile.bin", "w") do f serialize(f, Profile.retrieve()) end
end
main()
#=
function sendto(p::Int; args...)
    for (nm, val) in args
        @spawnat(p, eval(Main, Expr(:(=), nm, val)))
    end
end

function sendto(ps::Vector{Int}; args...)
    for p in ps
        sendto(p; args...)
    end
end
=#
