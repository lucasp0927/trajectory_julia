push!(LOAD_PATH, "./Fields")
push!(LOAD_PATH, "./TrajSolver")
using Fields
using TrajSolver
include("fileio.jl")
include("benchmark.jl")
include("parse.jl")

function buildfields!(rb_field_s::SharedArray,lb_field_s::SharedArray,gm_field_s::SharedArray)
    rb = Fields.VectorField{Complex{Float64},2}(rb_field_s,[0.0,0.0],[6666*15.0,3333*15.0],scaling= t->exp(2*pi*im*t))
    lb = Fields.VectorField{Complex{Float64},2}(lb_field_s,[0.0,0.0],[6666*15.0,3333*15.0],scaling= t-> 1.0+0.0im)
    vfn = Fields.VectorFieldNode{2}([lb,rb],scaling= t->1.0+0.0im)
    gm = Fields.ScalarField{Float64,2}(gm_field_s,[(6666*15.0)/2.0-1854.625,(3333*15.0)/2.0-1850.0],[3690.75,3690.75],scaling=  t-> 10.0)
#    gm = Fields.func2field(ScalarField{Float64,2},(x,y)->sin(x+y),[400,400],[(6666*15.0)/2.0-1854.625,(3333*15.0)/2.0-1850.0],[3690.75,3690.75],scaling=  t-> 10.0)
    sfn = Fields.ScalarFieldNode{2}([vfn,gm])
    info("aligning...")
    Fields.align_field_tree!(sfn)
    Fields.set_geometry!(sfn)
    Fields.set_typeof!(sfn)
    return sfn
end

function test(sfn::ScalarFieldNode)
    Fields.init_parallel!(sfn)
    ######
    # Fields.align_field_tree!(Fields.fields)
    # Fields.set_geometry!(Fields.fields)
    # Fields.set_typeof!(Fields.fields)

    # println("geometry of lb")
    # println(Fields.geometry(lb))
    # println("geometry of rb")
    # println(Fields.geometry(rb))
    # println("geometry of gm")
    # println(Fields.geometry(gm))
    # println("geometry of sfn")
    # println(Fields.geometry(sfn))
    ######gradient
    #=
    @time    Fields.sample(gm,[49140.0,24147.0],1.0)
    @time    output1 = Fields.sample(vfn,[1000.0,1000.0],1.0)
    @time    output2 = Fields.sample(vfn1,[1000.0,1000.0],1.0)

    println("benchmark sampling vfn")
    @time    benchmark_smp(vfn)
    println("benchmark sampling sfn")
=#
    # @time  benchmark_smp(Fields.fields)
    # println("benchmark value2 sfn")

    # @time  benchmark_value()
    #=
    println("diff: ",mean(output1-output2))
    =#
@time    begin
        println("value")
        r =@spawnat 2 itp_test()
        println("gradient")
        r_grad = @spawnat 3 itp_grad_test2()
        output = fetch(r)
        grad = fetch(r_grad)
        savemat("out.mat",output,"output")
        savemat("grad_out.mat",grad,"grad")
    end
    ######composite
# @time    f_out = Fields.composite(sfn,0.0)
#     file = matopen("comp_0.0.mat", "w")
#     write(file, "field", f_out.field)
#     close(file)
# @time    f_out = Fields.composite(sfn,0.25)
#     file = matopen("comp_0.25.mat", "w")
#     write(file, "field", f_out.field)
#     close(file)
# @time    f_out = Fields.composite(sfn,0.5)
#     file = matopen("comp_0.5.mat", "w")
#     write(file, "field", f_out.field)
#     close(file)
end

function main()
    println(nprocs()," processes running.")
    println("######################################")
#    println("building sf1")
#    sf1 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->sin(x/900.0*2pi),[1000,1000],[0.0,0.0],[900.0,900.0])
#    println("building sf2")
    #    sf2 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->-sin(x/900.0*2pi),[101,101],[299.0,299.0],[300.0,300.0])
    #=
    println("building lattice beams")
    println("right beam")
    rb_field_s = mat2sharedarray("lattice_right.mat","beam_right")
    println("left beam")
    lb_field_s = mat2sharedarray("lattice_left.mat","beam_left")
    println("gm")
    gm_field_s = mat2sharedarray("D2_TE.mat","gm")
    sfn = buildfields!(rb_field_s,lb_field_s,gm_field_s)
    =#
    fields_config,trajsolver_config = parse_config("config.yml",false)
    field_config = fields_config["field"]
    sfn = Fields.build_field(field_config,"field",0,true)
    println("aligning...")
    Fields.align_field_tree!(sfn)
    Fields.set_geometry!(sfn)
    Fields.set_typeof!(sfn)
#    Profile.init(delay=0.01)
    test(sfn)
#    Profile.clear()
    gc()
    test(sfn)
#    open("profile.bin", "w") do f serialize(f, Profile.retrieve()) end
    ##########3D test

#     println("######################################")
#     println("3D align test")
#     sf1 = Fields.zero(ScalarField{Float64,3},(1001,1001,11),(0,0,0),(1000,1000,10))
#     sf2 = Fields.zero(ScalarField{Float64,3},(101,101,11),(299,299,0),(300,300,10))
#     vf1 = Fields.func2field(VectorField{Complex{Float64},3},(x,y,z)->[sin(x)+cos(x)im,1.0im,1.0im],(101,101,11),(299,299,0),(300,300,10),scaling=t->t)
#     #########
#     println("aligning")
#     sfn = Fields.ScalarFieldNode{3}([sf1,sf2,vf1])
#     println("geometry of sf1")
#     println(Fields.geometry(sf1))
#     println("geometry of sf2")
#     println(Fields.geometry(sf2))
#     println("geometry of vf1")
#     println(Fields.geometry(vf1))
#     println(mean(vf1.field))
#     println("geometry of sfn")
#     println(Fields.geometry(sfn))
#     println("align sfn")
# @time    Fields.align_field_tree!(sfn)
#     println("geometry of sf1")
#     println(Fields.geometry(sf1))
#     println("geometry of sf2")
#     println(Fields.geometry(sf2))
#     println("geometry of vf1")
#     println(Fields.geometry(vf1))
#     println(mean(vf1.field))
#     println("geometry of sfn")
#     println(Fields.geometry(sfn))
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
