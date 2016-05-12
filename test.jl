push!(LOAD_PATH, "./Fields")
#using FastAnonymous
using Fields
using MAT
#using Plots
#plotly()
function readfields(filename::AbstractString,variable::AbstractString)
    matfile = matopen(filename)
    if exists(matfile, variable)
#        println("Reading $variable from $filename...")
        var = read(matfile, variable)
    else
        error("Can't read $variable from $filename")
    end
    close(matfile)
    return var
end

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

function buildfields!(rb_field_s::SharedArray,lb_field_s::SharedArray,gm_field_s::SharedArray)
    rb = Fields.VectorField{Complex{Float64},2}(rb_field_s,[0.0,0.0],[6666*15.0,3333*15.0],scaling=t-> exp(1.0im*t*2pi))
    lb = Fields.VectorField{Complex{Float64},2}(lb_field_s,[0.0,0.0],[6666*15.0,3333*15.0],scaling= t-> 1.0+0.0im)
    vfn = Fields.VectorFieldNode{2}([lb,rb],scaling= t->1.0+0.0im)
    gm = Fields.ScalarField{Float64,2}(gm_field_s,[(6666*15.0)/2.0-1854.625,(3333*15.0)/2.0-1850.0],[3690.75,3690.75],scaling= t-> 10.0)
    sfn = Fields.ScalarFieldNode{2}([vfn,gm])
    info("aligning...")
    Fields.align_field_tree!(sfn)
    Fields.set_geometry!(sfn)
    Fields.set_typeof!(sfn)
    return sfn
end


function test(sfn::ScalarFieldNode)
    Fields.init_parallel!(sfn)
    sfn = 0
    gc()
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
    @time  benchmark_smp(Fields.fields)
    println("benchmark value2 sfn")
    #=
    @time    benchmark_value(sfn)
    println("diff: ",mean(output1-output2))
    =#
    println("output")
    @time  r = @spawnat 2 itp_test()
    output = fetch(r)
    file = matopen("out.mat", "w")
    write(file, "itp", output)
    close(file)
#    heatmap(output)
 #   png("output")
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

function benchmark_smp(f)
    for i = 1:1000000
        Fields.sample2!(f,[50000.0+rand()*10.0,25000.0+rand()*10.0],1.0*rand())
    end
end

function benchmark_value(f)
    for i = 1:1000000
        Fields.value2(f,[50000.0+rand()*10.0,25000.0+rand()*10.0],1.0*rand())
    end
end

@everywhere function itp_test()
    info("test")
    N = 1000
    xx = linspace(48000-1000,52000+1000,N)
    yy = linspace(23000,27000,N)
#    xx = linspace(450-400,450+400,N)
#    yy = linspace(23000,27000,N)
    output = zeros(Float64,(N,N))
    for x in enumerate(xx)
        for y in enumerate(yy)
            output[x[1],y[1]] = Fields.value3([x[2],y[2]],0.75)
        end
    end
    return output;
end

function itp_grad_test(sfn)
    N = 1000
    xx = linspace(48000,52000,N)
    yy = linspace(23000,27000,N)
    output = zeros(Float64,(2,N,N))
    for x in enumerate(xx)
        for y in enumerate(yy)
            output[:,x[1],y[1]] = Fields.gradient(sfn,[x[2],y[2]],0.5)
        end
    end
    file = matopen("out.mat", "w")
    write(file, "itp", output)
    close(file)
    return output
end

function main()
    println(nprocs()," processes running.")
    println("######################################")
#    println("building sf1")
#    sf1 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->sin(x/900.0*2pi),[1000,1000],[0.0,0.0],[900.0,900.0])
#    println("building sf2")
#    sf2 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->-sin(x/900.0*2pi),[101,101],[299.0,299.0],[300.0,300.0])
    println("building lattice beams")
    println("right beam")
    file = matopen("lattice_right.mat")
    rb_field = read(file, "beam_right") # note that this does NOT introduce a variable ``varname`` into scope
    rb_field_s = Fields.copy_to_sharedarray!(rb_field)
    rb_field = []
    close(file)

    println("left beam")
    file = matopen("lattice_left.mat")
    lb_field = read(file, "beam_left") # note that this does NOT introduce a variable ``varname`` into scope
    lb_field_s = Fields.copy_to_sharedarray!(lb_field)
    lb_field = []
    close(file)

    println("gm")
    file = matopen("D2_TE.mat")
    gm_field = read(file, "gm") # note that this does NOT introduce a variable ``varname`` into scope
    gm_field_s = Fields.copy_to_sharedarray!(gm_field)
    gm_field = []
    close(file)
    gc()
    #########
    sfn = buildfields!(rb_field_s,lb_field_s,gm_field_s)
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
whos()
#@time main()
