push!(LOAD_PATH, "./Fields")
using Fields
using MAT
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

function test()
    println("######################################")
    println("2D composite test")
    println("building sf1")
    sf1 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->sin(x/900.0*2pi),[1000,1000],[0.0,0.0],[900.0,900.0])
    println("building sf2")
    sf2 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->-sin(x/900.0*2pi),[101,101],[299.0,299.0],[300.0,300.0])
    println("building lattice beams")
    println("left beam")
    file = matopen("lattice_left.mat")
    lb_field = read(file, "beam_left") # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    lb = VectorField{Complex{Float64},2}(lb_field,[0.0,0.0],[6666*15.0,3333*15.0],scaling=t->1.0)
    println("right beam")
    file = matopen("lattice_right.mat")
    rb_field = read(file, "beam_right") # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    rb = VectorField{Complex{Float64},2}(rb_field,[0.0,0.0],[6666*15.0,3333*15.0],scaling=t->exp(1.0im*t*2pi))
    vfn = Fields.VectorFieldNode{2}([lb,rb])
    vfn1 = Fields.VectorFieldNode{2}([vfn])
    println("building gm")
    file = matopen("D2_TE.mat")
    gm_field = read(file, "gm") # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    gm = ScalarField{Float64,2}(gm_field,[(6666*15.0)/2.0-1854.625,(3333*15.0)/2.0-1850.0],[3690.75,3690.75],scaling=t->10.0)
    #########
    sfn = Fields.ScalarFieldNode{2}([vfn,gm])
    println("aligning...")
    @time Fields.align_field_tree!(sfn)
    @time Fields.set_geometry!(sfn)
    @time Fields.set_typeof!(sfn)
    gc()
    println("geometry of lb")
    println(Fields.geometry(lb))
    println("geometry of rb")
    println(Fields.geometry(rb))
    println("geometry of gm")
    println(Fields.geometry(gm))
    println("geometry of sfn")
    println(Fields.geometry(sfn))
    ######gradient
    #=
    @time    Fields.sample(gm,[49140.0,24147.0],1.0)
    @time    output1 = Fields.sample(vfn,[1000.0,1000.0],1.0)
    @time    output2 = Fields.sample(vfn1,[1000.0,1000.0],1.0)

    println("benchmark sampling vfn")
    @time    benchmark_smp(vfn)
    println("benchmark sampling sfn")
=#
@time    @profile    benchmark_smp(sfn)    
    println("benchmark value2 sfn")
    #=
    @time    benchmark_value(sfn)
    println("diff: ",mean(output1-output2))
    =#
    println("output")
    @time   itp_test(sfn)
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

function itp_test(sfn)
    N = 1000
    xx = linspace(48000,52000,N)
    yy = linspace(23000,27000,N)
    output = zeros(Float64,(N,N))
    for x in enumerate(xx)
        for y in enumerate(yy)
            output[x[1],y[1]] = Fields.value2(sfn::ScalarFieldNode,[x[2],y[2]],0.25)
        end
    end
    file = matopen("out.mat", "w")
    write(file, "itp", output)
    close(file)
end

function itp_grad_test(sfn)
    N = 500
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
end

function main()
#    Profile.init(delay=0.01)
    test()
    Profile.clear()
    test()
    open("profile.bin", "w") do f serialize(f, Profile.retrieve()) end
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
#@time main()
