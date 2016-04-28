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

function main()
#    Fields.initialize((100,100),(100,100))
    ########2D test
    vf1 = Fields.zero(VectorField{Float64,2},(10,10),(2,2),(9,9),scaling=x->sin(x))
    vf2 = Fields.zero(VectorField{Float64,2},(10,10),(1,1),(9,9))
    vf3 = Fields.func2field(VectorField{Float64,2},(x,y)->[x,y,0],(2,2),(0,0),(2,2))
    vf4 = Fields.func2field(VectorField{Float64,2},(x,y)->[x,y,0],(2,2),(2,2),(4,12))
    #########################################
    vfn = Fields.VectorFieldNode{2}([vf3,vf4])
    vfn2 = Fields.VectorFieldNode{2}([vf1,vf2])
    vfn3 = Fields.VectorFieldNode{2}([vfn,vfn2])
    ###############
    sf1 = Fields.zero(ScalarField{Float64,2},(10,10),(2,2),(9,9))
    sf2 = Fields.zero(ScalarField{Float64,2},(10,10),(1,1),(9,9))
    sf3 = Fields.func2field(ScalarField{Float64,2},(x,y)->x,(2,2),(0,0),(2,2))
    sf4 = Fields.func2field(ScalarField{Float64,2},(x,y)->x,(2,2),(2,2),(4,12))
    ############
    sfn = Fields.ScalarFieldNode{2}([sf3,sf4])
    sfn2 = Fields.ScalarFieldNode{2}([sf1,sf2])
    sfn3 = Fields.ScalarFieldNode{2}([sfn,sfn2])
    ##########3D test
    println("######################################")
    println("2D align test")
    sf1 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->1.0+0im,(1001,1001),(0,0),(1000,1000))
    sf2 = Fields.func2field(ScalarField{Complex{Float64},2},(x,y)->2.0+0im,(101,101),(299,299),(300,300))
    vf1 = Fields.func2field(VectorField{Complex{Float64},2},(x,y)->[2.0+0im,2.0+0im,2.0+0im],(101,101),(299,299),(300,300))
    #########
    println("aligning")
    sfn = Fields.ScalarFieldNode{2}([sf1,sf2,vf1])
    println("geometry of sf1")
    println(Fields.geometry(sf1))
    println("geometry of sf2")
    println(Fields.geometry(sf2))
    println("geometry of vf1")
    println(Fields.geometry(vf1))
    println("geometry of sfn")
    println(Fields.geometry(sfn))
    println("align sfn")
    Fields.align_field_tree!(sfn)
    println("geometry of sf1")
    println(Fields.geometry(sf1))
    println("geometry of sf2")
    println(Fields.geometry(sf2))
    println("geometry of vf1")
    println(Fields.geometry(vf1))
    println("geometry of sfn")
    println(Fields.geometry(sfn))
    ##########3D test
    println("######################################")
    println("3D align test")
    sf1 = Fields.zero(ScalarField{Float64,3},(1001,1001,11),(0,0,0),(1000,1000,10))
    sf2 = Fields.zero(ScalarField{Float64,3},(101,101,11),(299,299,0),(300,300,10))
    vf1 = Fields.func2field(VectorField{Complex{Float64},3},(x,y,z)->[sin(x)+cos(x)im,1.0im,1.0im],(101,101,11),(299,299,0),(300,300,10))
    #########
    println("aligning")
    sfn = Fields.ScalarFieldNode{3}([sf1,sf2,vf1])
    println("geometry of sf1")
    println(Fields.geometry(sf1))
    println("geometry of sf2")
    println(Fields.geometry(sf2))
    println("geometry of vf1")
    println(Fields.geometry(vf1))
    println(mean(vf1.field))
    println("geometry of sfn")
    println(Fields.geometry(sfn))
    println("align sfn")
    Fields.align_field_tree!(sfn)
    println("geometry of sf1")
    println(Fields.geometry(sf1))
    println("geometry of sf2")
    println(Fields.geometry(sf2))
    println("geometry of vf1")
    println(Fields.geometry(vf1))
    println(mean(vf1.field))
    println("geometry of sfn")
    println(Fields.geometry(sfn))    
end
main()
@time main()
