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
    Fields.initialize((100,100),(100,100))
    ########2D test
    vf1 = Fields.zero(VectorField{Float64,2},(10,10),(2,2),(9,9),scaling=x->sin(x))
    vf2 = Fields.zero(VectorField{Float64,2},(10,10),(1,1),(9,9))
    vf3 = Fields.func2field(VectorField{Float64,2},(x,y)->[x,y,0],(2,2),(0,0),(2,2))
    vf4 = Fields.func2field(VectorField{Float64,2},(x,y)->[x,y,0],(2,2),(2,2),(4,12))
    #########################################
    vfn = Fields.VectorFieldNode{2}([vf3,vf4])
    vfn2 = Fields.VectorFieldNode{2}([vf1,vf2])
    println(Fields.geometry(vfn))
    println(Fields.geometry(vfn2))
    vfn3 = Fields.VectorFieldNode{2}([vfn,vfn2])
    println(Fields.geometry(vfn3))
    ###############
    sf1 = Fields.zero(ScalarField{Float64,2},(10,10),(2,2),(9,9))
    sf2 = Fields.zero(ScalarField{Float64,2},(10,10),(1,1),(9,9))
    sf3 = Fields.func2field(ScalarField{Float64,2},(x,y)->x,(2,2),(0,0),(2,2))
    sf4 = Fields.func2field(ScalarField{Float64,2},(x,y)->x,(2,2),(2,2),(4,12))
    ############
    sfn = Fields.ScalarFieldNode{2}([sf3,sf4])
    sfn2 = Fields.ScalarFieldNode{2}([sf1,sf2])
    println(Fields.geometry(sfn))
    println(Fields.geometry(sfn2))
    sfn3 = Fields.ScalarFieldNode{2}([sfn,sfn2])
    println(Fields.geometry(sfn3))
    ##########3D test
    sf1 = Fields.zero(ScalarField{Float64,3},(10,10,10),(2,2,0),(9,9,0.9))
    sf2 = Fields.zero(ScalarField{Float64,3},(10,10,10),(1,1,0),(9,9,0.9))
    sf3 = Fields.func2field(ScalarField{Float64,3},(x,y,z)->x,(2,2,2),(0,0,0),(2,2,1))
    sf4 = Fields.func2field(ScalarField{Float64,3},(x,y,z)->x,(2,2,2),(2,2,2),(4,12,2))
    #########
    sfn = Fields.ScalarFieldNode{3}([sf3,sf4])
    sfn2 = Fields.ScalarFieldNode{3}([sf1,sf2])
    println(Fields.geometry(sfn))
    println(Fields.geometry(sfn2))
    sfn3 = Fields.ScalarFieldNode{3}([sfn,sf1,sf2])
    println(Fields.geometry(sfn3))
    Fields.composite(sfn3)
end
main()
@time main()

