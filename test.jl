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
    f3 = ScalarField{Float64,3}(zeros(Float64,(10,10,10)),(0,0,0),(10,10,10),scaling = 3.0)
    f1 = Fields.zero(VectorField{Float64,3},(10,10,10),(0,0,0),(10,10,10))
    f2 = Fields.zero(ScalarField{Float64,2},(10,10),(2,2),(10,10))
    f4 = Fields.func2field(ScalarField{Float64,3},(x,y,z)->sin(x),(5,5,2),(0,0,0),(2pi,2pi,1))
    f5 = Fields.func2field(VectorField{Float64,3},(x,y,z)->[x,y,z],(2,2,2),(0,0,0),(2,2,2))
    println(f4.field)
    f6 = Fields.composite(f5)
    println(typeof(f6))
    # a = Fields.VectorFieldNode([f1,f2])
    # println(map((x)->getfield(x,:position)[1],[f1,f2]))
    # println(Fields.composite(a))
end
main()