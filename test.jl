push!(LOAD_PATH, pwd())
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
    println(f3.position)
    println(f3.dim)
    println(f3.scaling)
    f1 = Fields.zero(VectorField{Float64,3},(10,10,10),(0,0,0),(10,10,10))
    f2 = Fields.zero(ScalarField{Float64,2},(10,10),(2,2),(10,10))
    # a = Fields.VectorFieldNode([f1,f2])
    # println(map((x)->getfield(x,:position)[1],[f1,f2]))
    # println(Fields.composite(a))
end
main()