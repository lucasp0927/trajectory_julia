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
    E1S = readfields("/Data/Lucas/lattice_240nm.mat","E1S")
#    @time E1S_f = [(E1S["Ex"][i,j]::Float64,E1S["Ey"][i,j]::Float64,E1S["Ez"][i,j]::Float64) for i=1:size(E1S["Ex"],1), j = 1:size(E1S["Ex"],2)]::VectorFieldArray
    @time E1S_f2 = map(Complex{Float64},cat(3,E1S["Ex"],E1S["Ey"],E1S["Ez"]))
    f2 = VectorField{Complex{Float64}}(E1S_f2,(0,0),(15*6667,15*3334),1)
    f = ScalarField{Complex{Float64}}(E1S["Ex"],(0,0),(15*6667,15*3334),1)
    println(typeof(f2) <: VectorField)
    println(size(f.field))
    println(f.position)
    println(f.size)
end
main()