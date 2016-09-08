addprocs(4)
push!(LOAD_PATH, "./Fields")
using Fields

function test()
    zero_f = Fields.zero_field(ScalarField{Float64,2},[10,10],[0.0,0.0],[100.0,100.0])
    sfn = Fields.ScalarFieldNode{2}([zero_f])
    return sfn
end

function main()
    sfn = test()
    println("test")
    Fields.init_parallel!(sfn)
end
main()
