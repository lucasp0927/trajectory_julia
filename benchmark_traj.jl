@info "benchmark trajectory_julia"
using Distributed
using BenchmarkTools
using Profile
#addprocs(1)

@everywhere push!(LOAD_PATH, "./Fields")
@everywhere push!(LOAD_PATH, "./TrajSolver")
@everywhere push!(LOAD_PATH, "./TrajAnalyzer")
using Fields
using TrajSolver
using TrajAnalyzer
using YAML
include("fileio.jl")

config = YAML.load(open("./test.yml"))
fields_config = config["fields-config"]

@info "Read test field."
sfn = Fields.buildAndAlign(fields_config["field"],0,name=ascii([k for k in keys(fields_config)][1]))

Fields.init_parallel!(sfn)
range = [-150.0, 150.0, -150.0, 150.0]
#range = [-1.0, 1.0, -1.0, 1.0]
#print(Fields.fields.res)
#@btime Fields.composite(range,0.0)
function benchmark_value(range,t)
    xx = range[1]:1.0:range[2]
    yy = range[3]:1.0:range[4]
    output = zeros(Float64,(length(xx),length(yy)))
    for x in enumerate(xx), y in enumerate(yy)
        output[x[1],y[1]] = Fields.value([x[2],y[2]],t)
    end
    return output
end

function benchmark_gradient(range,t)
    xx = range[1]:1.0:range[2]
    yy = range[3]:1.0:range[4]
    output = zeros(Float64,(4,length(xx),length(yy)))
    for x in enumerate(xx), y in enumerate(yy)
        posvel = Float64[x[2],y[2],0.0,0.0]
        grad = zeros(Float64,4)
        Fields.gradient!(t,posvel,grad,Fields.fields)
        output[:,x[1],y[1]] = grad
    end
    return output
end

#benchmark_value(range,0.0)
benchmark_gradient(range,0.0)
@time benchmark_gradient(range,0.0)
Profile.clear()
Profile.clear_malloc_data()
#@profile benchmark_value(range,0.0)
@info "profiling..."
#benchmark_gradient(range,0.0)
@profile benchmark_gradient(range,0.0)
#using ProfileView
#ProfileView.view()
#@time benchmark_gradient(range,0.0)

#@btime Fields.value([10.0,5.0],0.0)
#@btime Fields.gradient!(0.0,[10.0,5.0,0.0,0.0],zeros(Float64,4),sfn)
