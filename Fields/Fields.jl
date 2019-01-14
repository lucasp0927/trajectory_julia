module Fields
using Distributed
using Interpolations
using Base.Cartesian
using Test
using Logging
#TODO: align use Flat() boundary condition, add more
#TODO: getindex overload for fields.
#TODO: support different resolution for different area
global fields
global value
include("../fileio.jl")
include("Fields_type.jl")

# include all Fields files.
include("Fields_function.jl")

function init_parallel!(sfn::ScalarFieldNode)
    @info "initializing Fields module..."
    @sync begin
        for p = 1:Distributed.nprocs()     #initialize Fields module on each processe
            #set all scaling to t->0.0
            clean_scaling!(sfn)
            @async Distributed.remotecall_fetch(init!,p,sfn)
        end
    end
    eval_scaling!(sfn)
end

function init!(sfn::ScalarFieldNode)
    global fields
    fields = 0
    GC.gc()
    fields = copyfield(sfn)
    GC.gc()
end

function reset!()
    global fields
    fields = 0
    GC.gc()
end

function test()
    @info "testing zero_field"
    test_zero_field()
    @info "testing typeof"
    test_typeof()
    @info "testing find"
    test_find()
    @info "testing geometry"
    test_geometry()
    @info "testing align"
    test_align()
#####
#    @info "testing build"
#    test_build()
#####
    @info "testing interpolate"
    test_interpolate()
    @info "testing gradient"
    test_gradient()
end

end
