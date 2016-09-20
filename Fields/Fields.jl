module Fields
using Interpolations
using Devectorize
using Base.Cartesian
using Base.Test
using Logging
#using FastAnonymous
#TODO: align use Flat() boundary condition, add more
#TODO: getindex overload for fields.
#TODO: support different resolution for different area
global fields
global value
include("../fileio.jl")
include("Fields_type.jl")
include("Fields_function.jl")
# variables

function init_parallel!(sfn::ScalarFieldNode)
    info("initializing Fields module...")
    @sync begin
        for p = 1:nprocs()     #initialize Fields module on each processe
            #set all scaling to t->0.0
            clean_scaling!(sfn)
            @async remotecall_fetch(init!,p,sfn)
        end
    end
    eval_scaling!(sfn)
end

function init!(sfn::ScalarFieldNode)
    global fields
    Lumberjack.debug("in init!")
    fields = 0
    gc()
    fields = copyfield(sfn)
    gc()
end

function reset!()
    global fields
    fields = 0
    gc()
end

function test()
    info("testing zero_field")
    test_zero_field()
    info("testing typeof")
    test_typeof()
    info("testing find")
    test_find()
    info("testing geometry")
    test_geometry()
    info("testing align")
    test_align()
#    info("testing build")
#    test_build()
    info("testing interpolate")
    test_interpolate()
    info("testing gradient")
    test_gradient()
end

end
