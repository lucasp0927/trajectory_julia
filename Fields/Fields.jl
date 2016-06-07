module Fields
using Lumberjack
using Interpolations
using Devectorize
using Base.Cartesian
#using FastAnonymous
#TODO: align use Flat() boundary condition, add more
#TODO: getindex overload for fields.
#TODO: support different resolution for different area
global fields
include("Fields_type.jl")
include("Fields_function.jl")
# variables

function init_parallel!(sfn::ScalarFieldNode)
    Lumberjack.info("initializing Fields module...")
    @sync begin
        for p = 1:nprocs()     #initialize Fields module on each processe
        @async remotecall_fetch(p,init!,sfn)
        end
    end
end

function init!(sfn::ScalarFieldNode)
    global fields
    fields = 0
    gc()
    fields = copyfield(sfn)
    gc()
end

end
