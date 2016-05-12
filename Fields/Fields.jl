module Fields
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
    println("start initialization")
    @sync begin
        @async begin
            for p = 1:nprocs()     #initialize Fields module on each processe
                remotecall_fetch(p,Fields.init!,sfn)
            end
        end
    end
end

function init!(sfn::ScalarFieldNode)
    global fields
    fields = 0
    gc()
    info("copy fields...")
    fields = copyfield(sfn)
    gc()
end

end
