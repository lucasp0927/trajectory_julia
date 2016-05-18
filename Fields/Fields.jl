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
    println("start initialization Fields module")
    @sync begin
        @async begin
            for p = 2:nprocs()     #initialize Fields module on each processe
                remotecall_fetch(p,init!,sfn)
            end
        end
    end
end

function init!(sfn::ScalarFieldNode)
    println("initialize Field module on process ", myid())
    global fields
    fields = 0
    gc()
    fields = copyfield(sfn)
    gc()
end

end
