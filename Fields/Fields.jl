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
global material
global value
include("../fileio.jl")
include("Fields_type.jl")

# include all Fields files.
include("Fields_function.jl")

function init_parallel!(sfn::ScalarFieldNode, mat_sfn::ScalarFieldNode)
    @info "initializing Fields module..."
    @sync begin
        for p = 1:Distributed.nprocs()     #initialize Fields module on each processe
            #set all scaling to t->0.0
            clean_scaling!(sfn)
            clean_scaling!(mat_sfn)
            @async Distributed.remotecall_fetch(init!,p,sfn,mat_sfn)
        end
    end
    eval_scaling!(sfn)
    eval_scaling!(mat_sfn)
end

function init!(sfn::ScalarFieldNode,mat_sfn::ScalarFieldNode)
    global fields
    global material
    fields = 0
    material= 0
    GC.gc()
    fields = copyfield(sfn)
    material = copyfield(mat_sfn)
    GC.gc()
end

function reset!()
    global fields
    fields = 0
    GC.gc()
end

function save_field(f::ScalarFieldNode)
    @info "save field to debug_save_filed.h5!"
    dict = Dict()
    path = f.name
    field_to_dict(f,dict,path)
    dicttoh5("debug_save_field.h5",dict)
end

function field_to_dict(f::T, dict::Dict, path::String) where {T<:FieldNode}
    for ff in f.fields
        field_to_dict(ff,dict,path*"/"*f.name)
    end
end

function field_to_dict(f::T, dict::Dict, path::String) where {T <: Union{VectorField,ScalarField}}
    dict[path*"/"*f.name] = sdata(f.field)
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
