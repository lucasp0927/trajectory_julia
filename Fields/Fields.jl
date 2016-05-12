module Fields
using Interpolations
using Devectorize
using Base.Cartesian
using FastAnonymous
#TODO: align use Flat() boundary condition, add more
#TODO: getindex overload for fields.
#TODO: support different resolution for different area
global fields
include("Fields_type.jl")
include("Fields_function.jl")
# variables

function buildfields_parallel!(rb_field_s::SharedArray,lb_field_s::SharedArray,gm_field_s::SharedArray)
    @sync begin
        @async begin
            for p = 1:nprocs()     #initialize Fields module on each processe
                remotecall_fetch(p,Fields.buildfields!,rb_field_s,lb_field_s,gm_field_s)
            end
        end
    end
end

function buildfields!(rb_field_s::SharedArray,lb_field_s::SharedArray,gm_field_s::SharedArray)
    info("test")
    global fields
    rb = Fields.VectorField{Complex{Float64},2}(rb_field_s,[0.0,0.0],[6666*15.0,3333*15.0],scaling=t-> exp(1.0im*t*2pi))
    lb = Fields.VectorField{Complex{Float64},2}(lb_field_s,[0.0,0.0],[6666*15.0,3333*15.0],scaling= t-> 1.0+0.0im)
    vfn = Fields.VectorFieldNode{2}([lb,rb],scaling= t->1.0+0.0im)
    gm = Fields.ScalarField{Float64,2}(gm_field_s,[(6666*15.0)/2.0-1854.625,(3333*15.0)/2.0-1850.0],[3690.75,3690.75],scaling= t-> 10.0)
    fields = Fields.ScalarFieldNode{2}([vfn,gm])
    info("aligning...")
    align_field_tree!(fields)
    set_geometry!(fields)
    set_typeof!(fields)    
end

end
