module Fields
using Interpolations
using Devectorize
using Base.Cartesian
#TODO: align use Flat() boundary condition, add more
#TODO: getindex overload for fields.
include("Fields_type.jl")
# variables
export fields, U_total, resolution, size
global fields, U_total
include("Fields_function.jl")
end
