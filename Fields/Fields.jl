module Fields
using Interpolations
using Base.Cartesian
#TODO: use Base.Cartesian to write generic loops for both 2d and 3d
#TODO: align use Flat() boundary condition.
#TODO: resample field, for faster interpolation
include("Fields_type.jl")
# variables
export fields, U_total, resolution, size
global fields, U_total
include("Fields_function.jl")
end
