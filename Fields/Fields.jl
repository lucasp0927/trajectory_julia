module Fields
#TODO: resample field, for faster interpolation
include("Fields_type.jl")
# variables
export fields, U_total, resolution, size
global fields, U_total
include("Fields_function.jl")
end
