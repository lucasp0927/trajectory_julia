push!(LOAD_PATH, "./Fields")
using Fields
using MAT
include("parse.jl")
fields_config,trajsolver_config = parse_config("config.yml",false)
field_config = fields_config["field"]
file = matopen("D2_TE.mat")
var = read(file, "gm")
sfn = Fields.build_field(field_config,"field",0,true)
#scaling_fun = eval(parse(scaling))
