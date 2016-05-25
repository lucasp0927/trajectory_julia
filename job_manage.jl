function single_scan_scaling(config::Dict,sfn::ScalarFieldNode)
    range = config["range"]
    field_name = config["field"]
    scaling = config["scaling"]
    i = 1
    println(scaling,eval(parse(scaling))(0))
    #scaling = eval(parse(field_config["scaling"]))
end
