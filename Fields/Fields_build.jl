function mat2sharedarray(filename,variable)
    file = matopen(filename)
    var = read(file, variable) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    var_s = copy_to_sharedarray!(var)
    var = 0
    gc()
    return var_s
end

padding(level) = repeat("    ",level)
function build_field_file(field_config::Dict,level::Integer;name::ASCIIString="field")
    D_type = Dict("Complex" => Complex{Float64}, "Float" => Float64)
    F_type = Dict("ScalarField" => ScalarField, "VectorField" => VectorField)
    filename = ascii(field_config["filename"])
    var = ascii(field_config["variable"])
    ft = F_type[field_config["field-type"]]::DataType
    dt = D_type[ascii(field_config["D-type"])]::DataType
    dim = round(field_config["dim"])::Integer
    pos = convert(Vector{Float64},field_config["pos"])
    sz = convert(Vector{Float64},field_config["size"])
    scaling = eval(parse(field_config["scaling"]))
    Lumberjack.info(padding(level),"building ",field_config["field-type"]," ",name," from file")
    Lumberjack.info(padding(level),"    datatype: $dt")
    Lumberjack.info(padding(level),"    dimension: $dim")
    Lumberjack.info(padding(level),"    position: $pos")
    Lumberjack.info(padding(level),"    size: $sz")
    Lumberjack.info(padding(level),"    scaling: ",field_config["scaling"])
    Lumberjack.info(padding(level),"    reading ",var," from ",filename,"...")
    field_s = mat2sharedarray(filename,var)
    return ft{dt,dim}(field_s,pos,sz,scaling=scaling,name=name)
end

function build_field_zero(field_config::Dict,level::Integer;name::ASCIIString="field")
    F_type = Dict("ScalarField" => ScalarField, "VectorField" => VectorField)
    D_type = Dict("Complex" => Complex{Float64}, "Float" => Float64)
    ft = F_type[field_config["field-type"]]::DataType
    dt = D_type[ascii(field_config["D-type"])]::DataType
    dim = round(field_config["dim"])::Integer
    res = convert(Vector{Int64},field_config["res"])
    pos = convert(Vector{Float64},field_config["pos"])
    sz = convert(Vector{Float64},field_config["size"])
    scaling = eval(parse(field_config["scaling"]))
    Lumberjack.info(padding(level),"building ",field_config["field-type"]," ",name," type: zero")
    Lumberjack.info(padding(level),"    datatype: $dt")
    Lumberjack.info(padding(level),"    dimension: $dim")
    Lumberjack.info(padding(level),"    resolution: $res")
    Lumberjack.info(padding(level),"    position: $pos")
    Lumberjack.info(padding(level),"    size: $sz")
    Lumberjack.info(padding(level),"    scaling: ",field_config["scaling"])
    return Fields.zero_field(ft{dt,dim},res,pos,sz,scaling=scaling,name=name)
end

function build_field_func(field_config::Dict,level::Integer;name::ASCIIString="field")
    F_type = Dict("ScalarField" => ScalarField, "VectorField" => VectorField)
    D_type = Dict("Complex" => Complex{Float64}, "Float" => Float64)
    ft = F_type[field_config["field-type"]]::DataType
    dt = D_type[ascii(field_config["D-type"])]::DataType
    dim = round(field_config["dim"])::Integer
    res = convert(Vector{Int64},field_config["res"])
    pos = convert(Vector{Float64},field_config["pos"])
    sz = convert(Vector{Float64},field_config["size"])
    func = eval(parse(field_config["func"]))
    scaling = eval(parse(field_config["scaling"]))
    Lumberjack.info(padding(level),"building ",field_config["field-type"]," ",name," type: func")
    Lumberjack.info(padding(level),"    datatype: $dt")
    Lumberjack.info(padding(level),"    dimension: $dim")
    Lumberjack.info(padding(level),"    function: ",field_config["func"])
    Lumberjack.info(padding(level),"    resolution: $res")
    Lumberjack.info(padding(level),"    position: $pos")
    Lumberjack.info(padding(level),"    size: $sz")
    Lumberjack.info(padding(level),"    scaling: ",field_config["scaling"])

    return Fields.func2field(ft{dt,dim},func,res,pos,sz,scaling=scaling,name=name)
end

function build_field(field_config::Dict,level::Integer;name::ASCIIString = "field")
    #TODO: if cant find scaling, use default.
    #TODO: use polymorphism to reduce use of "if"
    if field_config["field-type"] == "ScalarFieldNode"
        Lumberjack.info(padding(level),"building ScalarFieldNode ", name)
        Lumberjack.info(padding(level),"scaling:", field_config["scaling"])
        f_arr = map(field_config["fields"])do x
            build_field(x[2],level+1,name=ascii(x[1]))
        end
        f_arr = [promote(f_arr...)...]
        dim = round(field_config["dim"])::Integer
        scaling = eval(parse(field_config["scaling"]))
        return ScalarFieldNode{dim}(f_arr,scaling=scaling,name=name)
    elseif field_config["field-type"] == "VectorFieldNode"
        Lumberjack.info(padding(level),"building VectorFieldNode ",name)
        Lumberjack.info(padding(level),"scaling:", field_config["scaling"])
        f_arr = map(field_config["fields"])do x
            build_field(x[2],level+1,name=ascii(x[1]))
        end
        f_arr = [promote(f_arr...)...]
        dim = round(field_config["dim"])::Integer
        scaling = eval(parse(field_config["scaling"]))
        return VectorFieldNode{dim}(f_arr,scaling=scaling,name=name)
    elseif field_config["field-type"] == "ScalarField" || field_config["field-type"] == "VectorField"
        if field_config["init-type"] == "file"
            return build_field_file(field_config,level,name=name)
        elseif field_config["init-type"] == "zero"
            return build_field_zero(field_config,level,name=name)
        elseif field_config["init-type"] == "func"
            return build_field_func(field_config,level,name=name)
        else
            Lumberjack.error("Unrecognized init-type")
        end
    else
        Lumberjack.error("Unrecognized field-type")
    end
end

function buildAndAlign(field_config::Dict,level::Integer;name::ASCIIString = "field")
    sfn = build_field(field_config,level,name=name)
    Lumberjack.info("aligning fields...")
    align_field_tree!(sfn)
    set_geometry!(sfn)
    set_typeof!(sfn)
    return sfn
end
