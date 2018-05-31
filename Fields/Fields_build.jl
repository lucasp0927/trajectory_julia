function buildAndAlign(field_config::Dict,level::Integer;name::String = "field")
    sfn = build_field(field_config,level,name=name)
    set_typeof!(sfn)
    @info "aligning fields..."
    align_field_tree!(sfn)
    return sfn
end

padding(level) = repeat("----",level)
function build_field_file(field_config::Dict,level::Integer;name::String="field")
    D_type = Dict("Complex" => Complex{Float64}, "Float" => Float64)
    F_type = Dict("ScalarField" => ScalarField, "VectorField" => VectorField)
    filename = ascii(field_config["filename"])
    var = ascii(field_config["variable"])
    ft = F_type[field_config["field-type"]]::Type
    dt = D_type[ascii(field_config["D-type"])]::Type
    dim = round(field_config["dim"])::Integer
    pos = convert(Vector{Float64},field_config["pos"])
    sz = convert(Vector{Float64},field_config["size"])
    scaling_expr = parse(field_config["scaling"])
    @info padding(level),"building "*field_config["field-type"]*" "*name*" from file"
    @info padding(level),"    datatype: $dt"
    @info padding(level),"    dimension: $dim"
    @info padding(level),"    position: $pos"
    @info padding(level),"    size: $sz"
    @info padding(level),"    scaling: ",field_config["scaling"]
    @info padding(level),"    reading ",var," from ",filename,"..."
    field_s = mat2sharedarray(filename,var)
    if ft == ScalarField
        @assert ndims(field_s) == dim
    else
        @assert ndims(field_s) == dim+1
        @assert size(field_s,1) == 3
    end
    return ft{dt,dim}(field_s,pos,sz,scaling_expr=scaling_expr,name=name)
end

function build_field_zero(field_config::Dict,level::Integer;name::String="field")
    F_type = Dict("ScalarField" => ScalarField, "VectorField" => VectorField)
    D_type = Dict("Complex" => Complex{Float64}, "Float" => Float64)
    ft = F_type[field_config["field-type"]]::DataType
    dt = D_type[ascii(field_config["D-type"])]::DataType
    dim = round(field_config["dim"])::Integer
    res = convert(Vector{Int64},field_config["res"])
    pos = convert(Vector{Float64},field_config["pos"])
    sz = convert(Vector{Float64},field_config["size"])
    scaling_expr = parse(field_config["scaling"])
    @info padding(level)*"building "*field_config["field-type"]*" "*name*" type: zero"
    @info padding(level)*"    datatype: $dt"
    @info padding(level)*"    dimension: $dim"
    @info padding(level)*"    resolution: $res"
    @info padding(level)*"    position: $pos"
    @info padding(level)*"    size: $sz"
    @info padding(level)*"    scaling: "*field_config["scaling"]
    return Fields.zero_field(ft{dt,dim},res,pos,sz,scaling_expr=scaling_expr,name=name)
#    return Fields.zero_field(ft{dt,dim},res,pos,sz,name=name)
end

function build_field_func(field_config::Dict,level::Integer;name::String="field")
    F_type = Dict("ScalarField" => ScalarField, "VectorField" => VectorField)
    D_type = Dict("Complex" => Complex{Float64}, "Float" => Float64)
    ft = F_type[field_config["field-type"]]::DataType
    dt = D_type[ascii(field_config["D-type"])]::DataType
    dim = round(field_config["dim"])::Integer
    res = convert(Vector{Int64},field_config["res"])
    pos = convert(Vector{Float64},field_config["pos"])
    sz = convert(Vector{Float64},field_config["size"])
    func = eval(parse(field_config["func"]))
    scaling_expr = eval(parse(field_config["scaling"]))
    @info padding(level)*"building "*field_config["field-type"]*" "*name*" type: func"
    @info padding(level)*"    datatype: $dt"
    @info padding(level)*"    dimension: $dim"
    @info padding(level)*"    function: "*field_config["func"]
    @info padding(level)*"    resolution: $res"
    @info padding(level)*"    position: $pos"
    @info padding(level)*"    size: $sz"
    @info padding(level)*"    scaling: "*field_config["scaling"]
    return Fields.func2field(ft{dt,dim},func,res,pos,sz,scaling_expr=scaling_expr,name=name)
end

function build_field(field_config::Dict,level::Integer;name::String = "field")
    #TODO: if cant find scaling, use default.
    #TODO: use polymorphism to reduce use of "if"
    if field_config["field-type"] == "ScalarFieldNode"
        @info padding(level)*"building ScalarFieldNode "* name
        @info padding(level)*"scaling:"* field_config["scaling"]
#        f_arr = Array(Any,length(field_config["fields"]))
        f_arr = Array{Any}(length(field_config["fields"]))
        for (i,x) in enumerate(field_config["fields"])
            f_arr[i] = build_field(x[2],level+1,name=ascii(x[1]))
        end
        f_arr = [promote(f_arr...)...]
        dim = round(field_config["dim"])::Integer
        scaling_expr = parse(field_config["scaling"])
        return ScalarFieldNode{dim}(f_arr,scaling_expr=scaling_expr,name=name)
        #return ScalarFieldNode{dim}(f_arr,name=name)
    elseif field_config["field-type"] == "VectorFieldNode"
        @info padding(level)*"building VectorFieldNode "*name
        @info padding(level)*"scaling:"* field_config["scaling"]
        f_arr = Array{Any}(length(field_config["fields"]))
        for (i,x) in enumerate(field_config["fields"])
            f_arr[i] = build_field(x[2],level+1,name=ascii(x[1]))
        end
        f_arr = [promote(f_arr...)...]
        dim = round(field_config["dim"])::Integer
        scaling_expr = parse(field_config["scaling"])
        return VectorFieldNode{dim}(f_arr,scaling_expr=scaling_expr,name=name)
    elseif field_config["field-type"] == "ScalarField" || field_config["field-type"] == "VectorField"
        if field_config["init-type"] == "file"
            return build_field_file(field_config,level,name=name)
        elseif field_config["init-type"] == "zero"
            return build_field_zero(field_config,level,name=name)
        elseif field_config["init-type"] == "func"
            return build_field_func(field_config,level,name=name)
        else
            err("Unrecognized init-type")
        end
    else
        err("Unrecognized field-type")
    end
end
