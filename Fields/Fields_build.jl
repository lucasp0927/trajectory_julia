function buildAndAlign(field_config::Dict,level::Integer;name::String = "field")
    sfn = build_field(field_config,level,name=name)
    set_typeof!(sfn)
    Lumberjack.info("aligning fields...")
    align_field_tree!(sfn)
    return sfn
end

padding(level) = repeat("    ",level)
function build_field_file(field_config::Dict,level::Integer;name::String="field")
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
    if ft == ScalarField
        @assert ndims(field_s) == dim
    else
        @assert ndims(field_s) == dim+1
        @assert size(field_s,1) == 3
    end
    return ft{dt,dim}(field_s,pos,sz,scaling=scaling,name=name)
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

function build_field(field_config::Dict,level::Integer;name::String = "field")
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

function test_build()
    #2D vector
    fconfig = Dict{Any,Any}("scaling"=>"t->1.0","field-type"=>"ScalarFieldNode","fields"=>Dict{Any,Any}("lattice-beam-wrap"=>Dict{Any,Any}("scaling"=>"t->-1.42884e-5","field-type"=>"ScalarFieldNode","fields"=>Dict{Any,Any}("lattice-beam"=>Dict{Any,Any}("scaling"=>"t->1.0+0.0im","field-type"=>"VectorFieldNode","fields"=>Dict{Any,Any}("right-beam"=>Dict{Any,Any}("scaling"=>"t->1.0+0.0im","field-type"=>"VectorField","D-type"=>"Complex","init-type"=>"file","dim"=>2,"filename"=>"/home/lucaspeng/Desktop/trajectory_julia/Exy_60u_0_right.mat","size"=>Any[70000,50000],"pos"=>Any[0.0,0.0],"variable"=>"field"),"left-beam"=>Dict{Any,Any}("scaling"=>"t->exp(2*pi*im*0.8*t)","field-type"=>"VectorField","D-type"=>"Complex","init-type"=>"file","dim"=>2,"filename"=>"/home/lucaspeng/Desktop/trajectory_julia/Exy_60u_0_left.mat","size"=>Any[70000,50000],"pos"=>Any[0.0,0.0],"variable"=>"field")),"dim"=>2)),"dim"=>2)),"dim"=>2)
    buildAndAlign(fconfig,0,name = "field")
    #3D
end
