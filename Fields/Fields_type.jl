# type definitions
export ComplexOrFloat
export Field, AbstractVectorField, AbstractScalarField, VectorField, ScalarField, VectorFieldNode, ScalarFieldNode, FieldNode
ComplexOrFloat = Union{Complex{Float64},Float64}
abstract Field
abstract AbstractVectorField <: Field
abstract AbstractScalarField <: Field
#TODO: right now sample is only for 2d, also the ScalarFieldNode.sample datatype is Float64
type VectorField{T <: ComplexOrFloat, N} <: AbstractVectorField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
    sample::Array{T,3}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Complex{Float64}
    name::ASCIIString
    function VectorField(f::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling::Function =  t->1.0, name::ASCIIString = "VectorField")
        res = sz./(collect(size(f))[2:N+1]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)-1?new(f,pos,sz,res,scaling,N,zeros(T,(3,4,4)),[0.0,0.0],[0,0,0,0],zero(Complex{Float64}),name):error("dimension error!")
    end
end

type ScalarField{T <: ComplexOrFloat,N} <: AbstractScalarField
    field::SharedArray{T}
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    scaling::Function
    dim::Integer
    sample::Array{T,2}
    rel_pos::Vector{Float64}
    pidx::Vector{Int64}
    s::Float64
    name::ASCIIString
    function ScalarField(f::SharedArray{T,N},pos::Vector{Float64},sz::Vector{Float64};scaling::Function =  t->1.0, name::ASCIIString = "ScalarField")
        res = sz./(collect(size(f))[1:N]-1)
        @assert all(x->x!=0,res) "zero resolution!"
        length(pos)==length(sz)==N==ndims(f)?new(f,pos,sz,res,scaling,N,zeros(T,(4,4)),[0.0,0.0],[0,0,0,0],zero(Float64),name):error("dimension error!")
    end
end


type VectorFieldNode{N} <: AbstractVectorField
    fields::Vector{AbstractVectorField}
    scaling::Function
    dim::Integer
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    typeof::DataType
    sample::Array{Complex{Float64},3}
    s::Complex{Float64}
    name::ASCIIString
    function VectorFieldNode{T<:AbstractVectorField}(f::Vector{T};scaling::Function  =  t->1.0+0.0im,name::ASCIIString="VectorFieldNode")
        @assert all(x->x.dim==N,f) "dimension error!"
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Complex{Float64},(3,4,4)),zero(Complex{Float64}),name)
    end
end

type ScalarFieldNode{N} <: AbstractScalarField
    fields::Vector{Field}
    scaling::Function
    dim::Integer
    position::Vector{Float64}
    size::Vector{Float64}
    res::Vector{Float64}
    typeof::DataType
    sample::Array{Float64,2}
    vf_sample::Array{Complex{Float64},3}
    s::Float64
    one_vf_flag::Bool
    name::ASCIIString
    function ScalarFieldNode{T<:Field}(f::Vector{T};scaling::Function =  t->1.0, name::ASCIIString="ScalarFieldNode")
        @assert all(x->x.dim==N,f) "dimension error!"
        flag = length(f) == 1 && typeof(f[1])<:AbstractVectorField
        new(f,scaling,N,[],[],[],Complex{Float64},zeros(Float64,(4,4)),zeros(Complex{Float64},(3,4,4)),zero(Float64),flag, name)
    end
end

FieldNode = Union{VectorFieldNode,ScalarFieldNode}

function setfield!{T<:ComplexOrFloat,N}(f::ScalarField{T,N},A::SharedArray{T},pos::Vector{Float64},sz::Array{Float64};scaling::Function = t->1.0)
    res = sz./(collect(size(A))[1:N]-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A) "dimension error!"
    f.field = Array(T,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

function setfield!{T<:ComplexOrFloat,N}(f::VectorField{T,N},A::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling::Function = t->1.0)
    res = sz./(collect(size(A))[2:N+1]-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A)-1 "dimension error!"
    f.field = Array(T,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

function setscaling!(f::Field,scaling::Function)
    f.scaling = scaling
end
#=
function find_field(criteria::Function)
    if fields.
end
function find_field{T<:FieldNode}(criteria::Function,f::T)
    if criteria(f)
        return f
    else
        for ff in f.fields
            find_field(criteria,ff)
        end
    end
end

function find_field{T<:Union{ScalarField,VectorField}}(criteria::Function,f::T)
    if criteria(f)
        return f
    end
end
=#
function copyfield{T<:ComplexOrFloat,N}(f::ScalarField{T,N})
    return ScalarField{T,N}(f.field,f.position,f.size,scaling=f.scaling)
end

function copyfield{T<:ComplexOrFloat,N}(f::VectorField{T,N})
    return VectorField{T,N}(f.field,f.position,f.size,scaling=f.scaling)
end

function copyfield{N}(f::ScalarFieldNode{N})
    f_arr = map(copyfield,f.fields)
    fn = ScalarFieldNode{N}(f_arr,scaling = f.scaling)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end

function copyfield{N}(f::VectorFieldNode{N})
    f_arr = map(copyfield,f.fields)
    fn = VectorFieldNode{N}(f_arr,scaling = f.scaling)
    fn.dim = f.dim
    fn.position = f.position
    fn.size = f.size
    fn.res = f.res
    fn.typeof = f.typeof
    return fn
end

function mat2sharedarray(filename,variable)
    file = matopen(filename)
    var = read(file, variable) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    var_s = copy_to_sharedarray!(var)
    return var_s
end

include("../fileio.jl")
padding(level) = repeat("    ",level)
function build_field_file(field_config::Dict,level::Integer,verbose::Bool;name::ASCIIString="field")
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
    if verbose
        println(padding(level),"building ",field_config["field-type"]," ",name," from file")
        println(padding(level),"    datatype: ",dt)
        println(padding(level),"    dimension: ",dim)
        println(padding(level),"    position: ",pos)
        println(padding(level),"    size: ",sz)
        println(padding(level),"    scaling: ",field_config["scaling"])
        println(padding(level),"    reading ",var," from ",filename,"...")
    end
    field_s = mat2sharedarray(filename,var)
        return ft{dt,dim}(field_s,pos,sz,scaling=scaling,name=name)
end

function build_field_zero(field_config::Dict,level::Integer,verbose::Bool;name::ASCIIString="field")
    F_type = Dict("ScalarField" => ScalarField, "VectorField" => VectorField)
    D_type = Dict("Complex" => Complex{Float64}, "Float" => Float64)
    ft = F_type[field_config["field-type"]]::DataType
    dt = D_type[ascii(field_config["D-type"])]::DataType
    dim = round(field_config["dim"])::Integer
    res = convert(Vector{Int64},field_config["res"])
    pos = convert(Vector{Float64},field_config["pos"])
    sz = convert(Vector{Float64},field_config["size"])
    scaling = eval(parse(field_config["scaling"]))
    if verbose
        println(padding(level),"building ",field_config["field-type"]," ",name," type: zero")
        println(padding(level),"    datatype: ",dt)
        println(padding(level),"    dimension: ",dim)
        println(padding(level),"    resolution: ",res)
        println(padding(level),"    position: ",pos)
        println(padding(level),"    size: ",sz)
        println(padding(level),"    scaling: ",field_config["scaling"])
    end
    return Fields.zero_field(ft{dt,dim},res,pos,sz,scaling=scaling,name=name)
end

function build_field_func(field_config::Dict,level::Integer,verbose::Bool;name::ASCIIString="field")
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
    if verbose
        println(padding(level),"building ",field_config["field-type"]," ",name," type: func")
        println(padding(level),"    datatype: ",dt)
        println(padding(level),"    dimension: ",dim)
        println(padding(level),"    function: ",field_config["func"])
        println(padding(level),"    resolution: ",res)
        println(padding(level),"    position: ",pos)
        println(padding(level),"    size: ",sz)
        println(padding(level),"    scaling: ",field_config["scaling"])
    end
    return Fields.func2field(ft{dt,dim},func,res,pos,sz,scaling=scaling,name=name)
end

function build_field(field_config::Dict,level::Integer,verbose::Bool;name::ASCIIString = "field")
    #TODO: if cant find scaling, use default.
    if field_config["field-type"] == "ScalarFieldNode"
        if verbose
            println(padding(level),"building ScalarFieldNode ", name)
            println(padding(level),"scaling:", field_config["scaling"])
        end
        f_arr = map(keys(field_config["fields"]))do k
            build_field(field_config["fields"][k],level+1,verbose,name=ascii(k))
        end
        f_arr = [promote(f_arr...)...]
        dim = round(field_config["dim"])::Integer
        scaling = eval(parse(field_config["scaling"]))
        return ScalarFieldNode{dim}(f_arr,scaling=scaling,name=name)
    elseif field_config["field-type"] == "VectorFieldNode"
        if verbose
            println(padding(level),"building VectorFieldNode ",name)
            println(padding(level),"scaling:", field_config["scaling"])
        end
        f_arr = map(keys(field_config["fields"]))do k
            build_field(field_config["fields"][k],level+1,verbose,name=ascii(k))
        end
        f_arr = [promote(f_arr...)...]
        dim = round(field_config["dim"])::Integer
        scaling = eval(parse(field_config["scaling"]))
        return VectorFieldNode{dim}(f_arr,scaling=scaling,name=name)
    elseif field_config["field-type"] == "ScalarField" || field_config["field-type"] == "VectorField"
        if field_config["init-type"] == "file"
            return build_field_file(field_config,level,verbose,name=name)
        elseif field_config["init-type"] == "zero"
            return build_field_zero(field_config,level,verbose,name=name)
        elseif field_config["init-type"] == "func"
            return build_field_func(field_config,level,verbose,name=name)
        else
            error("Unrecognized init-type")
        end
    else
        error("Unrecognized field-type")
    end
end
