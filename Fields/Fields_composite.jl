#####composition: for FieldNode object return composite field at time t, replace scaling with x->1.0
function composite{T<:Union{VectorField,ScalarField}}(f::T,t::Real)
    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end

function composite{N}(f::ScalarFieldNode{N},t::Real)
    #remember scaling
    output_type = typeoffield(f)
    geo = geometry(f)
    res = collect(geo["res"])
    pos = collect(geo["pos"])
    sz = collect(geo["size"])
    arr_sz = round(Int64,sz./res+1)
    output = zeros(output_type,(arr_sz...))
    ####handle output
    ff = map(x->composite(x,t),f.fields)
    vf_arr = filter(x->typeof(x)<:AbstractVectorField,ff)
    sf_arr = filter(x->typeof(x)<:AbstractScalarField,ff)
    for sf in sf_arr ##loop over scalar fields
        sf_geo = geometry(sf)
        sf_pos = collect(sf_geo["pos"])
        sf_sz = collect(sf_geo["size"])
        sf_start_idx = map(x->convert(Int64,x),((sf_pos.-pos)./res)+1)
        sf_end_idx = map(x->convert(Int64,x),((sf_pos.+sf_sz.-pos)./res)+1)
        fill_field!(output,sf.field,sf_start_idx,sf_end_idx)
    end
    #TODO: vf_output dont have to be as large as output
    vf_output = zeros(output_type,(3,arr_sz...))
    for vf in vf_arr ##loop over vector fields
        vf_geo = geometry(vf)
        vf_pos = collect(vf_geo["pos"])
        vf_sz = collect(vf_geo["size"])
        vf_start_idx = map(x->convert(Int64,x),((vf_pos.-pos)./res)+1)
        vf_end_idx = map(x->convert(Int64,x),((vf_pos.+vf_sz.-pos)./res)+1)
        ### add all vector fields
        fill_field_vec!(vf_output,vf.field,vf_start_idx,vf_end_idx)
    end
    #abs2 vf_output
    vf_output_abs2::Array{output_type,N} = squeeze(sumabs2(vf_output,1),1)
    output = (output.+vf_output_abs2)*f.scaling(t)
    println(mean(output))
    if sum(imag(output)) == 0
        output = real(output)
        return ScalarField{Float64,N}(output::Array{Float64,N},tuple(pos...),tuple(sz...);scaling = t->1.0)
    else
        return ScalarField{Complex{Float64},N}(output::Array{Complex{Float64},N},tuple(pos...),tuple(sz...);scaling = t->1.0)        
    end
end

@generated function fill_field_vec!{T<:ComplexOrFloat,K<:ComplexOrFloat,N}(output::Array{T,N},field::Array{K,N},vf_start_idx::Vector,vf_end_idx::Vector)
    ex_str = "output[:,"
    for i = 1:N-1
        ex_str = ex_str*"vf_start_idx[$i]:vf_end_idx[$i],"
    end
    parse(ex_str[1:end-1]*"]+=field")
end

@generated function fill_field!{T<:ComplexOrFloat,K<:ComplexOrFloat,N}(output::Array{T,N},field::Array{K,N},sf_start_idx::Vector,sf_end_idx::Vector)
    ex_str = "output["
    for i = 1:N
        ex_str = ex_str*"sf_start_idx[$i]:sf_end_idx[$i],"
    end
    parse(ex_str[1:end-1]*"]+=field")
end

