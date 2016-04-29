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
    println(mean(output))
    ####handle output
    ff = map(x->composite(x,t),f.fields)
    vf_arr = filter(x->typeof(x)<:AbstractVectorField,ff)
    sf_arr = filter(x->typeof(x)<:AbstractScalarField,ff)
    #####process scalar fields
    for sf in sf_arr
        sf_geo = geometry(sf)
        sf_pos = collect(sf_geo["pos"])
        sf_sz = collect(sf_geo["size"])
        sf_start_idx = map(x->convert(Int64,x),((sf_pos.-pos)./res)+1)
        sf_end_idx = map(x->convert(Int64,x),((sf_pos.+sf_sz.-pos)./res)+1)
        fill_field!(output,sf.field,sf_start_idx,ndims(sf.field))
        # if N ==2
        #     output[sf_start_idx[1]:sf_end_idx[1],sf_start_idx[2]:sf_end_idx[2]] += sf.field
        # elseif N==3
        #     output[sf_start_idx[1]:sf_end_idx[1],sf_start_idx[2]:sf_end_idx[2],sf_start_idx[3]:sf_end_idx[3]] += sf.field
        # end
        end
    println(mean(output))        
    #ScalarField{Float64,N}(output,pos,new_sz;scaling = t->1.0)
end

@generated function fill_field!(output::Array,field::Array,sf_start_idx::Vector,dim::Int64)
    N::Int64=ndims(field)
    quote
        @nloops $N i field begin
            idx = collect((@ntuple $N i)).+sf_start_idx-1
            output[idx...]=(@nref $N field i)
        end
    end
end
