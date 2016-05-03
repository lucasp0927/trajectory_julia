#1.use cubic spline interpolation
#2.use bicubic interpolation 4x4 pixels
function sample{T<:ComplexOrFloat}(f::ScalarField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    fpos = collect(f.position)
    fres = collect(f.res)
    fsize = collect(f.size)
    rel_pos = pos-fpos
    pidx = round(Int64,div(rel_pos,fres)+1)
    #2D
    output = f.field[pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2]
    output *= f.scaling(t)
    return output
#    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end

function sample{T<:ComplexOrFloat}(f::VectorField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    fpos = collect(f.position)
    fres = collect(f.res)
    fsize = collect(f.size)
    rel_pos = pos-fpos
    pidx = round(Int64,div(rel_pos,fres)+1)
    #2D
    output = f.field[:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2]
    output *= f.scaling(t)
    return output
#    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end

function sample_rel{T<:ComplexOrFloat}(f::VectorField{T,2},rel_pos::Vector{Float64},t::Real;order::Integer = 3)
    fres = collect(f.res)
    pidx = round(Int64,div(rel_pos,fres)+1)
    #2D
    output = f.field[:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2]
    output *= f.scaling(t)
    return output
#    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end

function sample(f::VectorFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type = typeoffield(f)
    output = zeros(output_type,(3,4,4))
    for vf in f.fields
        vf_geo = geometry(vf)
        fpos = collect(vf_geo["pos"])
        fsize = collect(vf_geo["size"])
        fres = collect(vf_geo["res"])
        rel_pos = pos-fpos
        if all(i->fres[i]<rel_pos[i]<fsize[i]-fres[i],1:length(fsize))
            output += sample(vf,pos,t;order=order)
            # if typeof(vf) <: VectorField
            #     println("vectorfield")
            #     output += sample_rel(vf,rel_pos,t;order = order)
            # else
            #     println("vectorfield node")
            #     output += sample(vf,pos,t;order = order)
            # end
        end
    end
    output *= f.scaling(t)
    return output
end


function sample(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type = typeoffield(f)
    output = zeros(output_type,(4,4))
    vf_arr = filter(x->typeof(x)<:AbstractVectorField,f.fields)
    sf_arr = filter(x->typeof(x)<:AbstractScalarField,f.fields)
    vf_output = zeros(output_type,(3,4,4))
    for sf in sf_arr
        sf_geo = geometry(sf)
        fpos = collect(sf_geo["pos"])
        fsize = collect(sf_geo["size"])
        fres = collect(sf_geo["res"])
        rel_pos = pos-fpos
        if all(i->fres[i]<rel_pos[i]<fsize[i]-fres[i],1:length(fsize))
            output += sample(sf,pos,t;order=order)
            # if typeof(vf) <: VectorField
            #     println("vectorfield")
            #     output += sample_rel(vf,rel_pos,t;order = order)
            # else
            #     println("vectorfield node")
            #     output += sample(vf,pos,t;order = order)
            # end
        end
    end
    for vf in vf_arr
        vf_geo = geometry(vf)
        fpos = collect(vf_geo["pos"])
        fsize = collect(vf_geo["size"])
        fres = collect(vf_geo["res"])
        rel_pos = pos-fpos
        if all(i->fres[i]<rel_pos[i]<fsize[i]-fres[i],1:length(fsize))
            vf_output += sample(vf,pos,t;order=order)
            # if typeof(vf) <: VectorField
            #     println("vectorfield")
            #     output += sample_rel(vf,rel_pos,t;order = order)
            # else
            #     println("vectorfield node")
            #     output += sample(vf,pos,t;order = order)
            # end
        end
    end
    vf_output_abs2::Array{Float64,2} = squeeze(sumabs2(vf_output,1),1)
    output = output + vf_output_abs2 #TODO: type not stable
    output *= f.scaling(t)
    return output
end

function value(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    A = sample(f,pos,t;order=order)
    A = float(real(A)) #TODO: change this
    geo = geometry(f)
    res = collect(geo["res"])
    x_1 = rem(pos[1],res[1])/res[1]
    x_2 = rem(pos[2],res[2])/res[2]
    itp = interpolate(A, BSpline(Cubic(Line())), OnGrid())
    return itp[2.0+x_1,2.0+x_2]
end
