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
    output_type = f.typeof
    output = zeros(output_type,(3,4,4))
    loop_field!(f.fields,pos,t,output)
#    println(mean(output))
#=
    for vf in f.fields
#        vf_geo = geometry(vf)
        fpos = collect(vf.position)
        fsize = collect(vf.size)
        fres = collect(vf.res)
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
=#
    output *= f.scaling(t)
    return output
end

function loop_field!{T<:AbstractVectorField,K<:ComplexOrFloat}(f_arr::Vector{T},pos::Vector{Float64},t::Real,output::Array{K,3};order::Integer = 3)

    for vf in f_arr
#        vf_geo = geometry(vf)
        fpos = collect(vf.position)
        fsize = collect(vf.size)
        fres = collect(vf.res)
        rel_pos = pos-fpos
        if all(i->fres[i]<rel_pos[i]<fsize[i]-fres[i],1:length(fsize))
            #fill!(output,one(K))
            output[:,:] += sample(vf,pos,t;order=order)[:,:]
            # if typeof(vf) <: VectorField
            #     println("vectorfield")
            #     output += sample_rel(vf,rel_pos,t;order = order)
            # else
            #     println("vectorfield node")
            #     output += sample(vf,pos,t;order = order)
            # end
        end
    end
end

function sample(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type = f.typeof
    output = zeros(output_type,(4,4))
    vf_arr = filter(x->typeof(x)<:AbstractVectorField,f.fields)
    sf_arr = filter(x->typeof(x)<:AbstractScalarField,f.fields)
    vf_output = zeros(output_type,(3,4,4))
    for sf in sf_arr
        fpos = collect(sf.position)
        fsize = collect(sf.size)
        fres = collect(sf.res)
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
        fpos = collect(vf.position)
        fsize = collect(vf.size)
        fres = collect(vf.res)
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
    res = f.res
    x_1 = rem(pos[1],res[1])/res[1]
    x_2 = rem(pos[2],res[2])/res[2]
    return itp_bicubic(A,(x_1,x_2))
#    return itp_spline(A,(2.0+x_1,2.0+x_2))
end

function itp_spline(A::Array,pos::Tuple)
    itp = interpolate(A, BSpline(Cubic(Line())), OnGrid())
    return itp[pos...]
end

function itp_bicubic(A::Array,pos::Tuple)
    return bicubicInterpolate(A,pos[1],pos[2])
end

@fastmath @inbounds function cubicInterpolate{T<:ComplexOrFloat}(p::Array{T,1},x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
    return p[2] + 0.5 * x*(p[3] - p[1] + x*(2.0*p[1] - 5.0*p[2] + 4.0*p[3] - p[4] + x*(3.0*(p[2] - p[3]) + p[4] - p[1])));
end

function bicubicInterpolate{T<:ComplexOrFloat}(p::Array{T,2},x::Float64,y::Float64)
#    @assert size(p) == (4,4) "wrong size"
    arr = Array(T,4)
    @nexprs 4 j->(arr[j] = cubicInterpolate(p[:,j], x);)
    # arr[1] = cubicInterpolate(vec(p[1,:]), y);
    # arr[2] = cubicInterpolate(vec(p[2,:]), y);
    # arr[3] = cubicInterpolate(vec(p[3,:]), y);
    # arr[4] = cubicInterpolate(vec(p[4,:]), y);
    return cubicInterpolate(arr, y);
end
#=
double cubicInterpolate (double p[4], double x) {
return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}

double bicubicInterpolate (double p[4][4], double x, double y) {
double arr[4];
arr[0] = cubicInterpolate(p[0], y);
arr[1] = cubicInterpolate(p[1], y);
arr[2] = cubicInterpolate(p[2], y);
arr[3] = cubicInterpolate(p[3], y);
return cubicInterpolate(arr, x);
}

double tricubicInterpolate (double p[4][4][4], double x, double y, double z) {
double arr[4];
arr[0] = bicubicInterpolate(p[0], y, z);
arr[1] = bicubicInterpolate(p[1], y, z);
arr[2] = bicubicInterpolate(p[2], y, z);
arr[3] = bicubicInterpolate(p[3], y, z);
return cubicInterpolate(arr, x);
}
=#
