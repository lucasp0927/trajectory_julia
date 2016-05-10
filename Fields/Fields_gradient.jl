#1.use cubic spline interpolation
#2.use bicubic interpolation 4x4 pixels
@inbounds function sample{T<:ComplexOrFloat}(f::ScalarField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
    pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64})+1)
    return (sub(f.field,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)*f.scaling(t)::ComplexOrFloat)::Array{T,2}
end

@inbounds function sample{T<:ComplexOrFloat}(f::VectorField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
    pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64})+1)
#    return (f.field[:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2])::Array{T,3}
    return (sub(f.field,:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)*f.scaling(t)::ComplexOrFloat)::Array{T,3}
end

function sample(f::VectorFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type::DataType = f.typeof
    return sample_inner(output_type,f,pos,t,order=order)::Array{output_type,3}
end

function sample_inner{T<:ComplexOrFloat}(::Type{T},f::VectorFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    output::Array{T,3} = zeros(T,(3,4,4))
    loop_field!(f.fields,pos,t,output)
    return (output*f.scaling(t)::ComplexOrFloat)::Array{T,3}
end

function sample(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type::DataType = f.typeof
    return sample_inner(output_type,f,pos,t,order=order)::Array{Float64,2}
end

function sample_inner{T<:ComplexOrFloat}(::Type{T},f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    sfoutput::Array{Float64,2} = zeros(Float64,(4,4))#TODO: Float64 should be T
    vfoutput::Array{T,3} = zeros(T,(3,4,4))
    for ff in f.fields
        if in_field(ff,pos)
            sum_field(ff,pos,t,sfoutput,vfoutput,order = order)
        end
    end
    vfoutput_abs2::Array{Float64,2} = squeeze(sumabs2(vfoutput,1),1)
    sfoutput .+= vfoutput_abs2::Array{Float64,2}
    return (sfoutput*f.scaling(t)::ComplexOrFloat)::Array{Float64,2}
#    vf_arr::Vector{AbstractVectorField} = filter(x->isa(x,AbstractVectorField),f.fields)
#    sf_arr::Vector{AbstractScalarField} = filter(x->isa(x,AbstractScalarField),f.fields)

    # if ~isempty(sf_arr)
    #     loop_field!(sf_arr,pos,t,output)
    # end
    # if ~isempty(vf_arr)
    #     loop_field!(vf_arr,pos,t,vf_output)
    # end
end

@inbounds @fastmath function in_field{T<:Field}(f::T,pos::Vector{Float64})
    # fpos = collect(f.position)
    # fsize = collect(f.size)
    # fres = collect(f.res)
    # rel_pos = pos-fpos
    # return all(i->fres[i]<rel_pos[i]<fsize[i]-fres[i],1:length(fsize))
    return (f.res[1]::Float64<(pos[1]-f.position[1]::Float64)<(f.size[1]::Float64-f.res[1]::Float64) && f.res[2]::Float64<(pos[2]-f.position[2]::Float64)<(f.size[2]::Float64-f.res[2]::Float64))
end

function sum_field{T<:AbstractVectorField,K<:ComplexOrFloat}(f::T,pos::Vector{Float64},t::Real,sf_output::Array{Float64,2},vf_output::Array{K,3};order::Integer = 3)
    tmp = sample(f,pos,t;order=order)::Array{K}
    vf_output[:] += tmp[:]
end
function sum_field{T<:AbstractScalarField,K<:ComplexOrFloat}(f::T,pos::Vector{Float64},t::Real,sf_output::Array{Float64,2},vf_output::Array{K,3};order::Integer = 3)
    tmp = sample(f,pos,t;order=order)::Array{Float64}
    sf_output[:] += tmp[:]
end


@inbounds @fastmath function loop_field!{T<:Field,K<:ComplexOrFloat}(f_arr::Vector{T},pos::Vector{Float64},t::Real,output::Array{K};order::Integer = 3)
    ff = filter(x->in_field(x,pos),f_arr)
    for vf in ff
        tmp = sample(vf,pos,t;order=order)::Array{K}
        output[:] += tmp[:]
    end
end

function value(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    A = sample(f,pos,t;order=order)::Array{Float64,2}
#    A = float(real(A)) #TODO: change this
    res = f.res::Vector{Float64}
    @nexprs 2 j->x_j = rem(pos[j],res[j])/res[j]
    return itp_bicubic(A,[x_1,x_2])
#    return itp_spline(A,(2.0+x_1,2.0+x_2))
end

function gradient(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
#        A = Array(Float64,4,4)
        A = sample(f,pos,t;order=order)::Array{Float64,2}
        #    A = float(real(A)) #TODO: change this
        res = f.res::Vector{Float64}
        @nexprs 2 j->x_j = rem(pos[j],res[j])/res[j]
        return itp_bicubic_grad(A,[x_1,x_2],res)
    #    return itp_spline(A,(2.0+x_1,2.0+x_2))
end

function itp_spline(A::Array,pos::Vector{Float64})
    itp = interpolate(A, BSpline(Cubic(Line())), OnGrid())
    return itp[pos...]
end

function itp_bicubic(A::Array,pos::Vector{Float64})
    return bicubicInterpolate(A,pos[1],pos[2])
end

function itp_bicubic_grad{T<:ComplexOrFloat}(A::Array{T},pos::Vector{Float64},res::Vector{Float64})
    arr = Array(T,4)
    @nexprs 4 j->(arr[j] = cubicInterpolate_grad(A[:,j], pos[1]);)
    gx = cubicInterpolate(arr, pos[2])/res[1];
    @nexprs 4 j->(arr[j] = cubicInterpolate_grad(vec(A[j,:]), pos[2]);)
    gy = cubicInterpolate(arr, pos[1])/res[2];
    return [gx,gy]
end

@fastmath @inbounds function cubicInterpolate{T<:ComplexOrFloat}(p::Array{T,1},x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
#    return p[2]+0.5*(p[3]-p[1])*x+(p[1]-2.5*p[2]+2.0*p[3]-0.5*p[4])*x^2 +0.5*(-p[1]+3.0*(p[2]-p[3])+p[4])*x^3
    return p[2] + 0.5 * x*(p[3] - p[1] + x*(2.0*p[1] - 5.0*p[2] + 4.0*p[3] - p[4] + x*(3.0*(p[2] - p[3]) + p[4] - p[1])));
end

@fastmath @inbounds function cubicInterpolate_grad{T<:ComplexOrFloat}(p::Array{T,1},x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
    return -0.5p[1]+0.5p[3]+2p[1]*x-5p[2]*x+4p[3]*x-p[4]*x-1.5(p[1]-3p[2]+3p[3]-p[4])*x^2
end


@generated function bicubicInterpolate{T<:ComplexOrFloat}(p::Array{T,2},x::Vector{Float64})
    quote
        #    @assert size(p) == (4,4) "wrong size"
        arr = $(Array(T,4))
        @nexprs 4 j->(arr[j] = cubicInterpolate(p[:,j], x[1]);)
        # arr[1] = cubicInterpolate(vec(p[1,:]), y);
        # arr[2] = cubicInterpolate(vec(p[2,:]), y);
        # arr[3] = cubicInterpolate(vec(p[3,:]), y);
        # arr[4] = cubicInterpolate(vec(p[4,:]), y);
        return cubicInterpolate(arr, x[2]);
    end
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
