#1.use cubic spline interpolation
#2.use bicubic interpolation 4x4 pixels
@inbounds function sample{T<:ComplexOrFloat}(f::ScalarField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    fpos = collect(f.position)
    fres = collect(f.res)
    fsize = collect(f.size)
    rel_pos = pos-fpos
    pidx = round(Int64,div(rel_pos,fres)+1)
    #2D
    output = sub(f.field,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)
#    output = f.field[pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2]
#    output *= f.scaling(t)
    return output*f.scaling(t)
#    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end

@inbounds function sample{T<:ComplexOrFloat}(f::VectorField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    fpos = collect(f.position)
    fres = collect(f.res)
    fsize = collect(f.size)
    rel_pos = pos-fpos
    pidx = round(Int64,div(rel_pos,fres)+1)
    #2D
    output = sub(f.field,:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)
#    output *= f.scaling(t)
    return output*f.scaling(t)
#    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end
#=
function sample_rel{T<:ComplexOrFloat}(f::VectorField{T,2},rel_pos::Vector{Float64},t::Real;order::Integer = 3)
    fres = collect(f.res)
    pidx = round(Int64,div(rel_pos,fres)+1)
    #2D
    output = f.field[:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2]
    output *= f.scaling(t)
    return output
#    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end
=#

function sample(f::VectorFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type = f.typeof
    return sample_inner(output_type,f,pos,t,order=order)
end

function sample_inner{T<:ComplexOrFloat}(::Type{T},f::VectorFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    output = zeros(T,(3,4,4))
    loop_field!(f.fields,pos,t,output)
    scal = f.scaling(t)
    output .*= scal
    return output    
end

function sample(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type = f.typeof
    return sample_inner(output_type,f,pos,t,order=order)
end

function sample_inner{T<:ComplexOrFloat}(::Type{T},f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    vf_arr = filter(x->isa(x,AbstractVectorField),f.fields)
    sf_arr = filter(x->isa(x,AbstractScalarField),f.fields)
    output = zeros(Float64,(4,4))#TODO: Float64 should be T
    vf_output = zeros(T,(3,4,4))
    if ~isempty(sf_arr)
        loop_field!(sf_arr,pos,t,output)
    end
    if ~isempty(vf_arr)
        loop_field!(vf_arr,pos,t,vf_output)
    end    
    vf_output_abs2::Array{Float64,2} = squeeze(sumabs2(vf_output,1),1)
    output = output + vf_output_abs2 #TODO: type not stable
    scal = f.scaling(t)
    output .*= scal
    return output    
end

@inbounds @fastmath function in_field{T<:Field}(f::T,pos::Vector{Float64})
    # fpos = collect(f.position)
    # fsize = collect(f.size)
    # fres = collect(f.res)
    # rel_pos = pos-fpos
    # return all(i->fres[i]<rel_pos[i]<fsize[i]-fres[i],1:length(fsize))
    return all(i->f.res[i]::Float64<(pos[i]-f.position[i]::Float64)<(f.size[i]::Float64-f.res[i]::Float64),1:length(pos))    
end

@inbounds @fastmath function loop_field!{T<:Field,K<:ComplexOrFloat}(f_arr::Vector{T},pos::Vector{Float64},t::Real,output::Array{K};order::Integer = 3)
    ff = filter(x->in_field(x,pos),f_arr)
    
    for vf in ff
        tmp = sample(vf,pos,t;order=order)[:]
        output[:] += tmp[:]
    end
end

function value(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    A = sample(f,pos,t;order=order)
    A = float(real(A)) #TODO: change this
    res = f.res
    @nexprs 2 j->x_j = rem(pos[j],res[j])/res[j]    
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

@fastmath @inbounds function bicubicInterpolate{T<:ComplexOrFloat}(p::Array{T,2},x::Float64,y::Float64)
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
