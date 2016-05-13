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

@fastmath @inbounds function cubicInterpolate{T<:ComplexOrFloat}(p::SubArray{T,1},x::Float64)
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

@fastmath @inbounds function cubicInterpolate_grad{T<:ComplexOrFloat}(p::SubArray{T,1},x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
    return -0.5p[1]+0.5p[3]+2p[1]*x-5p[2]*x+4p[3]*x-p[4]*x-1.5(p[1]-3p[2]+3p[3]-p[4])*x^2
end


@generated function bicubicInterpolate{T<:ComplexOrFloat}(p::Array{T,2},x::Vector{Float64})
    quote
        arr = $(Array(T,4))
        @nexprs 4 j->(arr[j] = cubicInterpolate(slice(p,:,j), x[1]);)
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
