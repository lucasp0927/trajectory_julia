@generated function itp_bicubic_grad{T<:ComplexOrFloat}(A::Array{T,2},pos::Vector{Float64},res::Vector{Float64})
    quote
        arr = $(Array(T,4))
        grad = $(Array(T,2))
        @nexprs 4 j->(arr[j] = cubicInterpolate_grad(view(A,:,j), pos[1]);)
        grad[1] = cubicInterpolate(arr, pos[2])/res[1];
        @nexprs 4 j->(arr[j] = cubicInterpolate_grad(view(A,j,:), pos[2]);)
        grad[2] = cubicInterpolate(arr, pos[1])/res[2];
        return grad
    end
end

@generated function itp_tricubic_grad{T<:ComplexOrFloat}(A::Array{T,3},pos::Vector{Float64},res::Vector{Float64})
    quote
        arr = $(Array(T,4,4))
        grad = $(Array(T,3))
        for k = 1:4
            @nexprs 4 j->(arr[j,k] = cubicInterpolate_grad(view(A,:,j,k), pos[1]);)
        end
        grad[1] = bicubicInterpolate(arr, pos[2:3])/res[1];
        for k = 1:4
            @nexprs 4 j->(arr[j,k] = cubicInterpolate_grad(view(A,j,:,k), pos[2]);)
        end
        grad[2] = bicubicInterpolate(arr, pos[1:2:3])/res[2];
        for k = 1:4
            @nexprs 4 j->(arr[j,k] = cubicInterpolate_grad(view(A,j,k,:), pos[3]);)
        end
        grad[3] = bicubicInterpolate(arr, pos[1:2])/res[3];
        return grad
    end
end


@fastmath @inbounds function cubicInterpolate{T<:Union{AbstractArray{Float64,1},AbstractArray{Complex{Float64},1}}}(p::T,x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
    #    return p[2]+0.5*(p[3]-p[1])*x+(p[1]-2.5*p[2]+2.0*p[3]-0.5*p[4])*x^2 +0.5*(-p[1]+3.0*(p[2]-p[3])+p[4])*x^3
#    @test_approx_eq (p[2] + 0.5 * x*(p[3] - p[1] + x*(2.0*p[1] - 5.0*p[2] + 4.0*p[3] - p[4] + x*(3.0*(p[2] - p[3]) + p[4] - p[1])))) 0.5*x*(x*(x*(-p[1]+3.0*p[2]-3.0*p[3]+p[4])+2.0*p[1]-5.0*p[2]+4.0*p[3]-p[4])-p[1]+p[3])+p[2]
    return 0.5*x*(x*(x*(-p[1]+3.0*p[2]-3.0*p[3]+p[4])+2.0*p[1]-5.0*p[2]+4.0*p[3]-p[4])-p[1]+p[3])+p[2]
#    return p[2] + 0.5 * x*(p[3] - p[1] + x*(2.0*p[1] - 5.0*p[2] + 4.0*p[3] - p[4] + x*(3.0*(p[2] - p[3]) + p[4] - p[1])))
#    return -0.5*p[1]*x^3+p[1]*x^2-0.5*p[1]*x+1.5*p[2]*x^3-2.5*p[2]*x^2+p[2]-1.5*p[3]*x^3+2*p[3]*x^2+0.5*p[3]*x+0.5*p[4]*x^3-0.5*p[4]*x^2
end
#=
@fastmath @inbounds function cubicInterpolate{T<:ComplexOrFloat}(p::SubArray{T,1},x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
#    return p[2]+0.5*(p[3]-p[1])*x+(p[1]-2.5*p[2]+2.0*p[3]-0.5*p[4])*x^2 +0.5*(-p[1]+3.0*(p[2]-p[3])+p[4])*x^3
    return p[2] + 0.5 * x*(p[3] - p[1] + x*(2.0*p[1] - 5.0*p[2] + 4.0*p[3] - p[4] + x*(3.0*(p[2] - p[3]) + p[4] - p[1])));
end
=#
@fastmath @inbounds function cubicInterpolate_grad{T<:Union{AbstractArray{Float64,1},AbstractArray{Complex{Float64},1}}}(p::T,x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
    return 0.5(p[3]-p[1])+(2p[1]-5p[2]+4p[3]-p[4])*x-1.5(p[1]-3p[2]+3p[3]-p[4])*x^2
end
#=
@fastmath @inbounds function cubicInterpolate_grad{T<:ComplexOrFloat}(p::SubArray{T,1},x::Float64)
    #    @assert length(p) == 4 "wrong length"
    #   p1 p2 p3 p4
    #x= -1 0  1  2
    return 0.5(p[3]-p[1])+(2p[1]-5p[2]+4p[3]-p[4])*x-1.5(p[1]-3p[2]+3p[3]-p[4])*x^2
#    return -0.5p[1]+0.5p[3]+2p[1]*x-5p[2]*x+4p[3]*x-p[4]*x-1.5(p[1]-3p[2]+3p[3]-p[4])*x^2
end
=#
@generated function bicubicInterpolate{T<:ComplexOrFloat}(p::Array{T,2},x::Vector{Float64})
    #p is a 4x4 array
    quote
        arr = $(Array(T,4))
        @nexprs 4 j->(arr[j] = cubicInterpolate(view(p,:,j), x[1]);)
        return cubicInterpolate(arr, x[2]);
    end
end
@generated function bicubicInterpolate{T<:ComplexOrFloat}(p::SubArray{T,2},x::Vector{Float64})
    #p is a 4x4 array
    quote
        arr = $(Array(T,4))
        @nexprs 4 j->(arr[j] = cubicInterpolate(view(p,:,j), x[1]);)
        return cubicInterpolate(arr, x[2]);
    end
end

@generated function tricubicInterpolate{T<:ComplexOrFloat}(p::Array{T,3},x::Vector{Float64})
    # p in a 4x4x4 array
    quote
        arr = $(Array(T,4))
        @nexprs 4 j->(arr[j] = bicubicInterpolate(view(p,:,:,j), x[1:2]);)
        return cubicInterpolate(arr, x[3]);
    end
end

function test_interpolate()
    @info "2D interpolation test"
    x_grid = linspace(0.0,1.0,11)
    y_grid = linspace(0.0,1.0,11)
    # sample = [1.0 1.1 1.2 1.3;
    #           1.1 1.2 1.4 1.5;
    #           1.2 1.5 1.7 1.8;
    #           1.5 1.8 1.9 2.1]
    # sample = sample + 0.2*rand(4,4)
    sample = [1.0+x/8+x^2/10+y/8+y^2/12 for x=1:4,y=1:4]
    sample_itp = interpolate(sample, BSpline(Cubic(Line())), OnGrid())
    sum_err = 0.0
    for x in x_grid,y in y_grid
        my_itp = bicubicInterpolate(sample,[x,y])
        sum_err += abs(sample_itp[2.0+x,2.0+y]-my_itp)/my_itp
    end
    err = sum_err/(length(x_grid)*length(y_grid))
    @info "err: "*string(err)
    @test err < 3e-2
    #test gradient
    sum_err = 0.0
    for x in x_grid,y in y_grid
        my_itp = itp_bicubic_grad(sample,[x,y],[1.0,1.0])
        sum_err += norm(gradient(sample_itp,2.0+x,2.0+y)-my_itp)/norm(my_itp)
    end
    err = sum_err/(length(x_grid)*length(y_grid))
    @info "err: "*string(err)
    @test err < 3e-2

    @info "3D interpolation test"
    x_grid = linspace(0.0,1.0,11)
    y_grid = linspace(0.0,1.0,11)
    z_grid = linspace(0.0,1.0,11)
    sample = [1.0+x/8+y/8+z/8+x^2/10+y^2/10+z^2/10 for x=1:4,y=1:4,z=1:4]
    sample_itp = interpolate(sample, BSpline(Cubic(Line())), OnGrid())
    sum_err = 0.0
    for x in x_grid,y in y_grid, z in z_grid
        my_itp = tricubicInterpolate(sample,[x,y,z])
        sum_err += abs(sample_itp[2.0+x,2.0+y,2.0+z]-my_itp)/my_itp
    end
    err = sum_err/(length(x_grid)*length(y_grid)*length(z_grid))
    @info "err: "*string(err)
    @test err < 3e-2

    sum_err = 0.0
    for x in x_grid,y in y_grid, z in z_grid
        my_itp = itp_tricubic_grad(sample,[x,y,z],[1.0,1.0,1.0])
        sum_err += norm(gradient(sample_itp,2.0+x,2.0+y,2.0+z)-my_itp)/norm(my_itp)
    end
    err = sum_err/(length(x_grid)*length(y_grid)*length(z_grid))
    @info "err: "*string(err)
    @test err < 3e-2
end
