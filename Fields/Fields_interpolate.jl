@generated function itp_bicubic_grad(A::Array{T, 2}, pos::Vector{Float64}, res::Vector{Float64}) where T <: ComplexOrFloat
    quote
        arr = $(Array{T}(undef,4))
        grad = $(Array{T}(undef,2))
        @nexprs 4 j->(arr[j] = cubicInterpolate_grad(view(A,:,j), pos[1]);)
        grad[1] = cubicInterpolate(arr, pos[2])/res[1];
        @nexprs 4 j->(arr[j] = cubicInterpolate_grad(view(A,j,:), pos[2]);)
        grad[2] = cubicInterpolate(arr, pos[1])/res[2];
        return grad
    end
end

@generated function itp_tricubic_grad(A::Array{T, 3}, pos::Vector{Float64}, res::Vector{Float64}) where T <: ComplexOrFloat
    quote
        arr = $(Array{T}(undef,4,4))
        grad = $(Array{T}(undef,3))
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


@fastmath @inbounds function cubicInterpolate(p::T, x::Float64) where T <: Union{AbstractArray{Float64, 1}, AbstractArray{Complex{Float64}, 1}}
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
@fastmath @inbounds function cubicInterpolate_grad(p::T, x::Float64) where T <: Union{AbstractArray{Float64, 1}, AbstractArray{Complex{Float64}, 1}}
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
@generated function bicubicInterpolate(p::Array{T, 2}, x::Vector{Float64}) where T <: ComplexOrFloat
    #p is a 4x4 array
    quote
        arr = $(Array{T}(undef,4))
        @nexprs 4 j->(arr[j] = cubicInterpolate(view(p,:,j), x[1]);)
        return cubicInterpolate(arr, x[2]);
    end
end
@generated function bicubicInterpolate(p::SubArray{T, 2}, x::Vector{Float64}) where T <: ComplexOrFloat
    #p is a 4x4 array
    quote
        arr = $(Array{T}(undef,4))
        @nexprs 4 j->(arr[j] = cubicInterpolate(view(p,:,j), x[1]);)
        return cubicInterpolate(arr, x[2]);
    end
end

@generated function tricubicInterpolate(p::Array{T, 3}, x::Vector{Float64}) where T <: ComplexOrFloat
    # p in a 4x4x4 array
    quote
        arr = $(Array{T}(undef,4))
        @nexprs 4 j->(arr[j] = bicubicInterpolate(view(p,:,:,j), x[1:2]);)
        return cubicInterpolate(arr, x[3]);
    end
end

function test_interpolate()
    @testset "test interpolation" begin
        @info "2D interpolation test"
        x_grid = range(0.0,stop=1.0,length=11)
        y_grid = range(0.0,stop=1.0,length=11)
        # sample = [1.0 1.1 1.2 1.3;
        #           1.1 1.2 1.4 1.5;
        #           1.2 1.5 1.7 1.8;
        #           1.5 1.8 1.9 2.1]
        # sample = sample + 0.2*rand(4,4)
        sample = [1.0+x/8+x^2/10+y/8+y^2/12 for x=1:4,y=1:4]
        large_sample = [1.0+x/8+x^2/10+y/8+y^2/12 for x=0:5,y=0:5]
        #sample_itp = interpolate(sample, BSpline(Cubic(Line())), OnGrid())
        sample_itp = extrapolate(interpolate(large_sample, BSpline(Cubic(Line(Interpolations.OnGrid())))),Line())
        sum_err = 0.0
        for x in x_grid,y in y_grid
            my_itp = bicubicInterpolate(sample,[x,y])
            itp_result = sample_itp(3.0+x,3.0+y)
            sum_err += abs(itp_result-my_itp)/abs(itp_result)
        end
        err = sum_err/(length(x_grid)*length(y_grid))
        @info "err: "*string(err)
        @test err < 1e-2
        #test gradient
        sum_err = 0.0
        for x in x_grid,y in y_grid
            my_itp = itp_bicubic_grad(sample,[x,y],[1.0,1.0])
            itp_result = Interpolations.gradient(sample_itp,3.0+x,3.0+y)
            sum_err += norm(itp_result-my_itp)/norm(itp_result)
        end
        err = sum_err/(length(x_grid)*length(y_grid))
        @info "err: "*string(err)
        @test err < 1e-2

        @info "3D interpolation test"
        x_grid = range(0.0,stop=1.0,length=11)
        y_grid = range(0.0,stop=1.0,length=11)
        z_grid = range(0.0,stop=1.0,length=11)
        sample = [1.0+x/8+y/8+z/8+x^2/10+y^2/10+z^2/10 for x=1:4,y=1:4,z=1:4]
        large_sample = [1.0+x/8+y/8+z/8+x^2/10+y^2/10+z^2/10 for x=0:5,y=0:5,z=0:5]
        #sample_itp = interpolate(sample, BSpline(Cubic(Line())), OnGrid())
        sample_itp = extrapolate(interpolate(large_sample, BSpline(Cubic(Line(Interpolations.OnGrid())))),Line())
        sum_err = 0.0
        for x in x_grid,y in y_grid, z in z_grid
            my_itp = tricubicInterpolate(sample,[x,y,z])
            itp_result = sample_itp(3.0+x,3.0+y,3.0+z)
            sum_err += abs(itp_result-my_itp)/abs(itp_result)
        end
        err = sum_err/(length(x_grid)*length(y_grid)*length(z_grid))
        @info "err: "*string(err)
        @test err < 1e-2

        sum_err = 0.0
        for x in x_grid,y in y_grid, z in z_grid
            my_itp = itp_tricubic_grad(sample,[x,y,z],[1.0,1.0,1.0])
            itp_result = Interpolations.gradient(sample_itp,3.0+x,3.0+y,3.0+z)
            sum_err += norm(itp_result-my_itp)/norm(itp_result)
        end
        err = sum_err/(length(x_grid)*length(y_grid)*length(z_grid))
        @info "err: "*string(err)
        @test err < 1e-2
    end
end
