# utility functions, for simple fields
using LinearAlgebra
using Statistics
function zero_field(::Type{ScalarField{T,N}},res::Vector{Float64},pos::Vector{Float64},size::Vector{Float64};scaling_expr::Expr = Meta.parse("t->1.0"), name::String="scalarfield") where {T <: ComplexOrFloat, N}
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    arr_size = @. floor(Int64,size/res)+1
    f_size = @. (arr_size-1)*res
	return ScalarField{T,N}(copy_to_sharedarray!(zeros(T,arr_size...)),pos,f_size,scaling_expr=scaling_expr,name=ascii(name))
end

function zero_field(::Type{VectorField{T, N}},res::Vector{Float64},pos::Vector{Float64},size::Vector{Float64};scaling_expr::Expr = Meta.parse("t->1.0"), name::String = "vectorfield") where {T <: ComplexOrFloat, N}
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    arr_size = @. floor(Int64,size/res)+1
    f_size = @. (arr_size-1)*res
    return VectorField{T,N}(copy_to_sharedarray!(zeros(T,(3,arr_size...))),pos,f_size,scaling_expr = scaling_expr,name=ascii(name))
end

@generated function func2field(::Type{ScalarField{T, N}}, func::Expr, gradx::Expr, grady::Expr, res::Vector{Float64}, pos::Vector{Float64}, size::Vector{Float64}; scaling_expr::Expr = Meta.parse("t->1.0"), name::String = "scalarfield") where {T <: ComplexOrFloat, N}
    quote
        @assert length(res)==length(pos)==length(size)==N
        arr_size = @. floor(Int64,size/res)+1
        f_size = @. (arr_size-1)*res
        for j = 1:$N
            @assert(isapprox(mean(diff(range(pos[j],stop=pos[j]+f_size[j],length=arr_size[j]))),res[j]))
        end
        @nexprs $N j->(x_j = range(pos[j],stop=pos[j]+f_size[j],length=arr_size[j]))
        f = zeros(T,arr_size...)
	# func_eval = eval(func)
        # @nloops $N i j->1:length(x_j) begin
        #     v = Base.invokelatest(func_eval,(@ntuple $N k->x_k[i_k])...)
        #     f[(@ntuple $N k->i_k)...] = v
        # end
        return ScalarFieldFunc{T,$N}(copy_to_sharedarray!(f::Array{T,$N}),pos,res,f_size,func,gradx,grady,scaling_expr=scaling_expr,name=ascii(name))
#        return ScalarField{T,$N}(copy_to_sharedarray!(f::Array{T,$N}),pos,f_size,scaling_expr=scaling_expr,name=ascii(name))        
    end
end

@generated function func2field(::Type{VectorField{T, N}}, func::Expr, res::Vector{Float64}, pos::Vector{Float64}, size::Vector{Float64}; scaling_expr::Expr = Meta.parse("t->1.0"), name::String = "vectorfield") where {T <: ComplexOrFloat, N}

    quote
        @assert length(res)==length(pos)==length(size)==N
        arr_size = @. floor(Int64,size/res)+1
        f_size = @. (arr_size-1)*res
        for j = 1:$N
            @assert(isapprox(mean(diff(range(pos[j],stop=pos[j]+f_size[j],length=arr_size[j]))),res[j]))
        end

        @nexprs $N j->(x_j = range(pos[j],stop=pos[j]+f_size[j],length=arr_size[j]))
        f = zeros(T,(3,arr_size...))
	func_eval = eval(func)
        @nloops $N i j->1:length(x_j) begin
            v = func_eval((@ntuple $N k->x_k[i_k])...)
            f[:,(@ntuple $N k->i_k)...] = v
        end
        return VectorField{T,$N}(copy_to_sharedarray!(f::Array{T,$N+1}),pos,f_size,scaling_expr=scaling_expr,name=ascii(name))
    end
end

function test_zero_field()
    res = [2.0,2.0,2.0]
    @testset "test zero_field" begin
        for dim in [2,3]
            for DT in [Float64, Complex{Float64}]
                zero_f = zero_field(ScalarField{DT,dim},res[1:dim],repeat([0.0],dim),repeat([100.0],dim))
                @test size(zero_f.field) == tuple(((repeat([100.0],dim)./res[1:dim]).+1)...)
                @test all(map(x->x==0.0,sdata(zero_f.field)))
                zero_f = zero_field(VectorField{DT,dim},res[1:dim],repeat([0.0],dim),repeat([100.0],dim))
                @test size(zero_f.field) == (3,((repeat([100.0],dim)./res[1:dim]).+1)...)
                @test all(map(x->x==0.0,sdata(zero_f.field)))
            end
        end
    end
end
