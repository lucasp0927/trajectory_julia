# utility functions, for simple fields
using LinearAlgebra
function zero_field(::Type{ScalarField{T,N}},res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling_expr::Expr = Meta.parse("t->1.0"), name::String="scalarfield") where {T <: ComplexOrFloat, N}
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    arr_size = @. floor(Int64,size/res)
	return ScalarField{T,N}(copy_to_sharedarray!(zeros(T,arr_size...)),pos,size,scaling_expr=scaling_expr,name=ascii(name))
end

function zero_field(::Type{VectorField{T, N}},res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling_expr::Expr = Meta.parse("t->1.0"), name::String = "vectorfield") where {T <: ComplexOrFloat, N}
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    arr_size = @. floor(Int64,size/res)
    return VectorField{T,N}(copy_to_sharedarray!(zeros(T,(3,arr_size...))),pos,size,scaling_expr = scaling_expr,name=ascii(name))
end

@generated function func2field(::Type{ScalarField{T, N}}, func::Function, res::Vector{Int64}, pos::Vector{Float64}, size::Vector{Float64}; scaling_expr::Expr = Meta.parse("t->1.0"), name::String = "scalarfield") where {T <: ComplexOrFloat, N}
    quote
        @assert length(res)==length(pos)==length(size)==N
        @nexprs $N j->(x_j = linspace(pos[j],pos[j]+size[j],res[j]))
        f = zeros(T,res...)
        @nloops $N i j->1:length(x_j) begin
            v = func((@ntuple $N k->x_k[i_k])...)
            f[(@ntuple $N k->i_k)...] = v
        end
        return ScalarField{T,$N}(copy_to_sharedarray!(f::Array{T,$N}),pos,size,scaling_expr=scaling_expr,name=ascii(name))
    end
end

@generated function func2field(::Type{VectorField{T, N}}, func::Function, res::Vector{Int64}, pos::Vector{Float64}, size::Vector{Float64}; scaling_expr::Expr = Meta.parse("t->1.0"), name::String = "vectorfield") where {T <: ComplexOrFloat, N}

    quote
        @assert length(res)==length(pos)==length(size)==N
        @nexprs $N j->(x_j = linspace(pos[j],pos[j]+size[j],res[j]))
        f = zeros(T,(3,res...))
        @nloops $N i j->1:length(x_j) begin
            v = func((@ntuple $N k->x_k[i_k])...)
            f[:,(@ntuple $N k->i_k)...] = v
        end
        return VectorField{T,$N}(copy_to_sharedarray!(f::Array{T,$N+1}),pos,size,scaling_expr=scaling_expr,name=ascii(name))
    end
end

function test_zero_field()
    res = [2,2,2]
    @testset "test zero_field" begin
        for dim in [2,3]
            for DT in [Float64, Complex{Float64}]
                zero_f = zero_field(ScalarField{DT,dim},res[1:dim],repeat([0.0],dim),repeat([100.0],dim))
                @test size(zero_f.field) == tuple((repeat([100.0],dim)./res[1:dim])...)
                @test all(map(x->x==0.0,sdata(zero_f.field)))
                zero_f = zero_field(VectorField{DT,dim},res[1:dim],repeat([0.0],dim),repeat([100.0],dim))
                @test size(zero_f.field) == (3,(repeat([100.0],dim)./res[1:dim])...)
                @test all(map(x->x==0.0,sdata(zero_f.field)))
            end
        end
    end
end
