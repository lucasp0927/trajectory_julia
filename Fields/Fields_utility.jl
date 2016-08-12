# utility functions, for simple fields
# function copy_to_sharedarray!{T<:ComplexOrFloat,N}(arr::Array{T,N})
#     arr_s = SharedArray(T,size(arr))
#     arr_s[:] = arr
#     return arr_s
# end
function zero_field{T<:ComplexOrFloat,N}(::Type{ScalarField{T,N}},res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling =  t->1.0, name="scalarfield")
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    return ScalarField{T,N}(copy_to_sharedarray!(zeros(T,res...)),pos,size,scaling=scaling,name=name)
end

function zero_field{T<:ComplexOrFloat,N}(::Type{VectorField{T,N}},res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling =  t->1.0, name="vectorfield")
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    return VectorField{T,N}(copy_to_sharedarray!(zeros(T,(3,res...))),pos,size,scaling = scaling,name=name)
end
#=
function test_multiple_dispatch{T<:ComplexOrFloat}(f::ScalarField{T,2})
    println("d=2")
end

function test_multiple_dispatch{T<:ComplexOrFloat}(f::ScalarField{T,3})
    println("d=3")
end
=#
function test_zero_field()
    res_x = rand(2:100)
    res_y = rand(2:100)
    res_z = rand(2:100)
    println(res_x," ",res_y," ",res_z)
    zero_f = zero_field(ScalarField{Float64,2},[res_x,res_y],[0.0,0.0],[100.0,100.0])
    @test size(zero_f.field) == (res_x,res_y)
    @test all(map(x->x==0.0,sdata(zero_f.field)))
    zero_f = zero_field(VectorField{Float64,2},[res_x,res_y],[0.0,0.0],[100.0,100.0])
    @test size(zero_f.field) == (3,res_x,res_y)
    @test all(map(x->x==0.0,sdata(zero_f.field)))
    zero_f = zero_field(ScalarField{Float64,3},[res_x,res_y,res_z],[0.0,0.0,0.0],[100.0,100.0,100.0])
    @test size(zero_f.field) == (res_x,res_y,res_z)
    @test all(map(x->x==0.0,sdata(zero_f.field)))
    zero_f = zero_field(VectorField{Float64,3},[res_x,res_y,res_z],[0.0,0.0,0.0],[100.0,100.0,100.0])
    @test all(map(x->x==0.0,sdata(zero_f.field)))
    @test size(zero_f.field) == (3,res_x,res_y,res_z)
end

@generated function func2field{T<:ComplexOrFloat,N}(::Type{ScalarField{T,N}},func::Function,res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling =  t->1.0,name="scalarfield")
    quote
        @assert length(res)==length(pos)==length(size)==N
        @nexprs $N j->(x_j = linspace(pos[j],pos[j]+size[j],res[j]))
        f = zeros(T,res...)
        @nloops $N i j->1:length(x_j) begin
            v = func((@ntuple $N k->x_k[i_k])...)
            f[(@ntuple $N k->i_k)...] = v
        end
        return ScalarField{T,$N}(copy_to_sharedarray!(f::Array{T,$N}),pos,size,scaling=scaling,name=name)
    end
end

@generated function func2field{T<:ComplexOrFloat,N}(::Type{VectorField{T,N}},func::Function,res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling = t->1.0,name="vectorfield")
    quote
        @assert length(res)==length(pos)==length(size)==N
        @nexprs $N j->(x_j = linspace(pos[j],pos[j]+size[j],res[j]))
        f = zeros(T,(3,res...))
        @nloops $N i j->1:length(x_j) begin
            v = func((@ntuple $N k->x_k[i_k])...)
            f[:,(@ntuple $N k->i_k)...] = v
        end
        return VectorField{T,$N}(copy_to_sharedarray!(f::Array{T,$N+1}),pos,size,scaling=scaling,name=name)
    end
end
