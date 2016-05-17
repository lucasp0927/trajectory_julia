# utility functions, for simple fields
function copy_to_sharedarray!{T<:ComplexOrFloat,N}(arr::Array{T,N})
    arr_s = SharedArray(T,size(arr))
    arr_s[:] = arr
    return arr_s
end

function zero_field{T<:ComplexOrFloat,N}(::Type{ScalarField{T,N}},res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling =  t->1.0)
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    return ScalarField{T,N}(copy_to_sharedarray!(zeros(T,res...)),pos,size,scaling=scaling)
end

function zero_field{T<:ComplexOrFloat,N}(::Type{VectorField{T,N}},res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling =  t->1.0)
    @assert length(res) ==  length(pos) == length(size) "dimension mismatch"
    return VectorField{T,N}(copy_to_sharedarray!(zeros(T,(3,res...))),pos,size,scaling = scaling)
end

@generated function func2field{T<:ComplexOrFloat,N}(::Type{ScalarField{T,N}},func::Function,res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling =  t->1.0)
    quote
        @assert length(res)==length(pos)==length(size)==N
        @nexprs $N j->(x_j = linspace(pos[j],pos[j]+size[j],res[j]))
        f = zeros(T,res...)
        @nloops $N i j->1:length(x_j) begin
            v = func((@ntuple $N k->x_k[i_k])...)
            f[(@ntuple $N k->i_k)...] = v
        end
        return ScalarField{T,2}(copy_to_sharedarray!(f::Array{T,$N}),pos,size,scaling=scaling)
    end
end


# function func2field{T<:ComplexOrFloat}(::Type{ScalarField{T,2}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Function = t->1.0)
#     xx = linspace(pos[1],pos[1]+size[1],res[1])
#     yy = linspace(pos[2],pos[2]+size[2],res[2])
#     f = [func(x,y)::T for x in xx, y in yy]
#     return ScalarField{T,2}(f::Array{T,2},pos,size,scaling=scaling)
# end

# function func2field{T<:ComplexOrFloat}(::Type{ScalarField{T,3}},func::Function,res::Tuple{Integer,Integer,Integer},pos::Tuple{Real,Real,Real},size::Tuple{Real,Real,Real};scaling::Function = t->1.0)
#     xx = linspace(pos[1],pos[1]+size[1],res[1])
#     yy = linspace(pos[2],pos[2]+size[2],res[2])
#     zz = linspace(pos[3],pos[3]+size[3],res[3])
#     f = [func(x,y,z)::T for x in xx, y in yy, z in zz]
#     return ScalarField{T,3}(f::Array{T,3},pos,size,scaling=scaling)
# end


@generated function func2field{T<:ComplexOrFloat,N}(::Type{VectorField{T,N}},func::Function,res::Vector{Int64},pos::Vector{Float64},size::Vector{Float64};scaling = t->1.0)
    quote
        @assert length(res)==length(pos)==length(size)==N
        @nexprs $N j->(x_j = linspace(pos[j],pos[j]+size[j],res[j]))
        f = zeros(T,(3,res...))
        @nloops $N i j->1:length(x_j) begin
            v = func((@ntuple $N k->x_k[i_k])...)
            f[:,(@ntuple $N k->i_k)...] = v
        end
        return VectorField{T,$N}(copy_to_sharedarray!(f::Array{T,$N+1}),pos,size,scaling=scaling)
    end
end


# function func2field{T<:ComplexOrFloat}(::Type{VectorField{T,2}},func::Function,res::Tuple{Integer,Integer},pos::Tuple{Real,Real},size::Tuple{Real,Real};scaling::Function = t->1.0)
#     xx = linspace(pos[1],pos[1]+size[1],res[1])
#     yy = linspace(pos[2],pos[2]+size[2],res[2])
#     f = zeros(T,(3,res...))
#     for x in enumerate(xx), y in enumerate(yy)
#         v = func(x[2],y[2])::Vector{T}
#         f[:,x[1],y[1]] = v
#     end
#     return VectorField{T,2}(f::Array{T,3},pos,size,scaling=scaling)
# end

# function func2field{T<:ComplexOrFloat}(::Type{VectorField{T,3}},func::Function,res::Tuple{Integer,Integer,Integer},pos::Tuple{Real,Real,Real},size::Tuple{Real,Real,Real};scaling::Function = t->1.0)
#     xx = linspace(pos[1],pos[1]+size[1],res[1])
#     yy = linspace(pos[2],pos[2]+size[2],res[2])
#     zz = linspace(pos[3],pos[3]+size[3],res[3])
#     f = zeros(T,(3,res...))
#     for x in enumerate(xx), y in enumerate(yy), z in enumerate(zz)
#         v = func(x[2],y[2],z[2])::Vector{T}
#         f[:,x[1],y[1],z[1]] = v
#     end
#     return VectorField{T,3}(f::Array{T,4},pos,size,scaling=scaling)
# end
# file2field?

#expiwt(t::Float64) = exp(1.0im*t*2pi)::Complex{Float64}
