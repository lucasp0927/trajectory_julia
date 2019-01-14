using SharedArrays

function setfield!(f::VectorField{T,N},A::SharedArray{T},pos::Vector{Float64}, sz::Vector{Float64};scaling::Function = t-> begin
1
end) where {T <: ComplexOrFloat, N}
    res = sz./(collect(size(A))[2:N+1].-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A)-1 "dimension error!"
    f.field = Array{T}(undef,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

function setfield!(f::ScalarField{T,N},A::SharedArray{T},pos::Vector{Float64},sz::Array{Float64};scaling::Function = t-> begin
1
end) where {T <: ComplexOrFloat, N}
    res = sz./(collect(size(A))[1:N].-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A) "dimension error!"
    f.field = Array{T}(undef,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

# function setscaling!(f::Field,scaling::Function)
#     Base.error("scaling_expr not changed.")
#     f.scaling = scaling
# end

function setscaling!(f::Field,scaling_expr::Expr)
    f.scaling_expr = scaling_expr
    f.scaling = eval(scaling_expr)
end

function clean_scaling!(f::ScalarField{T,N}) where {T <: ComplexOrFloat, N}
    f.scaling = t->0.0
end

function clean_scaling!(f::VectorField{T,N}) where {T <: ComplexOrFloat, N}
    f.scaling = t->0.0
end

function clean_scaling!(f::ScalarFieldNode{N}) where N
    f.scaling=t->0.0
    map(clean_scaling!,f.fields)
end

function clean_scaling!(f::VectorFieldNode{N}) where N
    f.scaling=t->0.0
    map(clean_scaling!,f.fields)
end

function eval_scaling!(f::ScalarField{T,N}) where {T <: ComplexOrFloat, N}
    f.scaling = eval(f.scaling_expr)
end

function eval_scaling!(f::VectorField{T,N}) where {T <: ComplexOrFloat, N}
    f.scaling = eval(f.scaling_expr)
end

function eval_scaling!(f::ScalarFieldNode{N}) where N
    f.scaling=eval(f.scaling_expr)
    map(eval_scaling!,f.fields)
end

function eval_scaling!(f::VectorFieldNode{N}) where N
    f.scaling=eval(f.scaling_expr)
    map(eval_scaling!,f.fields)
end
