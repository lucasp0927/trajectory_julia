function setfield!{T<:ComplexOrFloat,N}(f::VectorField{T,N},A::SharedArray{T},pos::Vector{Float64},sz::Vector{Float64};scaling::Function = t->1.0)
    res = sz./(collect(size(A))[2:N+1]-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A)-1 "dimension error!"
    f.field = Array(T,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

function setfield!{T<:ComplexOrFloat,N}(f::ScalarField{T,N},A::SharedArray{T},pos::Vector{Float64},sz::Array{Float64};scaling::Function = t->1.0)
    res = sz./(collect(size(A))[1:N]-1)
    @assert all(x->x!=0,res) "zero resolution!"
    @assert length(pos)==length(sz)==N==ndims(A) "dimension error!"
    f.field = Array(T,1)
    f.field = A
    f.position = pos
    f.size = sz
    f.res = res
    f.scaling = scaling
    f.dim = N
end

function setscaling!(f::Field,scaling::Function)
    Base.error("scaling_expr not changed.")
    f.scaling = scaling
end

function setscaling!(f::Field,scaling_expr::Expr)
    f.scaling_expr = scaling_expr
    f.scaling = eval(scaling_expr)
end

function clean_scaling!{T<:ComplexOrFloat,N}(f::ScalarField{T,N})
    f.scaling = t->0.0
end

function clean_scaling!{T<:ComplexOrFloat,N}(f::VectorField{T,N})
    f.scaling = t->0.0
end

function clean_scaling!{N}(f::ScalarFieldNode{N})
    f.scaling=t->0.0
    map(clean_scaling!,f.fields)
end

function clean_scaling!{N}(f::VectorFieldNode{N})
    f.scaling=t->0.0
    map(clean_scaling!,f.fields)
end

function eval_scaling!{T<:ComplexOrFloat,N}(f::ScalarField{T,N})
    f.scaling = eval(f.scaling_expr)
end

function eval_scaling!{T<:ComplexOrFloat,N}(f::VectorField{T,N})
    f.scaling = eval(f.scaling_expr)
end

function eval_scaling!{N}(f::ScalarFieldNode{N})
    f.scaling=eval(f.scaling_expr)
    map(clean_scaling!,f.fields)
end

function eval_scaling!{N}(f::VectorFieldNode{N})
    f.scaling=eval(f.scaling_expr)
    map(clean_scaling!,f.fields)
end
