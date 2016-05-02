#1.use cubic spline interpolation
#2.use bicubic interpolation 4x4 pixels
function sample{T<:ComplexOrFLoat,N}(f::ScalarField{T,N},pos::Vector{Float64},t::Real,order::Integer = 3)
    fpos = collect(f.position)
    fres = collect(f.res)
    fsize = collect(f.size)
    pidx = div(pos/fres)+1
    
#    return T(f.field*f.scaling(t),f.position,f.size,scaling = t->1.0)
end


function sample{N}(f::ScalarFieldNode{N},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    
end
