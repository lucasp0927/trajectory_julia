#1.use cubic spline interpolation
#2.use bicubic interpolation 4x4 pixels
include("Fields_interpolate.jl")
@inbounds function sample{T<:ComplexOrFloat}(f::ScalarField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
    pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64})+1)
    return (sub(f.field,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)*f.scaling(t)::ComplexOrFloat)::Array{T,2}
end

@inbounds function sample{T<:ComplexOrFloat}(f::VectorField{T,2},pos::Vector{Float64},t::Real;order::Integer = 3)
    rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
    pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64})+1)
#    return (f.field[:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2])::Array{T,3}
    return (sub(f.field,:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)*f.scaling(t)::ComplexOrFloat)::Array{T,3}
end

function sample(f::VectorFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type::DataType = f.typeof
    return sample_inner(output_type,f,pos,t,order=order)::Array{output_type,3}
end

function sample_inner{T<:ComplexOrFloat}(::Type{T},f::VectorFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    output::Array{T,3} = zeros(T,(3,4,4))
    loop_field!(f.fields,pos,t,output)
    return (output*f.scaling(t)::ComplexOrFloat)::Array{T,3}
end

function sample(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    # get (order+1)x(order+1) pixels of local field
    output_type::DataType = f.typeof
    return sample_inner(output_type,f,pos,t,order=order)::Array{Float64,2}
end

function sample_inner{T<:ComplexOrFloat}(::Type{T},f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    sfoutput::Array{Float64,2} = zeros(Float64,(4,4))#TODO: Float64 should be T
    vfoutput::Array{T,3} = zeros(T,(3,4,4))
    for ff in f.fields
        if in_field(ff,pos)
            sum_field(ff,pos,t,sfoutput,vfoutput,order = order)
        end
    end
    vfoutput_abs2::Array{Float64,2} = squeeze(sumabs2(vfoutput,1),1)
    sfoutput .+= vfoutput_abs2::Array{Float64,2}
    return (sfoutput*f.scaling(t)::ComplexOrFloat)::Array{Float64,2}
#    vf_arr::Vector{AbstractVectorField} = filter(x->isa(x,AbstractVectorField),f.fields)
#    sf_arr::Vector{AbstractScalarField} = filter(x->isa(x,AbstractScalarField),f.fields)

    # if ~isempty(sf_arr)
    #     loop_field!(sf_arr,pos,t,output)
    # end
    # if ~isempty(vf_arr)
    #     loop_field!(vf_arr,pos,t,vf_output)
    # end
end

@inbounds @fastmath function in_field{T<:Field}(f::T,pos::Vector{Float64})
    # fpos = collect(f.position)
    # fsize = collect(f.size)
    # fres = collect(f.res)
    # rel_pos = pos-fpos
    # return all(i->fres[i]<rel_pos[i]<fsize[i]-fres[i],1:length(fsize))
    return (f.res[1]::Float64<(pos[1]-f.position[1]::Float64)<(f.size[1]::Float64-f.res[1]::Float64) && f.res[2]::Float64<(pos[2]-f.position[2]::Float64)<(f.size[2]::Float64-f.res[2]::Float64))
end

function sum_field{T<:AbstractVectorField,K<:ComplexOrFloat}(f::T,pos::Vector{Float64},t::Real,sf_output::Array{Float64,2},vf_output::Array{K,3};order::Integer = 3)
    tmp = sample(f,pos,t;order=order)::Array{K}
    vf_output[:] += tmp[:]
end
function sum_field{T<:AbstractScalarField,K<:ComplexOrFloat}(f::T,pos::Vector{Float64},t::Real,sf_output::Array{Float64,2},vf_output::Array{K,3};order::Integer = 3)
    tmp = sample(f,pos,t;order=order)::Array{Float64}
    sf_output[:] += tmp[:]
end


@inbounds @fastmath function loop_field!{T<:Field,K<:ComplexOrFloat}(f_arr::Vector{T},pos::Vector{Float64},t::Real,output::Array{K};order::Integer = 3)
    ff = filter(x->in_field(x,pos),f_arr)
    for vf in ff
        tmp = sample(vf,pos,t;order=order)::Array{K}
        output[:] += tmp[:]
    end
end

function value(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
    A = sample(f,pos,t;order=order)::Array{Float64,2}
#    A = float(real(A)) #TODO: change this
    res = f.res::Vector{Float64}
    @nexprs 2 j->x_j = rem(pos[j],res[j])/res[j]
    return itp_bicubic(A,[x_1,x_2])
#    return itp_spline(A,(2.0+x_1,2.0+x_2))
end

function gradient(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real;order::Integer = 3)
#        A = Array(Float64,4,4)
        A = sample(f,pos,t;order=order)::Array{Float64,2}
        #    A = float(real(A)) #TODO: change this
        res = f.res::Vector{Float64}
        @nexprs 2 j->x_j = rem(pos[j],res[j])/res[j]
        return itp_bicubic_grad(A,[x_1,x_2],res)
    #    return itp_spline(A,(2.0+x_1,2.0+x_2))
end
