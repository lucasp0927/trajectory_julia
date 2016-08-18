#1.use cubic spline interpolation
#2.use bicubic interpolation 4x4 pixels
using Base.LinAlg.BLAS
include("Fields_interpolate.jl")
include("../constant.jl")

#@inbounds @fastmath function in_field{T<:Field}(f::T,pos::Vector{Float64})

@inbounds @fastmath function in_field{T<:Union{FieldNode2D,ComplexField2D,FloatField2D}}(f::T,pos::Vector{Float64})
    return (f.res[1]::Float64<(pos[1]-f.position[1]::Float64)<(f.size[1]::Float64-f.res[1]::Float64) && f.res[2]::Float64<(pos[2]-f.position[2]::Float64)<(f.size[2]::Float64-f.res[2]::Float64))
end

@inbounds @fastmath function in_field{T<:Union{FieldNode3D,ComplexField3D,FloatField3D}}(f::T,pos::Vector{Float64})
    # periodic in z direction?
    return (f.res[1]::Float64<(pos[1]-f.position[1]::Float64)<(f.size[1]::Float64-f.res[1]::Float64) && f.res[2]::Float64<(pos[2]-f.position[2]::Float64)<(f.size[2]::Float64-f.res[2]::Float64) && f.res[3]::Float64<(pos[3]-f.position[3]::Float64)<(f.size[3]::Float64-f.res[3]::Float64))
end


#2D Scalar sample
@fastmath @inbounds function sample2!{T<:ComplexOrFloat}(f::ScalarField{T,2},pos::Vector{Float64},t::Real)
    f.rel_pos[1] = pos[1]-f.position[1]::Float64
    f.rel_pos[2] = pos[2]-f.position[2]::Float64
    f.pidx[1] = round(Int64,div(f.rel_pos[1],f.res[1]::Float64))
    f.pidx[2] = f.pidx[1]+3
    f.pidx[3] = round(Int64,div(f.rel_pos[2],f.res[2]::Float64))
    f.pidx[4] = f.pidx[3]+3
#    rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
#    pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64}))
    # pidx = Array(Int64,2)
    # pidx[1] = round(Int64,cld(rel_pos[1],f.res[1]::Float64))
    # pidx[2] = round(Int64,cld(rel_pos[2],f.res[2]::Float64))
    f.s::Float64 = (f.scaling::Function)(t)::Float64
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field{T<:ComplexOrFloat,K<:ComplexOrFloat}(f::ScalarField{T,2},pidx::Vector{Int64},s::K)
#    blascopy!(16,f.field[pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2],1,f.sample,1)
    copy!(f.sample,sub(f.field,pidx[1]:pidx[2],pidx[3]:pidx[4]))
    scal!(16,s,f.sample,1)
    # for j=1:4,i=1:4
    #     f.sample[i,j] = f.field[pidx[1]-2+i,pidx[2]-2+j]
    # end
    # for i in eachindex(f.sample)
    #     f.sample[i] *=s
    # end
end
#3D scalar sample
@fastmath @inbounds function sample2!{T<:ComplexOrFloat}(f::ScalarField{T,3},pos::Vector{Float64},t::Real)
    f.rel_pos[1] = pos[1]-f.position[1]::Float64
    f.rel_pos[2] = pos[2]-f.position[2]::Float64
    f.rel_pos[3] = pos[3]-f.position[3]::Float64
    f.pidx[1] = round(Int64,div(f.rel_pos[1],f.res[1]::Float64))
    f.pidx[2] = f.pidx[1]+3
    f.pidx[3] = round(Int64,div(f.rel_pos[2],f.res[2]::Float64))
    f.pidx[4] = f.pidx[3]+3
    f.pidx[5] = round(Int64,div(f.rel_pos[3],f.res[3]::Float64))
    f.pidx[6] = f.pidx[5]+3
    f.s::Float64 = (f.scaling::Function)(t)::Float64
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field{T<:ComplexOrFloat,K<:ComplexOrFloat}(f::ScalarField{T,3},pidx::Vector{Int64},s::K)
    copy!(f.sample,sub(f.field,pidx[1]:pidx[2],pidx[3]:pidx[4]),pidx[5]:pidx[6])
    scal!(64,s,f.sample,1)
end
#2D vector sample
@fastmath @inbounds function sample2!{T<:ComplexOrFloat}(f::VectorField{T,2},pos::Vector{Float64},t::Real)
    f.rel_pos[1] = pos[1]-f.position[1]::Float64
    f.rel_pos[2] = pos[2]-f.position[2]::Float64
    f.pidx[1] = round(Int64,div(f.rel_pos[1],f.res[1]::Float64))
    f.pidx[2] = f.pidx[1]+3
    f.pidx[3] = round(Int64,div(f.rel_pos[2],f.res[2]::Float64))
    f.pidx[4] = f.pidx[3]+3
    # rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
    # pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64}))
    # pidx = Array(Int64,2)
    # pidx[1] = round(Int64,cld(rel_pos[1],f.res[1]::Float64))
    # pidx[2] = round(Int64,cld(rel_pos[2],f.res[2]::Float64))
#    s::Complex{Float64} = convert(Complex{Float64},f.scaling(t))
    f.s::Complex{Float64} = (f.scaling::Function)(t)::Complex{Float64}
#    s = f.scaling(t)
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field{T<:ComplexOrFloat,K<:ComplexOrFloat}(f::VectorField{T,2},pidx::Vector{Int64},s::K)
#    blascopy!(48,f.field[:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2],1,f.sample,1)
    #    f.sample[:,:,:] = sub(f.field,:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)
    copy!(f.sample,sub(f.field,:,pidx[1]:pidx[2],pidx[3]:pidx[4]))
    # for j=1:4,i=1:4,k=1:3
    #     f.sample[k,i,j] = f.field[k,pidx[1]-2+i,pidx[2]-2+j]
    # end
    scal!(48,s,f.sample,1)
    # for i in eachindex(f.sample)
    #     f.sample[i] *=s
    # end
end
#3D vector sample
@fastmath @inbounds function sample2!{T<:ComplexOrFloat}(f::VectorField{T,3},pos::Vector{Float64},t::Real)
    f.rel_pos[1] = pos[1]-f.position[1]::Float64
    f.rel_pos[2] = pos[2]-f.position[2]::Float64
    f.rel_pos[3] = pos[3]-f.position[3]::Float64
    f.pidx[1] = round(Int64,div(f.rel_pos[1],f.res[1]::Float64))
    f.pidx[2] = f.pidx[1]+3
    f.pidx[3] = round(Int64,div(f.rel_pos[2],f.res[2]::Float64))
    f.pidx[4] = f.pidx[3]+3
    f.pidx[5] = round(Int64,div(f.rel_pos[3],f.res[3]::Float64))
    f.pidx[6] = f.pidx[5]+3
    f.s::Complex{Float64} = (f.scaling::Function)(t)::Complex{Float64}
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field{T<:ComplexOrFloat,K<:ComplexOrFloat}(f::VectorField{T,3},pidx::Vector{Int64},s::K)
    copy!(f.sample,sub(f.field,:,pidx[1]:pidx[2],pidx[3]:pidx[4],pidx[5]:pidx[6]))
    scal!(192,s,f.sample,1)
end
#2D Vector Field Node
function sample2!(f::VectorFieldNode{2},pos::Vector{Float64},t::Real)
    fill!(f.sample,zero(Complex{Float64}))
    for ff in f.fields
        if in_field(ff,pos)
            sample2!(ff,pos,t)
            add_sample!(f.sample,ff.sample)
        end
    end
#    s::Complex{Float64} = convert(Complex{Float64},f.scaling(t))
    f.s::Complex{Float64} = (f.scaling::Function)(t)::Complex{Float64}
    scal!(48,f.s,f.sample,1)
#    scale_sample!(f.sample,s)
end

@fastmath @inbounds function add_sample!{T<:ComplexOrFloat}(sample1::Array{T,3},sample2::Array{T,3})
    for i in eachindex(sample1)
        sample1[i] += sample2[i]
    end
end
#3D Vector Field Node
function sample2!(f::VectorFieldNode{3},pos::Vector{Float64},t::Real)
    fill!(f.sample,zero(Complex{Float64}))
    for ff in f.fields
        if in_field(ff,pos)
            sample2!(ff,pos,t)
            add_sample!(f.sample,ff.sample)
        end
    end
    f.s::Complex{Float64} = (f.scaling::Function)(t)::Complex{Float64}
    scal!(192,f.s,f.sample,1)
end

@fastmath @inbounds function add_sample!{T<:ComplexOrFloat}(sample1::Array{T,4},sample2::Array{T,4})
    for i in eachindex(sample1)
        sample1[i] += sample2[i]
    end
end
#2D Scalar Field Node
@inbounds function sample2!(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real)
    fill!((f.sample)::Array{Float64,2},zero(Float64))
     if f.one_vf_flag
         sample2!(f.fields[1],pos,t)
         add_vector_field!(f.sample,f.fields[1].sample)
     else
        fill!((f.vf_sample)::Array{Complex{Float64},3},zero(Complex{Float64}))
        for ff in f.fields
            if in_field(ff,pos)
                sample2!(ff,pos,t)
                add_fields!(ff,f.sample,f.vf_sample)
            end
        end
        add_vector_field!(f.sample,f.vf_sample)
     end
        f.s = (f.scaling::Function)(t)::Float64
        scal!(16,f.s,f.sample,1)
#    scale_sample!(f.sample,s)
end

@fastmath @inbounds function add_vector_field!{T<:ComplexOrFloat}(sample::Array{Float64,2},vf_sample::Array{T,3})
    for j = 1:4, i = 1:4
        #        sample[i,j] += sumabs2(vf_sample[:,i,j])
        @simd        for k = 1:3
            sample[i,j] += abs2(vf_sample[k,i,j])
        end
    end
    #    f.sample[:,:] .+= squeeze(sumabs2(f.vf_sample,1),1)[:,:]
end

@fastmath @inbounds function add_fields!{T<:AbstractVectorField}(f::T,sample::Array{Float64,2},vf_sample::Array{Complex{Float64},3})
    for i in eachindex(vf_sample)
        vf_sample[i] += f.sample[i]
    end
#    vf_sample[:] .+= f.sample[:]
end

@fastmath @inbounds function add_fields!{T<:AbstractScalarField}(f::T,sample::Array{Float64,2},vf_sample::Array{Complex{Float64},3})
    for i in eachindex(sample)
        sample[i] += f.sample[i]
    end
#    sample[:] .+= f.sample[:]
end
#3D Scalar Field Node
@inbounds function sample2!(f::ScalarFieldNode{3},pos::Vector{Float64},t::Real)
    fill!((f.sample)::Array{Float64,3},zero(Float64))
    if f.one_vf_flag
        sample2!(f.fields[1],pos,t)
        add_vector_field!(f.sample,f.fields[1].sample)
    else
        fill!((f.vf_sample)::Array{Complex{Float64},4},zero(Complex{Float64}))
        for ff in f.fields
            if in_field(ff,pos)
                sample2!(ff,pos,t)
                add_fields!(ff,f.sample,f.vf_sample)
            end
        end
        add_vector_field!(f.sample,f.vf_sample)
    end
    f.s = (f.scaling::Function)(t)::Float64
    scal!(64,f.s,f.sample,1)
    #    scale_sample!(f.sample,s)
end

@fastmath @inbounds function add_vector_field!{T<:ComplexOrFloat}(sample::Array{Float64,3},vf_sample::Array{T,4})
    for k=1:4, j = 1:4, i = 1:4
        #        sample[i,j] += sumabs2(vf_sample[:,i,j])
        @simd        for s = 1:3
            sample[i,j,k] += abs2(vf_sample[s,i,j,k])
        end
    end
    #    f.sample[:,:] .+= squeeze(sumabs2(f.vf_sample,1),1)[:,:]
end

@fastmath @inbounds function add_fields!{T<:AbstractVectorField}(f::T,sample::Array{Float64,3},vf_sample::Array{Complex{Float64},4})
    for i in eachindex(vf_sample)
        vf_sample[i] += f.sample[i]
    end
#    vf_sample[:] .+= f.sample[:]
end

@fastmath @inbounds function add_fields!{T<:AbstractScalarField}(f::T,sample::Array{Float64,3},vf_sample::Array{Complex{Float64},4})
    for i in eachindex(sample)
        sample[i] += f.sample[i]
    end
#    sample[:] .+= f.sample[:]
end
#########################

#=
@generated function value2D(pos::Vector{Float64},t::Real)
    quote
        x = $(Array(Float64,2))
        res = $(Array(Float64,2))
        res[:] = (fields::ScalarFieldNode).res
        sample2!(fields::ScalarFieldNode,pos,t)
        @nexprs 2 j->x[j] = rem(pos[j],res[j])/res[j]
        return bicubicInterpolate((fields::ScalarFieldNode).sample,x)
    end
end

@generated function value3D(pos::Vector{Float64},t::Real)
    quote
        x = $(Array(Float64,3))
        res = $(Array(Float64,3))
        res[:] = (fields::ScalarFieldNode).res
        sample2!(fields::ScalarFieldNode,pos,t)
        @nexprs 3 j->x[j] = rem(pos[j],res[j])/res[j]
        return tricubicInterpolate((fields::ScalarFieldNode).sample,x)
    end
end
=#

function value(pos::Vector{Float64},t::Real)
    value(pos,t,fields)
end

@generated function value(pos::Vector{Float64},t::Real,sfn::ScalarFieldNode{2})
    quote
        x = $(Array(Float64,2))
        res = $(Array(Float64,2))
        res[:] = (sfn::ScalarFieldNode).res
        sample2!(sfn::ScalarFieldNode,pos,t)
        @nexprs 2 j->x[j] = rem(pos[j],res[j])/res[j]
        return bicubicInterpolate((sfn::ScalarFieldNode).sample,x)
    end
end

@generated function value(pos::Vector{Float64},t::Real,sfn::ScalarFieldNode{3})
    quote
        x = $(Array(Float64,3))
        res = $(Array(Float64,3))
        res[:] = (sfn::ScalarFieldNode).res
        sample2!(sfn::ScalarFieldNode,pos,t)
        @nexprs 3 j->x[j] = rem(pos[j],res[j])/res[j]
        return tricubicInterpolate((sfn::ScalarFieldNode).sample,x)
    end
end

function gradient!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64})
    gradient!(t,posvel,grad,fields)
end

@generated function gradient!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64},sfn::ScalarFieldNode{2})
    quote
        x = $(Array(Float64,2))
        res = $(Array(Float64,2))
        pos = $(Array(Float64,2))
        res[:] = (sfn::ScalarFieldNode).res
        pos[:] = posvel[1:2]
        @nexprs 2 j->x[j] = rem(pos[j],res[j])/res[j]
        sample2!(sfn::ScalarFieldNode,pos,t)
        grad[1] = posvel[3]
        grad[2] = posvel[4]
        grad[3:4] = -1.0*itp_bicubic_grad((sfn::ScalarFieldNode).sample,x,res)*KB/M_CS
    end
end

@generated function gradient!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64},sfn::ScalarFieldNode{3})
    quote
        x = $(Array(Float64,3))
        res = $(Array(Float64,3))
        pos = $(Array(Float64,3))
        res[:] = (sfn::ScalarFieldNode).res
        pos[:] = posvel[1:3]
        @nexprs 3 j->x[j] = rem(pos[j],res[j])/res[j]
        sample2!(sfn::ScalarFieldNode,pos,t)
        grad[1] = posvel[4]
        grad[2] = posvel[5]
        grad[3] = posvel[6]
        grad[4:6] = -1.0*itp_bicubic_grad((sfn::ScalarFieldNode).sample,x,res)*KB/M_CS
    end
end

function test_gradient()

end
