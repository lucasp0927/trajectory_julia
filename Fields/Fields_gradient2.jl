using Base.LinAlg.BLAS
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
    f.s::Float64 = f.scaling(t)::Float64
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field{T<:ComplexOrFloat}(f::ScalarField,pidx::Vector{Int64},s::T)
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
    f.s::Complex{Float64} = f.scaling(t)::Complex{Float64}
#    s = f.scaling(t)
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field{T<:ComplexOrFloat}(f::VectorField,pidx::Vector{Int64},s::T)
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

function sample2!(f::VectorFieldNode{2},pos::Vector{Float64},t::Real)
    fill!(f.sample,zero(Complex{Float64}))
    for ff in f.fields
        if in_field(ff,pos)
            sample2!(ff,pos,t)
            add_sample!(f.sample,ff.sample)
        end
    end
#    s::Complex{Float64} = convert(Complex{Float64},f.scaling(t))
    f.s::Complex{Float64} = f.scaling(t)::Complex{Float64}
    scal!(48,f.s,f.sample,1)
#    scale_sample!(f.sample,s)
end

@inbounds function add_sample!{T<:ComplexOrFloat}(sample1::Array{T,3},sample2::Array{T,3})
    for i in eachindex(sample1)
        sample1[i] += sample2[i]
    end
end

@inbounds function sample2!(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real)
    fill!(f.sample,zero(Float64))
    fill!(f.vf_sample,zero(Complex{Float64}))
    for ff in f.fields
        if in_field(ff,pos)
            sample2!(ff,pos,t)
            add_fields!(ff,f.sample,f.vf_sample)
        end
    end
    add_vector_field!(f.sample,f.vf_sample)
    f.s = f.scaling(t)::Float64
    scal!(16,f.s,f.sample,1)
#    scale_sample!(f.sample,s)
end

@inbounds function add_vector_field!{T<:ComplexOrFloat}(sample::Array{Float64,2},vf_sample::Array{T,3})
    for j = 1:4, i = 1:4
        #        sample[i,j] += sumabs2(vf_sample[:,i,j])
        @simd        for k = 1:3
            sample[i,j] += abs2(vf_sample[k,i,j])
        end
    end
    #    f.sample[:,:] .+= squeeze(sumabs2(f.vf_sample,1),1)[:,:]
end

# @inbounds function scale_sample!{T<:ComplexOrFloat,N,K<:ComplexOrFloat}(sample::Array{T,N},s::K)
#     for i in eachindex(sample)
#         sample[i] *=s
#     end
# end

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

@generated function value2(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real)
    quote
        x = $(Array(Float64,2))
        sample2!(f,pos,t)
        @nexprs 2 j->x[j] = rem(pos[j],f.res[j])/f.res[j]
        return bicubicInterpolate(f.sample,x)
#        return itp_bicubic(f.sample,x)
    end
    #    return itp_spline(A,(2.0+x_1,2.0+x_2))
end

@generated function value3(pos::Vector{Float64},t::Real)
    quote
        global fields
        x = $(Array(Float64,2))
        sample2!(fields::ScalarFieldNode,pos,t)
        @nexprs 2 j->x[j] = rem(pos[j],(fields::ScalarFieldNode).res[j])/(fields::ScalarFieldNode).res[j]
        return bicubicInterpolate((fields::ScalarFieldNode).sample,x)
#        return itp_bicubic(f.sample,x)
    end
    #    return itp_spline(A,(2.0+x_1,2.0+x_2))
end

