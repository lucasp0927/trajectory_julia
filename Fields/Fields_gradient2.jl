@inbounds function sample2!{T<:ComplexOrFloat}(f::ScalarField{T,2},pos::Vector{Float64},t::Real)
    rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
    pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64})+1)
    # scale = f.scaling(t)
    # for j = 1:4, i = 1:4
    #     f.sample[i:j] = f.field[pidx[1]-2+i,pidx[2]-2+j]*scale
    # end
    f.sample[:] = (sub(f.field,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)*f.scaling(t))[:]    
end

@inbounds function sample2!{T<:ComplexOrFloat}(f::VectorField{T,2},pos::Vector{Float64},t::Real)
    rel_pos::Vector{Float64} = pos-f.position::Vector{Float64}
    pidx::Vector{Int64} = round(Int64,div(rel_pos,f.res::Vector{Float64})+1)
#    scale = f.scaling(t)
    # for j = 1:4, i = 1:4,k=1:3
    #     f.sample[k,i:j] = f.field[k,pidx[1]-2+i,pidx[2]-2+j]*scale
    # end    
    f.sample[:] = (sub(f.field,:,pidx[1]-1:pidx[1]+2,pidx[2]-1:pidx[2]+2)*f.scaling(t))[:]
end

function sample2!(f::VectorFieldNode{2},pos::Vector{Float64},t::Real)
    fill!(f.sample,Base.zero(Complex{Float64}))
    for ff in f.fields
        if in_field(ff,pos)
            sample2!(ff,pos,t)
            f.sample[:] .+= ff.sample[:]
        end
    end    
    f.sample *= f.scaling(t) 
end

function sample2!(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real)
    fill!(f.sample,Base.zero(Float64))
    fill!(f.vf_sample,Base.zero(Complex{Float64}))
    for ff in f.fields
        if in_field(ff,pos)
            sample2!(ff,pos,t)
            add_fields(ff,f.sample,f.vf_sample)
        end
    end
    f.sample[:] .+= squeeze(sumabs2(f.vf_sample,1),1)[:]
end

function add_fields{T<:AbstractVectorField}(f::T,sample::Array{Float64,2},vf_sample::Array{Complex{Float64},3})
    vf_sample[:] .+= f.sample[:]
end

function add_fields{T<:AbstractScalarField}(f::T,sample::Array{Float64,2},vf_sample::Array{Complex{Float64},3})
    sample[:] .+= f.sample[:]
end

function value2(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real)
    sample2!(f,pos,t)
    A = f.sample
    res = f.res::Vector{Float64}
    @nexprs 2 j->x_j = rem(pos[j],res[j])/res[j]
    return itp_bicubic(A,[x_1,x_2])
#    return itp_spline(A,(2.0+x_1,2.0+x_2))
end

