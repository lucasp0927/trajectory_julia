#1.use cubic spline interpolation
#2.use bicubic interpolation 4x4 pixels
using LinearAlgebra
include("Fields_interpolate.jl")
include("../constant.jl")

#@inbounds @fastmath function in_field{T<:Field}(f::T,pos::Vector{Float64})

@inbounds @fastmath function in_field(f::T, pos::Vector{Float64}) where T <: Union{FieldNode2D, ComplexField2D, FloatField2D}
    return (
        3*f.res[1]::Float64<(pos[1]-f.position[1]::Float64)<(f.size[1]::Float64-3*f.res[1]::Float64) &&
        3*f.res[2]::Float64<(pos[2]-f.position[2]::Float64)<(f.size[2]::Float64-3*f.res[2]::Float64)
    )
end

@inbounds @fastmath function in_field(f::T, pos::Vector{Float64}) where T <: Union{FieldNode3D, ComplexField3D, FloatField3D}
    # periodic in z direction?
    return (
        3*f.res[1]::Float64<(pos[1]-f.position[1]::Float64)<(f.size[1]::Float64-3*f.res[1]::Float64) &&
        3*f.res[2]::Float64<(pos[2]-f.position[2]::Float64)<(f.size[2]::Float64-3*f.res[2]::Float64) &&
        3*f.res[3]::Float64<(pos[3]-f.position[3]::Float64)<(f.size[3]::Float64-3*f.res[3]::Float64)
    )
end

#2D Scalar sample
@fastmath @inbounds function sample2!(f::ScalarField{T, 2}, pos::Vector{Float64}, t::Real) where T <: ComplexOrFloat
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

    #f.s::Float64 = (f.scaling::Function)(t)::Float64
    f.s::Float64 = Base.invokelatest(f.scaling, t)::Float64
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field(f::ScalarField{T, 2}, pidx::Vector{Int64}, s::K) where {T <: ComplexOrFloat, K <: ComplexOrFloat}
    #f.sample[:,:] = view(f.field,pidx[1]:pidx[2],pidx[3]:pidx[4]).*s
    copyto!(f.sample::Array{T,2},view(f.field,pidx[1]:pidx[2],pidx[3]:pidx[4]))
    LinearAlgebra.BLAS.scal!(16,s,f.sample,1)
    #scale!(f.sample::Array{T,2},s)
#    scal!(16,s,f.sample,1)
end

#3D scalar sample
@fastmath @inbounds function sample2!(f::ScalarField{T, 3}, pos::Vector{Float64}, t::Real) where T <: ComplexOrFloat
    f.rel_pos[1] = pos[1]-f.position[1]::Float64
    f.rel_pos[2] = pos[2]-f.position[2]::Float64
    f.rel_pos[3] = pos[3]-f.position[3]::Float64
    f.pidx[1] = round(Int64,div(f.rel_pos[1],f.res[1]::Float64))
    f.pidx[2] = f.pidx[1]+3
    f.pidx[3] = round(Int64,div(f.rel_pos[2],f.res[2]::Float64))
    f.pidx[4] = f.pidx[3]+3
    f.pidx[5] = round(Int64,div(f.rel_pos[3],f.res[3]::Float64))
    f.pidx[6] = f.pidx[5]+3
    f.s::Float64 = Base.invokelatest(f.scaling::Function,t)::Float64
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field(f::ScalarField{T, 3}, pidx::Vector{Int64}, s::K) where {T <: ComplexOrFloat, K <: ComplexOrFloat}
    copyto!(f.sample::Array{T,3},view(f.field,pidx[1]:pidx[2],pidx[3]:pidx[4],pidx[5]:pidx[6]))
    LinearAlgebra.BLAS.scal!(64,s,f.sample,1)
    # copy!(f.sample,view(f.field,pidx[1]:pidx[2],pidx[3]:pidx[4],pidx[5]:pidx[6]))
    # scal!(64,s,f.sample,1)
end
#2D vector sample
@fastmath @inbounds function sample2!(f::VectorField{T, 2}, pos::Vector{Float64}, t::Real) where T <: ComplexOrFloat
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
    #f.s::Complex{Float64} = convert(Complex{Float64},(f.scaling::Function)(t))
    f.s = convert(Complex{Float64},Base.invokelatest(f.scaling::Function,t))::Complex{Float64}
#    s = f.scaling(t)
    sample_field(f,f.pidx,f.s)
end

#function sample_field(f::VectorField{T, 2}, pidx::Vector{Int64}, s::K) where {T <: ComplexOrFloat, K <: Float64}
#    f.sample[1:3,1:4,1:4] = f.field[:,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)]*s;
    #@views f.sample[1:3,1:4,1:4] .= @views f.field[:,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)].*s;
    #    @inbounds copyto!(f.sample,view(f.field,:,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)))
    #    @fastmath f.sample .*= s
#    @. @views f.sample[1:3,1:4,1:4] = (f.field::SharedArray{T,3})[:,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)]*s
#    f.sample[1:3,1:4,1:4] .= @views f.field[:,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)].*s
    #copy!(f.sample::Array{T,3},@view f.field[:,pidx[1]:pidx[2],pidx[3]:pidx[4]])
    #scale!(f.sample::Array{T,3},s)
#end

function sample_field(f::VectorField{T, 2}, pidx::Vector{Int64}, s::K) where {T <: Float64, K <: Complex{Float64}}
    @. @views f.sample[1:3,1:4,1:4] = (f.field::SharedArray{T,3})[:,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)]*s
end

function sample_field(f::VectorField{T, 2}, pidx::Vector{Int64}, s::K) where {T <: Complex{Float64}, K <: Complex{Float64}}
    #@. @views f.sample[1:3,1:4,1:4] = (f.field::SharedArray{T,3})[1:3,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)]*s
    #use more memory but slightly faster!
    @views f.sample = (f.field::SharedArray{T,3})[1:3,pidx[1]:(pidx[1]+3),pidx[3]:(pidx[3]+3)]*s
end

#function copy_scale!{T<:ComplexOrFloat}(dest::AbstractArray,src::AbstractArray,s::T)
#    @simd for i in eachindex(dest)
#        @inbounds dest[i] = src[i]*s
#    end
#end

#3D vector sample
@fastmath @inbounds function sample2!(f::VectorField{T, 3}, pos::Vector{Float64}, t::Real) where T <: ComplexOrFloat
    f.rel_pos[1] = pos[1]-f.position[1]::Float64
    f.rel_pos[2] = pos[2]-f.position[2]::Float64
    f.rel_pos[3] = pos[3]-f.position[3]::Float64
    f.pidx[1] = round(Int64,div(f.rel_pos[1],f.res[1]::Float64))
    f.pidx[2] = f.pidx[1]+3
    f.pidx[3] = round(Int64,div(f.rel_pos[2],f.res[2]::Float64))
    f.pidx[4] = f.pidx[3]+3
    f.pidx[5] = round(Int64,div(f.rel_pos[3],f.res[3]::Float64))
    f.pidx[6] = f.pidx[5]+3
    f.s = convert(Complex{Float64},Base.invokelatest(f.scaling::Function,t))::Complex{Float64}
    #f.s::Complex{Float64} = convert(Complex{Float64},(f.scaling::Function)(t))
    sample_field(f,f.pidx,f.s)
end

@inbounds function sample_field(f::VectorField{T, 3}, pidx::Vector{Int64}, s::K) where {T <: ComplexOrFloat, K <: ComplexOrFloat}
    #copy!(f.sample,view(f.field,:,pidx[1]:pidx[2],pidx[3]:pidx[4],pidx[5]:pidx[6]))
    #scal!(192,s,f.sample,1)
    @views f.sample = (f.field::SharedArray{T,4})[1:3,pidx[1]:pidx[2],pidx[3]:pidx[4],pidx[5]:pidx[6]]*s
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
    #    f.s::Complex{Float64} = convert(Complex{Float64},f.scaling(t))
    #f.s::Complex{Float64} = ((f.scaling::Function)(t))::Complex{Float64}
    f.s::Complex{Float64} = (Base.invokelatest(f.scaling::Function,t))
    f.sample .*= f.s
    #scale!(f.sample,f.s)
#    scal!(48,f.s,f.sample,1)
#    scale_sample!(f.sample,s)
end

@fastmath @inbounds function add_sample!(sample1::Array{T, 3}, sample2::Array{T, 3}) where T <: ComplexOrFloat
    sample1 .+= sample2
    # @simd for i in eachindex(sample1)
    #     sample1[i] += sample2[i]
    # end
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
    f.s = convert(Complex{Float64},Base.invokelatest(f.scaling::Function,t))::Complex{Float64}
    #f.s = (f.scaling::Function)(t)::Complex{Float64}
    f.sample .*= f.s
    #scal!(192,f.s,f.sample,1)
end

@fastmath @inbounds function add_sample!(sample1::Array{T, 4}, sample2::Array{T, 4}) where T <: ComplexOrFloat
    sample1 .+= sample2
    # for i in eachindex(sample1)
    #     sample1[i] += sample2[i]
    # end
end
#2D Scalar Field Node
@inbounds function sample2!(f::ScalarFieldNode{2},pos::Vector{Float64},t::Real)
    fill!((f.sample)::Array{Float64,2},zero(Float64))
    if f.one_vf_flag
        if in_field(f.fields[1],pos)
            sample2!(f.fields[1],pos,t)
            add_vector_field!(f.sample,f.fields[1].sample)
        end
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
    #f.s = (f.scaling::Function)(t)::Float64
    f.s = Base.invokelatest(f.scaling, t)::Float64
    f.sample .*= f.s
    #scale!(f.sample,f.s)
#        scal!(16,f.s,f.sample,1)
#    scale_sample!(f.sample,s)
end

@fastmath @inbounds function add_vector_field!(sample::Array{Float64, 2}, vf_sample::Array{T, 3}) where T <: ComplexOrFloat
    for j = 1:4, i = 1:4
        @simd for k = 1:3
            sample[i,j] += abs2(vf_sample[k,i,j])
        end
    end
end

@fastmath @inbounds function add_fields!(f::T, sample::Array{Float64, 2}, vf_sample::Array{Complex{Float64}, 3}) where T <: AbstractVectorField
    # @simd for i in eachindex(vf_sample)
    #     vf_sample[i] += f.sample[i]
    # end
    vf_sample .+= f.sample
end

@fastmath @inbounds function add_fields!(f::T, sample::Array{Float64, 2}, vf_sample::Array{Complex{Float64}, 3}) where T <: AbstractScalarField
    # @simd for i in eachindex(sample)
    #     sample[i] += f.sample[i]
    # end
    sample .+= f.sample
end
#3D Scalar Field Node
@inbounds function sample2!(f::ScalarFieldNode{3},pos::Vector{Float64},t::Real)
    fill!((f.sample)::Array{Float64,3},zero(Float64))
    if f.one_vf_flag
        if in_field(f.fields[1],pos)
            sample2!(f.fields[1],pos,t)
            add_vector_field!(f.sample,f.fields[1].sample)
        end
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
    f.s = Base.invokelatest(f.scaling::Function,t)::Float64
    #scal!(64,f.s,f.sample,1)
    f.sample .*= f.s
    #    scale_sample!(f.sample,s)
end

@fastmath @inbounds function add_vector_field!(sample::Array{Float64, 3}, vf_sample::Array{T, 4}) where T <: ComplexOrFloat
    # for j=1:4, i=1:4
    #     sample[i,j] += sumabs2(vf_sample[:,i,j])
    # end
    for k=1:4, j = 1:4, i = 1:4
        @simd        for s = 1:3
            sample[i,j,k] += abs2(vf_sample[s,i,j,k])
        end
    end
end

@fastmath @inbounds function add_fields!(f::T, sample::Array{Float64, 3}, vf_sample::Array{Complex{Float64}, 4}) where T <: AbstractVectorField
    @simd for i in eachindex(vf_sample)
        vf_sample[i] += f.sample[i]
    end
#    vf_sample[:] .+= f.sample[:]
end

@fastmath @inbounds function add_fields!(f::T, sample::Array{Float64, 3}, vf_sample::Array{Complex{Float64}, 4}) where T <: AbstractScalarField
    @simd for i in eachindex(sample)
        sample[i] += f.sample[i]
    end
#    sample[:] .+= f.sample[:]
end

function value(pos::Vector{Float64},t::Real)
    value(pos,t,fields_arr)
end

function value(pos::Vector{Float64},t::Real,f_arr::Vector{ScalarFieldNode{N}}) where N
    mapreduce(x->value(pos,t,x),+,f_arr,init=0.0)
end

@generated function value(pos::Vector{Float64},t::Real,sfn::ScalarFieldNode{2})
    quote
        if in_field(sfn,pos)
            x = $(Array{Float64}(undef,2))
            res = $(Array{Float64}(undef,2))
            res .= (sfn::ScalarFieldNode).res
            sample2!(sfn::ScalarFieldNode,pos,t)
            @nexprs 2 j->x[j] = rem((pos[j]-sfn.position[j]),res[j])/res[j]
            return bicubicInterpolate((sfn::ScalarFieldNode).sample,x)::Float64
        else
            return zero(Float64)
        end
    end
end

@generated function value(pos::Vector{Float64},t::Real,sfn::ScalarFieldNode{3})
    quote
        if in_field(sfn,pos)
            x = $(Array{Float64}(undef,3))
            res = $(Array{Float64}(undef,3))
            res[:] = (sfn::ScalarFieldNode).res
            sample2!(sfn::ScalarFieldNode,pos,t)
            @nexprs 3 j->x[j] = rem(pos[j]-sfn.position[j],res[j])/res[j]
            return tricubicInterpolate((sfn::ScalarFieldNode).sample,x)::Float64
        else
            return zero(Float64)
        end
    end
end

function gradient!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64})
    #gradient!(t,posvel,grad,fields_arr::Vector{ScalarFieldNode})
    gradient!(t,posvel,grad,fields_arr)
end

function gradient_odejl!(grad::Vector{Float64},posvel::Vector{Float64},p,t::Float64,)
    #gradient!(t,posvel,grad,fields_arr::Vector{ScalarFieldNode})
    gradient!(t,posvel,grad,fields_arr)
end

function gradient!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64},f_arr::Vector{ScalarFieldNode{N}}) where N
    fill!(grad,zero(Float64))
    for f in f_arr
        gradient_add!(t,posvel,grad,f)
    end
end

@generated function gradient_add!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64},sfn::ScalarFieldNode{2})
    quote
        x = $(Array{Float64}(undef,2))
        res = $(Array{Float64}(undef,2))
        pos = $(Array{Float64}(undef,2))
        res[:] = (sfn::ScalarFieldNode).res
        pos[:] = @view posvel[1:2]
        @nexprs 2 j->x[j] = rem((pos[j]-sfn.position[j]),res[j])/res[j]
        sample2!(sfn::ScalarFieldNode,pos,t)
        grad[1] = posvel[3]
        grad[2] = posvel[4]
        grad[3:4] += -1.0*itp_bicubic_grad((sfn::ScalarFieldNode).sample,x,res)*KB/M_CS
    end
end

@generated function gradient_add!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64},sfn::ScalarFieldNode{3})
    quote
        x = $(Array{Float64}(undef,3))
        res = $(Array{Float64}(undef,3))
        pos = $(Array{Float64}(undef,3))
        res[:] = (sfn::ScalarFieldNode).res
        pos[:] = posvel[1:3]
        @nexprs 3 j->x[j] = rem(pos[j]-sfn.position[j],res[j])/res[j]
        sample2!(sfn::ScalarFieldNode,pos,t)
        grad[1] = posvel[4]
        grad[2] = posvel[5]
        grad[3] = posvel[6]
        grad[4:6] += -1.0*itp_tricubic_grad((sfn::ScalarFieldNode).sample,x,res)*KB/M_CS
    end
end


@generated function gradient!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64},sfn::ScalarFieldNode{2})
    quote
        x = $(Array{Float64}(undef,2))
        res = $(Array{Float64}(undef,2))
        pos = $(Array{Float64}(undef,2))
        res[:] = (sfn::ScalarFieldNode).res
        pos[:] = @view posvel[1:2]
        @nexprs 2 j->x[j] = rem((pos[j]-sfn.position[j]),res[j])/res[j]
        sample2!(sfn::ScalarFieldNode,pos,t)
        grad[1] = posvel[3]
        grad[2] = posvel[4]
        grad[3:4] = -1.0*itp_bicubic_grad((sfn::ScalarFieldNode).sample,x,res)*KB/M_CS
    end
end

@generated function gradient!(t::Float64,posvel::Vector{Float64},grad::Vector{Float64},sfn::ScalarFieldNode{3})
    quote
        x = $(Array{Float64}(undef,3))
        res = $(Array{Float64}(undef,3))
        pos = $(Array{Float64}(undef,3))
        res[:] = (sfn::ScalarFieldNode).res
        pos[:] = posvel[1:3]
        @nexprs 3 j->x[j] = rem(pos[j]-sfn.position[j],res[j])/res[j]
        sample2!(sfn::ScalarFieldNode,pos,t)
        grad[1] = posvel[4]
        grad[2] = posvel[5]
        grad[3] = posvel[6]
        grad[4:6] = -1.0*itp_tricubic_grad((sfn::ScalarFieldNode).sample,x,res)*KB/M_CS
    end
end


################################################
#Test
################################################
function sa_init(S::SharedArray)
    for id in SharedArrays.localindices(S)
        S[id] = func4_3d2(Tuple(CartesianIndices(S)[id])...)
    end
end
func1_2d(x,y)=(sin(x/40.0)+cos(y/20.0))*exp(-((x-200.0)^2+(y-200.0)^2)/100^2)::Float64
func2_2d(x,y)=(sin(x/20.0)+cos(y/50.0))*exp(-((x-400.0)^2+(y-400.0)^2)/200^2)::Float64
func3_2d(x,y)=(sin(x/50.0)+cos(y/30.0))*exp(-((x-600.0)^2+(y-600.0)^2)/300^2)::Float64
func4_2d(x,y)=[(cos(x/20.0)+sin(y/35.0))*exp(-((x-500.0)^2+(y-500.0)^2)/400^2),(cos(x/20.0)+sin(y/35.0))*exp(-((x-500.0)^2+(y-500.0)^2)/400^2)*1im,(cos(x/19.0)+sin(y/25.0))*exp(-((x-500.0)^2+(y-500.0)^2)/300^2)*(1+1im)/sqrt(2.0)]::Vector{Complex{Float64}}
func1_3d(x,y,z)=(sin(x/40.0)+cos(y/20.0)+sin(z/30.0))*exp(-((x-200.0)^2+(y-200.0)^2)/100^2)*exp(-(z-50.0)^2/30^2)::Float64
func2_3d(x,y,z)=(sin(x/20.0)+cos(y/50.0)+sin(z/30.0))*exp(-((x-400.0)^2+(y-400.0)^2)/200^2)*exp(-(z-50.0)^2/30^2)::Float64
func3_3d(x,y,z)=(sin(x/50.0)+cos(y/30.0)+sin(z/30.0))*exp(-((x-600.0)^2+(y-600.0)^2)/300^2)*exp(-(z-50.0)^2/30^2)::Float64
func4_3d(x,y,z)=[(cos(x/20.0)+sin(y/35.0)+sin(z/30.0))*exp(-((x-500.0)^2+(y-500.0)^2)/400^2)*exp(-(z-50.0)^2/30^2),(cos(x/20.0)+sin(y/35.0)+sin(z/30.0))*exp(-((x-500.0)^2+(y-500.0)^2)/400^2)*exp(-(z-50.0)^2/30^2)*1im,(cos(x/19.0)+sin(y/25.0)+sin(z/30.0))*exp(-((x-500.0)^2+(y-500.0)^2)/300^2)*exp(-(z-50.0)^2/30^2)*(1+1im)/sqrt(2.0)]::Vector{Complex{Float64}}
function func4_3d2(i,x,y,z)
    if i ==1
        return (cos(x/20.0)+sin(y/35.0)+sin(z/30.0))*exp(-((x-500.0)^2+(y-500.0)^2)/400^2)*exp(-(z-50.0)^2/30^2)
    elseif i ==2
        return (cos(x/20.0)+sin(y/35.0)+sin(z/30.0))*exp(-((x-500.0)^2+(y-500.0)^2)/400^2)*exp(-(z-50.0)^2/30^2)*1im
    else
        return (cos(x/19.0)+sin(y/25.0)+sin(z/30.0))*exp(-((x-500.0)^2+(y-500.0)^2)/300^2)*exp(-(z-50.0)^2/30^2)*(1+1im)/sqrt(2.0)
    end
end


###################
# TEST
###################
function test_gradient_2d()
    test_num = 10000
    ######################
    #Test 2D value
    ######################
    @info "testing 2D value"
    #prepare interpolation
    f1 = convert(Array{Float64},[func1_2d(x,y) for x = 1:1000, y = 1:1000])
    f2 = convert(Array{Float64},[func2_2d(x,y) for x = 1:1000, y = 1:1000])
    f3 = convert(Array{Float64},[func3_2d(x,y) for x = 1:1000, y = 1:1000])
    f4 = convert(Array{Float64},[norm(func4_2d(x,y))^2 for x = 1:1000, y = 1:1000])
    f2d = f1.+f2.+f3.+f4
    f_itp_2d = interpolate(f2d, BSpline(Cubic(Line(Interpolations.OnGrid()))))
    #prepare field
    float_sf1_2d = func2field(ScalarField{Float64,2},func1_2d,repeat([1.0],2),repeat([0.5],2),repeat([998.0],2),name = "float_sf1")
    float_sf2_2d = func2field(ScalarField{Float64,2},func2_2d,repeat([2.1],2),repeat([1.0],2),repeat([998.0],2),name = "float_sf2")
    float_sf3_2d = func2field(ScalarField{Float64,2},func3_2d,repeat([1.5],2),repeat([1.2],2),repeat([998.0],2),name = "float_sf3")
    float_vf1_2d = func2field(VectorField{Complex{Float64},2},func4_2d,repeat([1.2],2),repeat([0.9],2),repeat([998.0],2),name = "float_vf1",scaling_expr=Meta.parse("t->1.0+0.0im"))
    sfn1_2d = ScalarFieldNode{2}([float_sf1_2d,float_sf3_2d,float_sf2_2d],name = ascii("sfn1"))
    sfn2_2d = ScalarFieldNode{2}([float_vf1_2d,sfn1_2d],name = ascii("sfn2"))
    align_field_tree!(sfn2_2d)
    sum_err = 0.0
    for i = 1:test_num
        x = 979.0*rand()+10.0
        y = 979.0*rand()+10.0
        ref = f_itp_2d(x,y)
        result = value([x,y],0.0,sfn2_2d)
        sum_err += abs(ref-result)/abs(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2
    ######################
    #Test 2D gradient
    ######################
    grad = zeros(Float64,4)
    sum_err = 0.0
    for i = 1:test_num
        x = 970.0*rand()+10.0
        y = 970.0*rand()+10.0
        ref = -1.0*Interpolations.gradient(f_itp_2d,x,y)*KB/M_CS
        gradient!(0.0,[x,y,0.0,0.0],grad,sfn2_2d)
        result = grad[3:4]
        sum_err += norm(ref-result)/norm(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2
    cleanupfield!(sfn2_2d)

    ######################
    #Test 2D value for fields_array
    ######################
    @info "testing 2D value for fields array"
    #prepare field
    float_sf1_2d = func2field(ScalarField{Float64,2},func1_2d,repeat([1.0],2),repeat([0.5],2),repeat([998.0],2),name = "float_sf1")
    float_sf2_2d = func2field(ScalarField{Float64,2},func2_2d,repeat([2.1],2),repeat([1.0],2),repeat([998.0],2),name = "float_sf2")
    float_sf3_2d = func2field(ScalarField{Float64,2},func3_2d,repeat([1.5],2),repeat([1.2],2),repeat([998.0],2),name = "float_sf3")
    float_vf1_2d = func2field(VectorField{Complex{Float64},2},func4_2d,repeat([1.2],2),repeat([0.9],2),repeat([998.0],2),name = "float_vf1",scaling_expr=Meta.parse("t->1.0+0.0im"))
    sfn1_2d = ScalarFieldNode{2}([float_sf1_2d,float_sf3_2d],name = ascii("sfn1"))
    sfn2_2d = ScalarFieldNode{2}([float_vf1_2d,float_sf2_2d],name = ascii("sfn2"))
    sfn_arr_2d = Vector{ScalarFieldNode{2}}([sfn1_2d,sfn2_2d])
    align_field_arr!(sfn_arr_2d)
    sum_err = 0.0
    for i = 1:test_num
        x = 979.0*rand()+10.0
        y = 979.0*rand()+10.0
        ref = f_itp_2d(x,y)
        result = value([x,y],0.0,sfn_arr_2d)
        sum_err += abs(ref-result)/abs(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2
    ######################
    #Test 2D gradient for fields_array
    ######################
    grad = zeros(Float64,4)
    sum_err = 0.0
    for i = 1:test_num
        x = 970.0*rand()+10.0
        y = 970.0*rand()+10.0
        ref = -1.0*Interpolations.gradient(f_itp_2d,x,y)*KB/M_CS
        gradient!(0.0,[x,y,0.0,0.0],grad,sfn_arr_2d)
        result = grad[3:4]
        sum_err += norm(ref-result)/norm(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2
    cleanupfield!(sfn_arr_2d)

end

function test_gradient_3d()
    test_num = 10000
    ######################
    #Test 3D value
    ######################
    @info "testing 3D value"
    #prepare interpolation
    @info "init f1_s"
    f1_s = SharedArray{Float64}((1000,1000,100), init = S -> S[SharedArrays.localindices(S)] = map(x->func1_3d(Tuple(CartesianIndices(S)[x])...),SharedArrays.localindices(S)))
    @info "init f2_s"
    f2_s = SharedArray{Float64}((1000,1000,100), init = S -> S[SharedArrays.localindices(S)] = map(x->func2_3d(Tuple(CartesianIndices(S)[x])...),SharedArrays.localindices(S)))
    @info "init f3_s"
    f3_s = SharedArray{Float64}((1000,1000,100), init = S -> S[SharedArrays.localindices(S)] = map(x->func3_3d(Tuple(CartesianIndices(S)[x])...),SharedArrays.localindices(S)))
    @info "init f4_s"
    f4_s = SharedArray{Complex{Float64}}((3,1000,1000,100), init=sa_init)
    # for i = 1:100
    #     x = rand(1:1000)
    #     y = rand(1:1000)
    #     z = rand(1:100)
    #     @assert f4_s[:,x,y,z] == func4_3d(x,y,z)
    # end
    @info "init f4_norm2_s"
    f4norm2_s = SharedArray{Float64}((1000,1000,100), init = S -> S[SharedArrays.localindices(S)] = map(x->norm(func4_3d(Tuple(CartesianIndices(S)[x])...))^2,SharedArrays.localindices(S)))
    f3d = f1_s.+f2_s.+f3_s.+f4norm2_s
    f_itp_3d = interpolate(f3d, BSpline(Cubic(Line(Interpolations.OnGrid()))))
    ### Debug
    f4_itp_3d = interpolate(f4norm2_s, BSpline(Cubic(Line(Interpolations.OnGrid()))))
    f123_itp_3d = interpolate(f1_s.+f2_s.+f3_s, BSpline(Cubic(Line(Interpolations.OnGrid()))))
    ###
    @info "build field array"
    float_sf1_3d = func2field(ScalarField{Float64,3},func1_3d,repeat([1.0],3),repeat([0.5],3),[999.0,999.0,99.0],name = "float_sf1")
    float_sf2_3d = func2field(ScalarField{Float64,3},func2_3d,repeat([2.1],3),repeat([1.0],3),[999.0,999.0,99.0],name = "float_sf2")
    float_sf3_3d = func2field(ScalarField{Float64,3},func3_3d,repeat([2.5],3),repeat([1.2],3),[999.0,999.0,99.0],name = "float_sf3")
    float_vf1_3d = func2field(VectorField{Complex{Float64},3},func4_3d,repeat([1.2],3),repeat([0.9],3),[999.0,999.0,99.0],name = "float_vf1",scaling_expr=Meta.parse("t->1.0+0.0im"))
    sfn1_3d = ScalarFieldNode{3}([float_sf1_3d,float_sf3_3d],name = ascii("sfn1"))
    sfn2_3d = ScalarFieldNode{3}([float_vf1_3d,sfn1_3d,float_sf2_3d],name = ascii("sfn2"))
    align_field_tree!(sfn2_3d)
    sum_err = 0.0
    for i = 1:test_num
        x = 970.0*rand()+10.0
        y = 970.0*rand()+10.0
        z = 90.0*rand()+5.0
        ref = f_itp_3d(x,y,z)
        result = value([x,y,z],0.0,sfn2_3d)
        sum_err += abs(ref-result)/abs(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2

    ######################
    #Test 3D gradient
    ######################
    grad = zeros(Float64,6)
    sum_err = 0.0
    @time for i = 1:test_num
        x = 970.0*rand()+10.0
        y = 970.0*rand()+10.0
        z = 90.0*rand()+5.0
        ref = -1.0*Interpolations.gradient(f_itp_3d,x,y,z)*KB/M_CS
        gradient!(0.0,[x,y,z,0.0,0.0,0.0],grad,sfn2_3d)
        result = grad[4:6]
        sum_err += norm(ref-result)/norm(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2
    cleanupfield!(sfn2_3d)

    ######################
    #Test 3D value for fields_array
    ######################
    @info "testing 3D value for fields array"
    float_sf1_3d = func2field(ScalarField{Float64,3},func1_3d,repeat([1.0],3),repeat([0.5],3),[999.0,999.0,99.0],name = "float_sf1")
    float_sf2_3d = func2field(ScalarField{Float64,3},func2_3d,repeat([2.1],3),repeat([1.0],3),[999.0,999.0,99.0],name = "float_sf2")
    float_sf3_3d = func2field(ScalarField{Float64,3},func3_3d,repeat([2.5],3),repeat([1.2],3),[999.0,999.0,99.0],name = "float_sf3")
    float_vf1_3d = func2field(VectorField{Complex{Float64},3},func4_3d,repeat([1.2],3),repeat([0.9],3),[999.0,999.0,99.0],name = "float_vf1",scaling_expr=Meta.parse("t->1.0+0.0im"))
    sfn1_3d = ScalarFieldNode{3}([float_sf1_3d,float_sf3_3d],name = ascii("sfn1"))
    sfn2_3d = ScalarFieldNode{3}([float_vf1_3d,float_sf2_3d],name = ascii("sfn2"))
    sfn_arr_3d = Vector{ScalarFieldNode{3}}([sfn1_3d,sfn2_3d])
    align_field_arr!(sfn_arr_3d)
    sum_err = 0.0
    for i = 1:test_num
        x = 979.0*rand()+10.0
        y = 979.0*rand()+10.0
        z = 90.0*rand()+5.0
        ref = f_itp_3d(x,y,z)
        result = value([x,y,z],0.0,sfn_arr_3d)
        sum_err += abs(ref-result)/abs(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2
    ######################
    #Test 3D gradient for fields_array
    #####################
    @info "testing 3D gradient for fields array"
    grad = zeros(Float64,6)
    sum_err = 0.0
    for i = 1:test_num
        x = 970.0*rand()+10.0
        y = 970.0*rand()+10.0
        z = 90.0*rand()+5.0
        ref = -1.0*Interpolations.gradient(f_itp_3d,x,y,z)*KB/M_CS
        gradient!(0.0,[x,y,z,0.0,0.0,0.0],grad,sfn_arr_3d)
        result = grad[4:6]
        sum_err += norm(ref-result)/norm(ref)
    end
    err = sum_err/test_num
    @info "err: "*string(err)
    @test err<1e-2
    cleanupfield!(sfn_arr_3d)

end

function test_gradient()
    @testset "test field gradient" begin
        test_gradient_2d()
        test_gradient_3d()
    end
end
