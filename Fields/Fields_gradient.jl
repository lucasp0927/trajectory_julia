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
    copy!(f.sample,sub(f.field,pidx[1]:pidx[2],pidx[3]:pidx[4],pidx[5]:pidx[6]))
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
    f.s::Complex{Float64} = convert(Complex{Float64},(f.scaling::Function)(t))
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
    f.s::Complex{Float64} = convert(Complex{Float64},(f.scaling::Function)(t))
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
    value(pos,t,fields::ScalarFieldNode)
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
        grad[4:6] = -1.0*itp_tricubic_grad((sfn::ScalarFieldNode).sample,x,res)*KB/M_CS
    end
end

################################################
#Test
################################################
function sa_init(S::SharedArray)
    for id in Base.localindexes(S)
        S[id] = func4_3d2(ind2sub(S,id)...)
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

function test_gradient()
    #test with both scalar and vector field
    test_num = 10000

    ######################
    #Test 2D value
    ######################
    Lumberjack.info("testing 2D value")
    #prepare interpolation
    f1 = convert(Array{Float64},[func1_2d(x,y) for x = 1:1000, y = 1:1000])
    f2 = convert(Array{Float64},[func2_2d(x,y) for x = 1:1000, y = 1:1000])
    f3 = convert(Array{Float64},[func3_2d(x,y) for x = 1:1000, y = 1:1000])
    f4 = convert(Array{Float64},[norm(func4_2d(x,y))^2 for x = 1:1000, y = 1:1000])
    f2d = f1.+f2.+f3.+f4
    f_itp_2d = interpolate(f2d, BSpline(Cubic(Line())), OnGrid())
    #prepare field
    float_sf1_2d = func2field(ScalarField{Float64,2},func1_2d,repmat([2500],2),repmat([1.0],2),repmat([999.0],2),name = "float_sf1")
    float_sf2_2d = func2field(ScalarField{Float64,2},func2_2d,repmat([2000],2),repmat([1.0],2),repmat([999.0],2),name = "float_sf2")
    float_sf3_2d = func2field(ScalarField{Float64,2},func3_2d,repmat([1001],2),repmat([1.0],2),repmat([999.0],2),name = "float_sf3")
    float_vf1_2d = func2field(VectorField{Complex{Float64},2},func4_2d,repmat([1001],2),repmat([1.0],2),repmat([999.0],2),name = "float_vf1")
    sfn1_2d = ScalarFieldNode{2}([float_sf1_2d,float_sf3_2d,float_sf2_2d],name = ascii("sfn1"))
    sfn2_2d = ScalarFieldNode{2}([float_vf1_2d,sfn1_2d],name = ascii("sfn2"))
    align_field_tree!(sfn2_2d)
    sum_err = 0.0
    @time    for i = 1:test_num
        x = 996.0*rand()+2.0
        y = 996.0*rand()+2.0
        ref = f_itp_2d[x,y]
        result = value([x,y],0.0,sfn2_2d)
        sum_err += abs(ref-result)/abs(ref)
    end
    err = sum_err/test_num
    Lumberjack.info("err: ",string(err))
    @test err<3e-2
    ######################
    #Test 3D value
    ######################
    Lumberjack.info("testing 3D value")
    #prepare interpolation
    f1_s = SharedArray(Float64, (1000,1000,100), init = S -> S[Base.localindexes(S)] = map(x->func1_3d(ind2sub(S,x)...),Base.localindexes(S)))
    f2_s = SharedArray(Float64, (1000,1000,100), init = S -> S[Base.localindexes(S)] = map(x->func2_3d(ind2sub(S,x)...),Base.localindexes(S)))
    f3_s = SharedArray(Float64, (1000,1000,100), init = S -> S[Base.localindexes(S)] = map(x->func3_3d(ind2sub(S,x)...),Base.localindexes(S)))
    f4_s = SharedArray(Complex{Float64}, (3,1000,1000,100), init=sa_init )
    for i = 1:100
        x = rand(1:1000)
        y = rand(1:1000)
        z = rand(1:100)
        @test f4_s[:,x,y,z] == func4_3d(x,y,z)
    end
    #=
    f4_s = SharedArray(Complex{Float64},3,1000,1000,100)
    @time    for x = 1:1000,y=1:1000,z=1:100
    f4_s[:,x,y,z] = func4_3d(x,y,z)
    end
    =#
    println(5)
    f4norm2_s = SharedArray(Float64, (1000,1000,100), init = S -> S[Base.localindexes(S)] = map(x->norm(func4_3d(ind2sub(S,x)...))^2,Base.localindexes(S)))
    f3d = f1_s.+f2_s.+f3_s.+f4norm2_s
    f_itp_3d = interpolate(f3d, BSpline(Cubic(Line())), OnGrid())
    #prepare field
    float_sf1_3d = ScalarField{Float64,3}(f1_s,[1.0,1.0,1.0],[999.0,999.0,99.0],name="float_sf1")
    float_sf2_3d = ScalarField{Float64,3}(f2_s,[1.0,1.0,1.0],[999.0,999.0,99.0],name="float_sf2")
    float_sf3_3d = ScalarField{Float64,3}(f3_s,[1.0,1.0,1.0],[999.0,999.0,99.0],name="float_sf3")
    float_vf1_3d = VectorField{Complex{Float64},3}(f4_s,[1.0,1.0,1.0],[999.0,999.0,99.0],name="float_vf1")
    sfn1_3d = ScalarFieldNode{3}([float_sf1_3d,float_sf3_3d],name = ascii("sfn1"))
    sfn2_3d = ScalarFieldNode{3}([float_vf1_3d,sfn1_3d,float_sf2_3d],name = ascii("sfn2"))
#    sfn2_3d = ScalarFieldNode{3}([float_vf1_3d,sfn1_3d],name = ascii("sfn2"))
    align_field_tree!(sfn2_3d)
    sum_err = 0.0
    @time    for i = 1:test_num
        x = 996.0*rand()+2.0
        y = 996.0*rand()+2.0
        z = 96.0*rand()+2.0
        ref = f_itp_3d[x,y,z]
        result = value([x,y,z],0.0,sfn2_3d)
        sum_err += abs(ref-result)/abs(ref)
    end
    err = sum_err/test_num
    Lumberjack.info("err: ",string(err))
    @test err<3e-2
    ######################
    #Test 2D gradient
    ######################
    grad = zeros(Float64,4)
    sum_err = 0.0
    @time    for i = 1:test_num
        x = 996.0*rand()+2.0
        y = 996.0*rand()+2.0
        ref = -1.0*gradient(f_itp_2d,x,y)*KB/M_CS
        gradient!(0.0,[x,y,0.0,0.0],grad,sfn2_2d)
        result = grad[3:4]
        sum_err += norm(ref-result)/norm(ref)
    end
    err = sum_err/test_num
    Lumberjack.info("err: ",string(err))
    @test err<3e-2
    ######################
    #Test 3D gradient
    ######################
    grad = zeros(Float64,6)
    sum_err = 0.0
    @time    for i = 1:test_num
        x = 998.0*rand()+1.0
        y = 998.0*rand()+1.0
        z = 98.0*rand()+1.0
        ref = -1.0*gradient(f_itp_3d,x,y,z)*KB/M_CS
        gradient!(0.0,[x,y,z,0.0,0.0,0.0],grad,sfn2_3d)
        result = grad[4:6]
        sum_err += norm(ref-result)/norm(ref)
    end
    err = sum_err/test_num
    Lumberjack.info("err: ",string(err))
    @test err<3e-2
end
