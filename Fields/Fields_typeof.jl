####type of field function
typeoffield(f::VectorField{T,N}) where {T<:ComplexOrFloat,N}= T

typeoffield(f::ScalarField{T,N}) where {T<:ComplexOrFloat,N} = T

function typeoffield(f::VectorFieldNode)
    f_member_type = map(typeoffield,f.fields)
    @assert all(x->x<:ComplexOrFloat,f_member_type) "Fields have to be either Complex{Float64} or Float64"
    return typeintersect(reduce((x,y)->Union{x,y},f_member_type),ComplexOrFloat)==Float64 ? Float64 : Complex{Float64}
end

function typeoffield(f::ScalarFieldNode)
    f_member_type = map(typeoffield,f.fields)
    @assert all(x->x<:ComplexOrFloat,f_member_type) "Fields have to be either Complex{Float64} or Float64"
    return typeintersect(reduce((x,y)->Union{x,y},f_member_type),ComplexOrFloat)==Float64 ? Float64 : Complex{Float64}
end

function set_typeof!(f::T) where {T <:FieldNode}
    #memoize typeof
    f.typeof = typeoffield(f)
    for ff in f.fields
        if typeof(ff) <: FieldNode
            set_typeof!(ff)
        end
    end
end

function test_typeof()
    @testset "test fields typeof" begin
        for dim in [2,3]
            float_sf1 = zero_field(ScalarField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim))
            float_sf2 = zero_field(ScalarField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim))
            complex_sf = zero_field(ScalarField{Complex{Float64},dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim))
            sfn1 = ScalarFieldNode{dim}([float_sf1,complex_sf])
            set_typeof!(sfn1)
            @test sfn1.typeof == Complex{Float64}
            sfn2 = ScalarFieldNode{dim}([float_sf1,float_sf2])
            set_typeof!(sfn2)
            @test sfn2.typeof == Float64
            sfn3 = ScalarFieldNode{dim}([sfn1,sfn2])
            set_typeof!(sfn3)
            @test sfn3.typeof == Complex{Float64}
            float_vf1 = zero_field(VectorField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim))
            float_vf2 = zero_field(VectorField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim))
            complex_vf = zero_field(VectorField{Complex{Float64},dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim))
            vfn1 = VectorFieldNode{dim}([float_vf1,complex_vf])
            set_typeof!(vfn1)
            @test vfn1.typeof == Complex{Float64}
            vfn2 = VectorFieldNode{dim}([float_vf1,float_vf2])
            set_typeof!(vfn2)
            @test vfn2.typeof == Float64
            vfn3 = VectorFieldNode{dim}([vfn1,vfn2])
            set_typeof!(vfn3)
            @test vfn3.typeof == Complex{Float64}
        end
    end
end
