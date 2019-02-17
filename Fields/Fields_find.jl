function find_field(criteria::Function,f_arr::Vector{ScalarFieldNode})
    #currently only find one result
    for sfn in f_arr
        try
            return find_field(criteria,sfn)
        catch
        end
    end
    error("Can't find field.")
end

function find_field(criteria::Function,f::T) where {T <: FieldNode}
    #currently only find one result
    if criteria(f)
        return f
    else
        for ff in f.fields
            if find_field_bool(criteria,ff)
                return find_field(criteria,ff)
            end
        end
        error("Can't find field.")
    end
end

find_field(criteria::Function,f::T) where {T<:Union{ScalarField, VectorField}} = criteria(f) ? f : error("Can't find field.")
#find_field_bool(criteria::Function)=find_field_bool(criteria,fields)
find_field_bool(criteria::Function,f::T) where {T <: FieldNode} =criteria(f) ? true : any(map(f->find_field_bool(criteria,f),f.fields))
find_field_bool(criteria::Function,f::T) where {T <:Union{ScalarField, VectorField}}=criteria(f)

function test_find()
    @testset "Fields.find_field test" begin
        for dim in [2,3]
            float_sf1 = zero_field(ScalarField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim),name = "float_sf1")
            float_sf2 = zero_field(ScalarField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim),name = "float_sf2")
            complex_sf = zero_field(ScalarField{Complex{Float64},dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim),name = "complex_sf")
            sfn1 = ScalarFieldNode{dim}([float_sf1,complex_sf],name = ascii("sfn1"))
            @test find_field(x->x.name==ascii("sfn1"),sfn1) == sfn1
            @test find_field(x->x.name==ascii("float_sf1"),sfn1) == float_sf1
            @test find_field(x->x.name==ascii("complex_sf"),sfn1) == complex_sf
            sfn2 = ScalarFieldNode{dim}([float_sf2,sfn1],name = ascii("sfn2"))
            @test find_field(x->x.name==ascii("sfn2"),sfn2) == sfn2
            @test find_field(x->x.name==ascii("sfn1"),sfn2) == sfn1
            @test find_field(x->x.name==ascii("float_sf1"),sfn2) == float_sf1
            @test find_field(x->x.name==ascii("float_sf2"),sfn2) == float_sf2
            @test find_field(x->x.name==ascii("complex_sf"),sfn2) == complex_sf
            float_vf1 = zero_field(VectorField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim),name="float_vf1")
            float_vf2 = zero_field(VectorField{Float64,dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim),name="float_vf2")
            complex_vf = zero_field(VectorField{Complex{Float64},dim},repeat([1.0],dim),repeat([0.0],dim),repeat([100.0],dim),name="complex_vf")
            vfn1 = VectorFieldNode{dim}([float_vf1,complex_vf],name=ascii("vfn1"))
            @test find_field(x->x.name==ascii("vfn1"),vfn1) == vfn1
            @test find_field(x->x.name==ascii("float_vf1"),vfn1) == float_vf1
            @test find_field(x->x.name==ascii("complex_vf"),vfn1) == complex_vf
            vfn2 = VectorFieldNode{dim}([float_vf2,vfn1],name=ascii("vfn2"))
            @test find_field(x->x.name==ascii("vfn2"),vfn2) == vfn2
            @test find_field(x->x.name==ascii("vfn1"),vfn2) == vfn1
            @test find_field(x->x.name==ascii("float_vf1"),vfn2) == float_vf1
            @test find_field(x->x.name==ascii("float_vf2"),vfn2) == float_vf2
            @test find_field(x->x.name==ascii("complex_vf"),vfn2) == complex_vf
            @test_throws ErrorException find_field(x->x.name==ascii("dummy"),vfn2) == complex_vf

            sfn3 = ScalarFieldNode{dim}([vfn2],name = ascii("sfn3"))
            sfn_arr = Vector{ScalarFieldNode}([sfn2,sfn3])
            #test vector of sfns
            @test find_field(x->x.name==ascii("vfn2"),sfn_arr) == vfn2
            @test find_field(x->x.name==ascii("sfn3"),sfn_arr) == sfn3
            @test find_field(x->x.name==ascii("sfn2"),sfn_arr) == sfn2
            @test find_field(x->x.name==ascii("complex_sf"),sfn_arr) == complex_sf
            @test find_field(x->x.name==ascii("float_sf1"),sfn_arr) == float_sf1
        end
    end
end
