function set_geometry!(f::T) where {T <: FieldNode}
    #memoize geo parameters
    geo = geometry(f)
    f.position = geo["pos"]
    f.size = geo["size"]
    f.res = geo["res"]
    for ff in f.fields
        if typeof(ff) <: FieldNode
            set_geometry!(ff)
        end
    end
end

function geometry(f::T;lowest_res=false) where {T <: FieldNode}
    #find the new optimal position and size
    #minimum(cat(2,[1,2,3],[2,3,4]),2)
    minpos::Vector{Float64} = vec(minimum(cat(map(x->geometry(x)["pos"],f.fields)...,dims=2),dims=2))
    maxpos::Vector{Float64} = vec(maximum(cat(map(x->geometry(x)["pos"].+geometry(x)["size"],f.fields)...,dims=2),dims=2))
    myminmax = (lowest_res ? maximum : minimum)::Function
    res::Vector{Float64} = vec(myminmax(cat(map(x->geometry(x)["res"],f.fields)...,dims=2),dims=2))
    @assert all(x->x!=0,res) "zero resolution!"
    return Dict("pos"=>minpos, "size"=>(maxpos-minpos), "res"=>res)
end

geometry(f::T) where {T <: Union{VectorField,ScalarField}}= Dict("pos"=>f.position,"size"=>f.size,"res"=>f.res)

function geometry()
    Dict("pos"=>fields.position,"size"=>fields.size,"res"=>fields.res)
end

function test_geometry()
    #2D
    #=
    sfn2
     |_____
     |    |
    fsf2  sfn1
          |_____
          |    |
        fsf1   csf
    =#
    for dim in [2,3]
        float_sf1 = zero_field(ScalarField{Float64,dim},repmat([201],dim), repmat([50.0],dim),repmat([50.0],dim),name = "float_sf1")
        @test geometry(float_sf1) == Dict("res"=>repmat([0.25],dim),"size"=>repmat([50.0],dim),"pos"=>repmat([50.0],dim))
        float_sf2 = zero_field(ScalarField{Float64,dim},repmat([101],dim),repmat([0.0],dim),repmat([100.0],dim),name = "float_sf2")
        @test geometry(float_sf2) == Dict("res"=>repmat([1.0],dim),"size"=>repmat([100.0],dim),"pos"=>repmat([0.0],dim))
        complex_sf = zero_field(ScalarField{Complex{Float64},dim},repmat([121],dim),repmat([20.0],dim),repmat([60.0],dim),name = "complex_sf")
        @test geometry(complex_sf) == Dict("res"=>repmat([0.5],dim),"size"=>repmat([60.0],dim),"pos"=>repmat([20.0],dim))
        sfn1 = ScalarFieldNode{dim}([float_sf1,complex_sf],name = ascii("sfn1"))
        @test geometry(sfn1) == Dict("res"=>repmat([0.5],dim),"size"=>repmat([80.0],dim),"pos"=>repmat([20.0],dim))
        sfn2 = ScalarFieldNode{dim}([float_sf2,sfn1],name = ascii("sfn2"))
        @test geometry(sfn2) == Dict("res"=>repmat([1.0],dim),"size"=>repmat([100.0],dim),"pos"=>repmat([0.0],dim))
        set_geometry!(sfn2)
        @test sfn1.res == repmat([0.5],dim)
        @test sfn1.size ==repmat([80.0],dim)
        @test sfn1.position == repmat([20.0],dim)
        @test sfn2.res == repmat([1.0],dim)
        @test sfn2.size == repmat([100.0],dim)
        @test sfn2.position == repmat([0.0],dim)
    end
    for dim in [2,3]
        float_vf1 = zero_field(VectorField{Float64,dim},repmat([201],dim), repmat([50.0],dim),repmat([50.0],dim),name = "float_vf1")
        @test geometry(float_vf1) == Dict("res"=>repmat([0.25],dim),"size"=>repmat([50.0],dim),"pos"=>repmat([50.0],dim))
        float_vf2 = zero_field(VectorField{Float64,dim},repmat([101],dim),repmat([0.0],dim),repmat([100.0],dim),name = "float_vf2")
        @test geometry(float_vf2) == Dict("res"=>repmat([1.0],dim),"size"=>repmat([100.0],dim),"pos"=>repmat([0.0],dim))
        complex_vf = zero_field(VectorField{Complex{Float64},dim},repmat([121],dim),repmat([20.0],dim),repmat([60.0],dim),name = "complex_vf")
        @test geometry(complex_vf) == Dict("res"=>repmat([0.5],dim),"size"=>repmat([60.0],dim),"pos"=>repmat([20.0],dim))
        vfn1 = VectorFieldNode{dim}([float_vf1,complex_vf],name = ascii("vfn1"))
        @test geometry(vfn1) == Dict("res"=>repmat([0.5],dim),"size"=>repmat([80.0],dim),"pos"=>repmat([20.0],dim))
        vfn2 = VectorFieldNode{dim}([float_vf2,vfn1],name = ascii("vfn2"))
        @test geometry(vfn2) == Dict("res"=>repmat([1.0],dim),"size"=>repmat([100.0],dim),"pos"=>repmat([0.0],dim))
        set_geometry!(vfn2)
        @test vfn1.res == repmat([0.5],dim)
        @test vfn1.size ==repmat([80.0],dim)
        @test vfn1.position == repmat([20.0],dim)
        @test vfn2.res == repmat([1.0],dim)
        @test vfn2.size == repmat([100.0],dim)
        @test vfn2.position == repmat([0.0],dim)
    end
end
