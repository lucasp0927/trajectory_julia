function align_field_arr!(f_arr::Vector{ScalarFieldNode{N}}) where N
    map(align_field_tree!,f_arr)
end

function align_field_tree!(f::T) where T<:FieldNode
    set_geometry!(f)
    #debug("align top field node $(f.name)")
    @debug "align top field node "*f.name
    unalign_geo = geometry(f)
    @debug "unaligned position $(unalign_geo["pos"])"
    @debug "unaligned size $(unalign_geo["size"])"
    @debug "unaligned resolution $(unalign_geo["res"])"
    arr_sz = floor.(Integer,unalign_geo["size"]./unalign_geo["res"])
    new_sz = arr_sz.*unalign_geo["res"]
    align_geo = unalign_geo
    align_geo["size"] = new_sz
    @debug "aligned position $(align_geo["pos"])"
    @debug "aligned size $(align_geo["size"])"
    @debug "aligned resolution $(align_geo["res"])"
    align_field!(f,align_geo["res"],align_geo["pos"])
    set_geometry!(f)
    set_typeof!(f)
    Base.GC.gc()
end

function align_field!(f::T,res::Vector{Float64},pos::Vector{Float64}) where T<:FieldNode
    @debug "align FieldNode $(f.name)"
    @assert length(res) == length(pos) "dimension mismatch!"
    for ff in f.fields
        align_field!(ff,res,pos)
    end
end

function align_field!(f::ScalarFieldFunc{T, N}, res::Vector{Float64}, pos::Vector{Float64}) where {T <: ComplexOrFloat, N}
    @info "align ScalarFieldFunc $(f.name)"
end


function align_field!(f::ScalarField{T, N}, res::Vector{Float64}, pos::Vector{Float64}) where {T <: ComplexOrFloat, N}
    @info "align ScalarField $(f.name)"
    @assert length(res) == length(pos) "dimension mismatch!"
    unalign_geo = geometry(f)
    new_pos = ceil.((unalign_geo["pos"].-pos)./res).*res.+pos
    @assert all(map(>=, new_pos, unalign_geo["pos"]))
    @assert all(map(<, new_pos-unalign_geo["pos"],res))
    old_end = unalign_geo["pos"].+unalign_geo["size"]
    new_end = floor.((old_end.-pos)./res).*res.+pos
    @assert all(map(<=, new_end, old_end))
    @assert all(map(<, old_end.-new_end,res))
    new_size = new_end.-new_pos
    align_geo = Dict("pos"=>new_pos,"size"=>new_size,"res"=>res)
    ##### interpolate field
    if (unalign_geo["pos"] == align_geo["pos"])&&(unalign_geo["res"] == align_geo["res"])
        @info "no need to align!"
    else
        new_arr_size = round.(Int64,new_size ./ res).+one(Int64)
        new_field = SharedArray{T}(new_arr_size...)
        itp_field!(new_field,f.field,unalign_geo,align_geo)
        @assert collect(size(new_field)) == new_arr_size
        setfield!(f,new_field,align_geo["pos"],align_geo["size"],scaling=f.scaling)
        #check
        @assert f.res ≈ align_geo["res"]
    end
end

function align_field!(f::VectorField{T, N}, res::Vector{Float64}, pos::Vector{Float64}) where {T <: ComplexOrFloat, N}
    @info "align VectorField $(f.name)"
    @assert length(res) == length(pos) "dimension mismatch!"
    unalign_geo = geometry(f)
    new_pos = ceil.((unalign_geo["pos"].-pos)./res).*res.+pos
    @assert all(map(>=, new_pos, unalign_geo["pos"]))
    @assert all(map(<, new_pos-unalign_geo["pos"],res))
    old_end = unalign_geo["pos"].+unalign_geo["size"]
    new_end = floor.((old_end.-pos)./res).*res.+pos
    @assert all(map(<=, new_end, old_end))
    @assert all(map(<, old_end.-new_end,res))
    new_size = new_end.-new_pos
    align_geo = Dict("pos"=>new_pos,"size"=>new_size,"res"=>res)
    ##### interpolate field
    if (unalign_geo["pos"] == align_geo["pos"])&&(unalign_geo["res"] == align_geo["res"])
        @info "no need to align!"
    else
        #loop over three components
        new_arr_size = round.(Int64,new_size ./ res).+one(Int64)
        new_field = SharedArray{T}(3,new_arr_size...)
        for i = 1:3
            itp_field!(myslice(new_field,i),myslice(f.field,i),unalign_geo,align_geo)
        end
        @assert collect(size(new_field)[2:end]) == new_arr_size
        setfield!(f,new_field,align_geo["pos"],align_geo["size"],scaling=f.scaling)
        #check
        @assert f.res ≈ align_geo["res"]
    end
end

#TODO clean this up with meta programming?
@generated function myslice(A::Array{T, N}, i::Integer) where {T <: ComplexOrFloat, N}
    ex_str = "view(A,i"*repeat(",:",N-1)*")"
    Meta.parse(ex_str)
end

@generated function myslice(A::SharedArray{T, N}, i::Integer) where {T <: ComplexOrFloat, N}
    ex_str = "view(A,i"*repeat(",:",N-1)*")"
    Meta.parse(ex_str)
end

@inbounds function transform_coordinate(apos::Vector{Float64},ares::Vector{Float64},uapos::Vector{Float64},uares::Vector{Float64},index::Vector{Int64})
    #2D
    #position at current index
    #@devec    idx_pos = apos.+ ares.*(index-1)
    idx_pos = apos.+ ares.*(index.-1)
    #calculate the old field index
    #@devec    old_idx = ((idx_pos.-uapos)./uares)+1
    old_idx = ((idx_pos.-uapos)./uares).+1
    #test
    #@devec    un_idx_pos = uapos.+ uares.*(old_idx-1)
    un_idx_pos = uapos.+ uares.*(old_idx.-1)
    @assert   un_idx_pos ≈ idx_pos "coordinate transform check"
    return old_idx
end

@generated function itp_field!(new_field::T, old_field::K, unalign_geo, align_geo) where {T <: AbstractArray, K <: AbstractArray}
    N::Int64 = ndims(new_field)
    quote
#        old_field_itp = interpolate(old_field, BSpline(Cubic(Flat())), OnGrid())
        old_field_itp = extrapolate(interpolate(old_field, BSpline(Cubic(Flat(Interpolations.OnGrid())))),Line())
        apos::Vector{Float64} = align_geo["pos"]
        ares::Vector{Float64} = align_geo["res"]
        uapos::Vector{Float64} = unalign_geo["pos"]
        uares::Vector{Float64} = unalign_geo["res"]
        start_idx = transform_coordinate(apos,ares,uapos,uares,ones(Int64,$N))
        end_idx = transform_coordinate(apos,ares,uapos,uares,collect(size(new_field)))
        @nexprs $N j-> x_j = range(start_idx[j],stop=end_idx[j],length=size(new_field,j))
        @nloops $N i new_field begin
            #(@nref $N new_field i) = (@nref $N old_field_itp j->x_j[i_j])
            (@nref $N new_field i) = (@ncall $N old_field_itp j->x_j[i_j])
        end
        old_field_itp = 0
    end
end

##### test functions
    function test_on_grid(f,fn,res)
        sizerem = map(rem,(f.position .- fn.position),res)
        if !(all(map(≈,sizerem,res) .| map(x->x<1e-10,sizerem)))
            return false
        end
        sizerem = map(rem,(f.position .+ f.size .- fn.position),res)
        if !(all(map(≈,sizerem,res) .| map(x->x<1e-10,sizerem)))
            return false
        end
        return true
    end

    function test_interpolate(f::ScalarField{T, 2}, func) where T <: ComplexOrFloat
        @info "testing $(f.name)"
        sz = size(f.field)
        pos = f.position
        res = f.res
        sum = 0.0
        for i = 1:sz[1], j = 1:sz[2]
            x_test = (i-1)*res[1]+pos[1]
            y_test = (j-1)*res[2]+pos[2]
            if abs(func(x_test,y_test)) != 0.0
                sum += abs(func(x_test,y_test) - f.field[i,j])/abs(func(x_test,y_test))
            end
        end
        err = sum/(sz[1]*sz[2])
        @info "err: "*string(err)
        return err
    end

    function test_interpolate(f::VectorField{T, 2}, func) where T <: ComplexOrFloat
        @info "testing $(f.name)"
        sz = size(f.field)
        pos = f.position
        res = f.res
        sum = 0.0
        for i = 1:sz[2], j = 1:sz[3]
            x_test = (i-1)*res[1]+pos[1]
            y_test = (j-1)*res[2]+pos[2]
            if norm(func(x_test,y_test)) != 0.0
                sum += norm(func(x_test,y_test) .- f.field[:,i,j])/norm(func(x_test,y_test))
            end
        end
        err = sum/(3*sz[2]*sz[3])
        @info "err: "*string(err)
        return err
    end

    function test_interpolate(f::ScalarField{T, 3}, func) where T <: ComplexOrFloat
        @info "testing $(f.name)"
        sz = size(f.field)
        pos = f.position
        res = f.res
        sum = 0.0
        for i = 1:sz[1], j = 1:sz[2], k = 1:sz[3]
            x_test = (i-1)*res[1]+pos[1]
            y_test = (j-1)*res[2]+pos[2]
            z_test = (k-1)*res[3]+pos[3]
            if abs(func(x_test,y_test,z_test)) != 0.0
                sum += abs(func(x_test,y_test,z_test) - f.field[i,j,k])/abs(func(x_test,y_test,z_test))
            end
        end
        err = sum/(sz[1]*sz[2]*sz[3])
        @info "err: "*string(err)
        return err
    end

    function test_interpolate(f::VectorField{T, 3}, func) where T <: ComplexOrFloat
        @info "testing $(f.name)"
        sz = size(f.field)
        pos = f.position
        res = f.res
        sum = 0.0
        for i = 1:sz[2], j = 1:sz[3], k = 1:sz[4]
            x_test = (i-1)*res[1]+pos[1]
            y_test = (j-1)*res[2]+pos[2]
            z_test = (k-1)*res[3]+pos[3]
            if norm(func(x_test,y_test,z_test)) != 0.0
                sum += norm(func(x_test,y_test,z_test) .- f.field[:,i,j,k])/norm(func(x_test,y_test,z_test))
            end
        end
        err = sum/(sz[2]*sz[3]*sz[4]*3)
        @info "err: "*string(err)
        return err
    end

function test_align()
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
    #2D scalar
    @testset "test Fields align" begin
        @info "2D Scalar"
        func1 = (x,y)->(sin(x/10.0)+cos(y/10.0))*exp(-((x-75.0)^2+(y-75.0)^2)/20^2)
        func2 = (x,y)->(sin(x/10.0)+cos(y/10.0))*exp(-((x-50.0)^2+(y-50.0)^2)/30^2)
        func3 = (x,y)->(sin(x/10.0)+cos(y/10.0))*exp(-((x-50.0)^2+(y-50.0)^2)/20^2)
        float_sf1 = func2field(ScalarField{Float64,2},func1,repeat([0.05],2), repeat([50.1],2),repeat([50.0],2),name = "float_sf1")
        float_sf2 = func2field(ScalarField{Float64,2},func2,repeat([1.0],2),repeat([0.0],2),repeat([100.0],2),name = "float_sf2")
        complex_sf = func2field(ScalarField{Complex{Float64},2},func3,repeat([0.06],2),repeat([19.9],2),repeat([60.0],2),name = "complex_sf")
        sfn1 = ScalarFieldNode{2}([float_sf1,complex_sf],name = ascii("sfn1"))
        sfn2 = ScalarFieldNode{2}([float_sf2,sfn1],name = ascii("sfn2"))
        align_field_tree!(sfn2)
        # test if all points are on grid
        @test float_sf1.res ≈ float_sf2.res ≈ complex_sf.res ≈ sfn1.res ≈ sfn2.res
        res = sfn2.res
        @test test_on_grid(float_sf1,sfn2,res)
        @test test_on_grid(float_sf2,sfn2,res)
        @test test_on_grid(complex_sf,sfn2,res)
        @test test_on_grid(sfn1,sfn2,res)
        @test test_on_grid(sfn2,sfn2,res)
        # test interpolation result
        @test test_interpolate(float_sf1,func1) <= 1e-2
        @test test_interpolate(float_sf2,func2) <= 1e-2
        @test test_interpolate(complex_sf,func3)<= 1e-2


        #2D vector
        @info "2D Vector"
        func1 = (x,y)->repeat([(sin(x/10.0)+cos(y/10.0))*exp(-((x-75.0)^2+(y-75.0)^2)/20^2)],3)
        func2 = (x,y)->repeat([(sin(x/10.0)+cos(y/10.0))*exp(-((x-50.0)^2+(y-50.0)^2)/30^2)],3)
        func3 = (x,y)->repeat([(sin(x/10.0)+cos(y/10.0))*exp(-((x-50.0)^2+(y-50.0)^2)/20^2)],3)
        float_vf1 = func2field(VectorField{Float64,2},func1,repeat([0.05],2), repeat([50.1],2),repeat([50.0],2),name = "float_vf1")
        float_vf2 = func2field(VectorField{Float64,2},func2,repeat([1.0],2),repeat([0.0],2),repeat([100.0],2),name = "float_vf2")
        complex_vf = func2field(VectorField{Complex{Float64},2},func3,repeat([0.06],2),repeat([19.9],2),repeat([60.0],2),name = "complex_vf")
        vfn1 = VectorFieldNode{2}([float_vf1,complex_vf],name = ascii("vfn1"))
        vfn2 = VectorFieldNode{2}([float_vf2,vfn1],name = ascii("vfn2"))
        align_field_tree!(vfn2)
        # test if all points ar on grid
        @test float_vf1.res ≈ float_vf2.res ≈ complex_vf.res ≈ vfn1.res ≈ vfn2.res
        res = vfn2.res
        @test test_on_grid(float_vf1,vfn2,res)
        @test test_on_grid(float_vf2,vfn2,res)
        @test test_on_grid(complex_vf,vfn2,res)
        @test test_on_grid(vfn1,vfn2,res)
        @test test_on_grid(vfn2,vfn2,res)
        # test interpolation result
        @test test_interpolate(float_vf1,func1) <= 1e-2
        @test test_interpolate(float_vf2,func2) <= 1e-2
        @test test_interpolate(complex_vf,func3)<= 1e-2

        #3D Scalar
        @info "3D Scalar"
        func1 = (x,y,z)->(sin(x/10.0)+cos(y/10.0)+sin(z/1.0))*exp(-((x-75.0)^2+(y-75.0)^2)/20^2)*exp(-(z-5.0)^2/(2.0^2))
        func2 = (x,y,z)->(sin(x/10.0)+cos(y/10.0)+sin(z/1.0))*exp(-((x-50.0)^2+(y-50.0)^2)/30^2)*exp(-(z-5.0)^2/(2.0^2))
        func3 = (x,y,z)->(sin(x/10.0)+cos(y/10.0)+sin(z/1.0))*exp(-((x-50.0)^2+(y-50.0)^2)/20^2)*exp(-(z-5.0)^2/(2.0^2))
        float_sf1 = func2field(ScalarField{Float64,3},func1,[0.05,0.05,1.0],[50.0,50.0,0.0],[50.0,50.0,10.0],name = "float_sf1")
        float_sf2 = func2field(ScalarField{Float64,3},func2,[1.0,1.0,1.0],[0.0,0.0,0.0],[100.0,100.0,10.0],name = "float_sf2")
        complex_sf = func2field(ScalarField{Complex{Float64},3},func3,[0.06,0.06,1.0],[20.0,20.0,0.0],[60.0,60.0,10.0],name = "complex_sf")
        sfn1 = ScalarFieldNode{3}([float_sf1,complex_sf],name = ascii("sfn1"))
        sfn2 = ScalarFieldNode{3}([float_sf2,sfn1],name = ascii("sfn2"))
        align_field_tree!(sfn2)
        # test if all points ar on grid
        @test float_sf1.res ≈ float_sf2.res ≈ complex_sf.res ≈ sfn1.res ≈ sfn2.res
        res = sfn2.res
        @test test_on_grid(float_sf1,sfn2,res)
        @test test_on_grid(float_sf2,sfn2,res)
        @test test_on_grid(complex_sf,sfn2,res)
        @test test_on_grid(sfn1,sfn2,res)
        @test test_on_grid(sfn2,sfn2,res)
        # test interpolation result
        @test test_interpolate(float_sf1,func1) <= 1e-2
        @test test_interpolate(float_sf2,func2) <= 1e-2
        @test test_interpolate(complex_sf,func3)<= 1e-2
        #3D Vector
        @info "3D Vector"
        func1 = (x,y,z)->repeat([(sin(x/10.0)+cos(y/10.0)+sin(z/1.0))*exp(-((x-75.0)^2+(y-75.0)^2)/20^2)*exp(-(z-5.0)^2/(2.0^2))],3)
        func2 = (x,y,z)->repeat([(sin(x/10.0)+cos(y/10.0)+sin(z/1.0))*exp(-((x-50.0)^2+(y-50.0)^2)/30^2)*exp(-(z-5.0)^2/(2.0^2))],3)
        func3 = (x,y,z)->repeat([(sin(x/10.0)+cos(y/10.0)+sin(z/1.0))*exp(-((x-50.0)^2+(y-50.0)^2)/20^2)*exp(-(z-5.0)^2/(2.0^2))],3)
        float_vf1 = func2field(VectorField{Float64,3},func1,[0.05,0.05,1.0],[50.0,50.0,0.0],[50.0,50.0,10.0],name = "float_vf1")
        float_vf2 = func2field(VectorField{Float64,3},func2,[1.0,1.0,1.0],[0.0,0.0,0.0],[100.0,100.0,10.0],name = "float_vf2")
        complex_vf = func2field(VectorField{Complex{Float64},3},func3,[0.06,0.06,1.0],[19.9,19.9,0.0],[60.0,60.0,10.0],name = "complex_vf")
        vfn1 = VectorFieldNode{3}([float_vf1,complex_vf],name = ascii("vfn1"))
        vfn2 = VectorFieldNode{3}([float_vf2,vfn1],name = ascii("vfn2"))
        align_field_tree!(vfn2)
        # test if all points ar on grid
        @test float_vf1.res ≈ float_vf2.res ≈ complex_vf.res ≈ vfn1.res ≈ vfn2.res
        res = vfn2.res
        @test test_on_grid(float_vf1,vfn2,res)
        @test test_on_grid(float_vf2,vfn2,res)
        @test test_on_grid(complex_vf,vfn2,res)
        @test test_on_grid(vfn1,vfn2,res)
        @test test_on_grid(vfn2,vfn2,res)
        # test interpolation result
        @test test_interpolate(float_vf1,func1) <= 1e-2
        @test test_interpolate(float_vf2,func2) <= 1e-2
        @test test_interpolate(complex_vf,func3)<= 1e-2
    end
end
