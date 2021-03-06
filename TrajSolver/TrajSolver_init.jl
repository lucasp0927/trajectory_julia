function allocate_jobs(totaljob)
    nchunks = nworkers()
    jobs = fill(fld(totaljob,nchunks),nchunks)
    jobs[1:totaljob%nchunks] .+= 1
    start = round(Int64,sum(jobs[1:myid()-2])+1)
    stop = round(Int64,sum(jobs[1:myid()-1]))
    return start,stop
end

function init_parallel(config::Dict)
#    println("start initialization TrajSolver module...")
    @sync begin
        for p = 1:Distributed.nprocs()
            @async Distributed.remotecall_wait(init!,p,config)
        end
    end
end

struct PeriodicCondition
    periodic_condition::Bool
    dim::Int64
    periodic_start::Float64
    periodic_end::Float64
end

struct BoundaryCondition
    boundary_condition::Bool
    xmin::Float64
    ymin::Float64
    zmin::Float64
    xmax::Float64
    ymax::Float64
    zmax::Float64
end

function init!(config::Dict)
    #    srand()
    Random.seed!()
    #    println("initialize TrajSolver module on process ", myid())
    global sim_type
    global my_trajnum
    global solver, reltol, abstol
    global trajnum, tspan, tdiv
    global radial_temperature, axial_temperature, init_speed, init_range
    global in_boundaries, out_boundaries
    #global result
    global trajsolver_config
    global periodic_condition
    global boundary_condition
    trajsolver_config = config
    #simulation-config
    sim_type = config["simulation-type"]::String
    @assert (sim_type == "2D") | (sim_type == "3D")
    trajnum = round(Int64,config["simulation-config"]["traj_num"])::Int64
    tstart = float(config["simulation-config"]["tstart"])::Float64
    tend = float(config["simulation-config"]["tend"])::Float64
    tdiv = float(config["simulation-config"]["tdiv"])::Float64
    tspan = collect(Float64,tstart:tdiv:tend)
    #solver-config
    solver = ascii(config["solver-config"]["solver"])
    reltol = float(config["solver-config"]["reltol"])::Float64
    abstol = float(config["solver-config"]["abstol"])::Float64
    #atom-config
    if trajsolver_config["atom-config"]["init-type"] == "fit-trap"
        radial_temperature = float(config["atom-config"]["radial-temperature"])::Float64
        axial_temperature = float(config["atom-config"]["axial-temperature"])::Float64
        init_speed = float(config["atom-config"]["init-speed"])::Float64
        init_range = convert(Dict{String,Vector{Float64}},config["atom-config"]["init-range"])
    end
    #boundary
    if sim_type == "2D"
        in_boundaries = map(values(config["in-boundary"]))do x
            return Polygon([promote(x...)...])
        end
        in_boundaries = [promote(in_boundaries...)...]::Vector{Polygon}
        out_boundaries = map(values(config["out-boundary"]))do x
            return Polygon([promote(x...)...])
        end
        out_boundaries = [promote(out_boundaries...)...]::Vector{Polygon}
    end
    #boundary "3D"
    if "boundary-condition" in keys(trajsolver_config)
        boundary_condition = BoundaryCondition(true,
                                               Float64(trajsolver_config["boundary-condition"]["xmin"]),
                                               Float64(trajsolver_config["boundary-condition"]["ymin"]),
                                               Float64(trajsolver_config["boundary-condition"]["zmin"]),
                                               Float64(trajsolver_config["boundary-condition"]["xmax"]),
                                               Float64(trajsolver_config["boundary-condition"]["ymax"]),
                                               Float64(trajsolver_config["boundary-condition"]["zmax"]))
    else
        boundary_condition = BoundaryCondition(false,0.0,0.0,0.0,0.0,0.0,0.0)
    end
    #periodic condition
    if "periodic-condition" in keys(trajsolver_config)
        periodic_condition = PeriodicCondition(true,
                                               Int64(trajsolver_config["periodic-condition"]["dim"]),
                                               Float64(trajsolver_config["periodic-condition"]["start"]),
                                               Float64(trajsolver_config["periodic-condition"]["end"]))

    else
        periodic_condition = PeriodicCondition(false,0,0.0,0.0)
    end

    #calculate my_trajnum
    jobs = allocate_jobs(trajnum)
    my_trajnum = jobs[2]-jobs[1]+1
    #result = Array{Float64}(undef,4,length(tspan),my_trajnum)
end

function range2sfn(range)
    @info "in range2sfn range: $(range)"
    #res = [(range[2]-range[1])/10000,(range[4]-range[3])/10000]
    #@info "in range2sfn res: $(res)"
    #init_U_data = Fields.composite_with_position(range,tspan[1],res)
    #TODO: init with all fields_arr
    res = Fields.fields_arr[1].res
    @info "res: $(res)"
    @info "range: $(range)"
    init_U_data = Fields.composite_with_position(range,tspan[1],res)
    prob,xstart,xend,ystart,yend = fit_trap(init_U_data,axial_temperature,radial_temperature)
    prob_s = copy_to_sharedarray!(prob)
    prob_f = ScalarFieldNode{2}([ScalarField{Float64,2}(prob_s,[xstart,ystart],[xend-xstart,yend-ystart])])
    Fields.set_geometry!(prob_f)
    Fields.set_typeof!(prob_f)
#    try
#       output_image_gp(tspan[1],range,"init_prob_"*string(iter)*string(range)*".png",prob_f)
#    end
    return prob_f
end

function prepare_U_prob()
    #calculate probablility distrubution, and construct a scalar field object.
    #calculate position probability
    sfns = map(range2sfn,values(init_range))
    sfns = [promote(sfns...)...]
    @sync begin
        for p = 1:nprocs()
            map(Fields.clean_scaling!,sfns)
            @async remotecall_fetch(init_U_prob!,p,sfns)
        end
    end
    map(Fields.eval_scaling!,sfns)
end

function init_U_prob!(sfns::Vector{ScalarFieldNode{2}})
    global U_prob
    U_prob = 0
    GC.gc()
    U_prob = map(Fields.copyfield,sfns)
    GC.gc()
end

function distribute_atoms_one_shot()
    pancake_id = rand(1:length(U_prob))
    result = distribute_atoms_inner(U_prob[pancake_id],1)
    init_xv = dropdims(result,dims=2)
    return init_xv
end

function distribute_atoms()
    #distribute my_trajnum atoms among pancakes.
    # pancake_num = length(U_prob)
    # traj_num = zeros(Int64,pancake_num)
    # d,r = divrem(my_trajnum,pancake_num)
    # traj_num[:] = d
    # traj_num[end] += r
    # @assert sum(traj_num) == my_trajnum
    presult = @showprogress pmap(i->distribute_atoms_one_shot(),1:trajnum)
    init_xv = cat(presult...,dims=2)
    # @info "presult: $(presult)"
    # @info "init_xv size: $(size(init_xv))"
    # @info "init_xv: $(init_xv)"
    return init_xv
#    init_xvs = map(i->distribute_atoms_inner(U_prob[i],traj_num[i]),1:pancake_num)
#    init_xv = cat(2,init_xvs...)
#    return init_xv
end

function distribute_atoms_inner(sfn::ScalarFieldNode{2},traj_num::Int64)
    t = tspan[1]
    x_range = [sfn.position[1],sfn.position[1]+sfn.size[1]]
    y_range = [sfn.position[2],sfn.position[2]+sfn.size[2]]
    @debug "distrubute atoms at t=$t, t_axial=$axial_temperature, t_radial=$radial_temperature, x_range=$x_range, y_range=$y_range,  traj_num=$traj_num"
    init_xv = zeros(Float64,(4,traj_num))
    vp_a = sqrt(2.0*KB*axial_temperature/M_CS)
    vp_r = sqrt(2.0*KB*radial_temperature/M_CS)
    for i = 1:traj_num
        while true
            x = (x_range[2]-x_range[1])*rand()+x_range[1]
            y = (y_range[2]-y_range[1])*rand()+y_range[1]
            # @info "x: $(x) y: $(y)"
            p_pos = Fields.value([x,y],t,sfn)

            vx = -4.0*vp_a+8.0*vp_a*rand()
            vy = -4.0*vp_r+8.0*vp_r*rand()
            vz = -4.0*vp_r+8.0*vp_r*rand()

            ek_a = 0.5*M_CS*(vx^2)/KB
            ek_r = 0.5*M_CS*(vy^2+vz^2)/KB
            p_vel = exp(-1.0*ek_a/axial_temperature)*exp(-1.0*ek_r/radial_temperature)
            p = p_pos*p_vel
            if p<0.0 || p > 1.0
                if abs(p)<1e-5
                    p = 0.0
                # else
                #     @warn "p=$p: out of range!"
                end
            end
            init_xv[1,i] = x
            init_xv[2,i] = y
            init_xv[3,i] = vx+init_speed
            init_xv[4,i] = vy
            if rand() <= p
                break
            end
        end
    end
    return init_xv
end
