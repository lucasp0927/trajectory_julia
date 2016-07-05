function allocate_jobs(totaljob)
    nchunks = nworkers()
    jobs = fill(fld(totaljob,nchunks),nchunks)
    jobs[1:totaljob%nchunks] += 1
    start = round(Int64,sum(jobs[1:myid()-2])+1)
    stop = round(Int64,sum(jobs[1:myid()-1]))
    return start,stop
end

function init_parallel(config::Dict)
#    println("start initialization TrajSolver module...")
    @sync begin
        for p = 1:nprocs()
            @async remotecall_wait(p,init!,config)
        end
    end
end

function init!(config::Dict)
    srand()
#    println("initialize TrajSolver module on process ", myid())
    global my_trajnum
    global solver, reltol, abstol
    global trajnum, tspan, tdiv
    global radial_temperature, axial_temperature, init_speed, init_range
    global in_boundaries, out_boundaries
    global result
    #simulation-config
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
    radial_temperature = float(config["atom-config"]["radial-temperature"])::Float64
    axial_temperature = float(config["atom-config"]["axial-temperature"])::Float64
    init_speed = float(config["atom-config"]["init-speed"])::Float64
    init_range = convert(Vector{Float64},config["atom-config"]["init-range"])
    #boundary
    in_boundaries = map(values(config["in-boundary"]))do x
        return Polygon([promote(x...)...])
    end
    in_boundaries = [promote(in_boundaries...)...]::Vector{Polygon}
    out_boundaries = map(values(config["out-boundary"]))do x
        return Polygon([promote(x...)...])
    end
    out_boundaries = [promote(out_boundaries...)...]::Vector{Polygon}
    #calculate my_trajnum
    jobs = allocate_jobs(trajnum)
    my_trajnum = jobs[2]-jobs[1]+1
    result = Array(Float64,4,length(tspan),my_trajnum)
end

function distribute_atoms()
    t = tspan[1]
    Lumberjack.debug("distrubute atoms at t=$t, t_axial=$axial_temperature, t_radial=$radial_temperature, range=$init_range")
    x_range = init_range[1:2]
    y_range = init_range[3:4]
    init_xv = zeros(Float64,(4,my_trajnum))
    vp_a = sqrt(2.0*KB*axial_temperature/M_CS)
    vp_r = sqrt(2.0*KB*radial_temperature/M_CS)
    for i = 1:my_trajnum::Int64
        while true
            x = (x_range[2]-x_range[1])*rand()+x_range[1]
            y = (y_range[2]-y_range[1])*rand()+y_range[1]
            p_pos = Fields.value3([x,y],t,U_prob)
            vx = -4.0*vp_a+8.0*vp_a*rand()
            vy = -4.0*vp_r+8.0*vp_r*rand()
            vz = -4.0*vp_r+8.0*vp_r*rand()
            ek_a = 0.5*M_CS*(vx^2)/KB
            ek_r = 0.5*M_CS*(vy^2+vz^2)/KB
            p_vel = exp(-1.0*ek_a/axial_temperature)*exp(-1.0*ek_r/radial_temperature)
            p = p_pos*p_vel
            if p<0.0 || p > 1.0
                if abs(p)<1e-10
                    p = 0.0
                else
                    Lumberjack.warn("p=$p: out of range!")
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
