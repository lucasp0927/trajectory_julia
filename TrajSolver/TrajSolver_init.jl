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
#    println("initialize TrajSolver module on process ", myid())
    global my_trajnum
    global solver, reltol, abstol
    global trajnum, tspan, tdiv
    global temperature, init_speed, init_range
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
    temperature = float(config["atom-config"]["temperature"])::Float64
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
#    boundary = reduce(hcat,map(d->convert(Array{Float64},d),values(config["boundary"])))
    #calculate my_trajnum
    jobs = allocate_jobs(trajnum)
    my_trajnum = jobs[2]-jobs[1]+1
    result = Array(Float64,4,length(tspan),my_trajnum)
end

function distribute_atoms(init_range::Vector{Float64},atom_temp::Float64,t::Float64)
#    println("initialize atoms...")
    srand()
    x_range = init_range[1:2]
    y_range = init_range[3:4]
    U_range = Fields.composite_slow(init_range,t)
    U_min = minimum(U_range)
    U_max = maximum(U_range)
    U_min -= (U_max-U_min)/1000
#    println("U_min= ",U_min)
#    println("U_max= ",U_max)
    init_xv = zeros(Float64,(4,my_trajnum))
    for i = 1:my_trajnum::Int64
        while true
            #randomize position
            x = (x_range[2]-x_range[1])*rand()+x_range[1]
            y = (y_range[2]-y_range[1])*rand()+y_range[1]
            eu = Fields.value3([x,y],t)
            #randomize speed
            vp = sqrt(2.0*KB*atom_temp/M_CS)
            vx = -4.0*vp+8.0*vp*rand()
            vy = -4.0*vp+8.0*vp*rand()
            vz = -4.0*vp+8.0*vp*rand()
            ek = 0.5*M_CS*(vx*vx+vy*vy+vz*vz)/KB
            etot = ek+eu-U_min
            p = exp(-1.0*etot/atom_temp)
            if p<0.0 || p > 1.0
                warn("p=$p: out of range!")
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
