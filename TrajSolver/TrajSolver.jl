module TrajSolver
using Sundials
using Fields
include("../constant.jl")
include("polygon.jl")
global my_trajnum
global solver, reltol, abstol
global trajnum, tspan, tdiv
global temperature, init_speed, init_range
global in_boundaries
global out_boundaries
global result

function get_tspan()
    return tspan
end

function get_result()
    return result
end

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
        poly = convert(Array{Float64,2},reshape(x,(round(Int64,length(x)/2)),2))'
        return Polygon(poly)
    end
    in_boundaries = [promote(in_boundaries...)...]::Vector{Polygon}
    out_boundaries = map(values(config["out-boundary"]))do x
        poly = convert(Array{Float64,2},reshape(x,(round(Int64,length(x)/2)),2))'
        return Polygon(poly)
    end
    out_boundaries = [promote(out_boundaries...)...]::Vector{Polygon}
#    boundary = reduce(hcat,map(d->convert(Array{Float64},d),values(config["boundary"])))
    #calculate my_trajnum
    jobs = allocate_jobs(trajnum)
    my_trajnum = jobs[2]-jobs[1]+1
    result = Array(Float64,4,length(tspan),my_trajnum)
end

@inbounds function boundary(pos::Vector{Float64})
    #return false if particle should be removed
    for p in in_boundaries::Vector{Polygon}
        if pointInPolygon(p,pos)==true
            return false
        end
    end
    for p in out_boundaries::Vector{Polygon}
        if pointInPolygon(p,pos)==false
            return false
        end
    end
    return true
end

function mycvode(mem, f::Function, y0::Vector{Float64}, t::Vector{Float64} , yout::SubArray; reltol::Float64=1e-8, abstol::Float64=1e-7)
    # f, Function to be optimized of the form f(y::Vector{Float64}, fy::Vector{Float64}, t::Float64)
    #    where `y` is the input vector, and `fy` is the
    # y0, Vector of initial values
    # t, Vector of time values at which to record integration results
    # reltol, Relative Tolerance to be used (default=1e-4)
    # abstol, Absolute Tolerance to be used (default=1e-6)
    flag = Sundials.CVodeInit(mem, cfunction(Sundials.cvodefun, Int32, (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Ref{Function})), t[1], Sundials.nvector(y0).ptr[1])
    flag = Sundials.CVodeSetUserData(mem, f)
    flag = Sundials.CVodeSStolerances(mem, reltol, abstol)
    flag = Sundials.CVDense(mem, length(y0))
    yout[1:2,1] = y0[1:2]
#    yout[3,1] = Fields.value3(y0[1:2],t[1])
    y = copy(y0)
    tout = [t[1]]
    for k in 2:length(t)
        flag = Sundials.CVode(mem, t[k], y, tout, Sundials.CV_NORMAL)
        yout[1:2,k] = y[1:2]
#        yout[3,k] = Fields.value3(y[1:2],t[k])
#        yout[4,k] = t[k]
        if boundary(yout[1:2,k]) == false
            break
        end
    end
end

function solve_traj()
    fill!(result,NaN)
    init_xv = distribute_atoms(init_range,temperature,tspan[1])
    #initialize sundials
    if solver == "ADAMS"
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_ADAMS, Sundials.CV_FUNCTIONAL))
    elseif solver == "BDF"
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_BDF, Sundials.CV_NEWTON))
    else
        error("Unknown ODE solver $solver.")
    end
    init = zeros(Float64,4)
    for i = 1:my_trajnum::Int64
        init[:] = init_xv[:,i]
        yout = slice(result::Array{Float64,3},:,:,i)
        mycvode(mem,Fields.gradient!,init,tspan,yout; reltol = reltol, abstol =abstol)
    end
    gc()
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
end
