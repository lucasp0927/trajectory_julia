module TrajSolver
using Sundials
using Fields
using Logging
#using Optim
using ProgressMeter

export calculate_traj
global my_trajnum
global solver, reltol, abstol
global trajnum, tspan, tdiv
global radial_temperature, axial_temperature, init_speed, init_range
global in_boundaries
global out_boundaries
global result
global U_prob #probability for atom distribution
global trajsolver_config

include("../constant.jl")
include("../fileio.jl")
include("polygon.jl")
include("TrajSolver_init.jl")
include("TrajSolver_get.jl")
include("../TrajAnalyzer/TrajAnalyzer_output.jl")
include("fit_trap.jl")

function calculate_traj()
    info("calculate trajectories...")
#    parallel_set_iter(i)
    prepare_U_prob()
    traj = Array(Float64,4,length(tspan),trajnum)
    # function to produce the next work item from the queue.
    # in this case it's just an index.
    i = 1
    pm = Progress(trajnum, 1)
    nextidx() = (next!(pm);idx=i; i+=1; idx)
#=
    Profile.init(delay=0.01)
    while true
        idx = nextidx()
        if idx>trajnum
            break
        end
        traj[:,:,idx] = solve_traj_one_shot()
    end
    Profile.clear()
    i = 1
    pm = Progress(trajnum, 1)
    @profile begin
        while true
            idx = nextidx()
            if idx>trajnum
                break
            end
            traj[:,:,idx] = solve_traj_one_shot()
        end
    end
    Profile.print()
    ProfileView.view()
    sleep(10000)
=#

    @time @sync begin
        for p = 2:nprocs()
            @async begin
                while true
                    idx = nextidx()
                    if idx>trajnum
                        break
                    end
                    traj[:,:,idx] = remotecall_fetch(solve_traj_one_shot,p)
                end
            end
        end
    end
    #truncate save time range
    t_idx_start = searchsortedlast(tspan,trajsolver_config["save-range"]["tstart"])
    t_idx_end = searchsortedfirst(tspan,trajsolver_config["save-range"]["tend"])
    traj_save = traj[:,t_idx_start:t_idx_end,:]
    tspan_save = tspan[t_idx_start:t_idx_end]
    # println("t_idx_start = ",t_idx_start)
    # println("t_idx_end = ",t_idx_end)
    # println("tspan[t_idx_start] = ", tspan[t_idx_start])
    # println("tspan[t_idx_end] = ", tspan[t_idx_end])
    result = Dict(
                  "traj"=>traj_save,
                  "tspan"=>tspan_save,
                  "pos"=>Fields.fields.position,
                  "siz"=>Fields.fields.size,
                  "radial_temperature"=>radial_temperature,
                  "axial_temperature"=>axial_temperature,
                  "init_speed"=>init_speed,
                  "reltol"=>reltol,
                  "abstol"=>abstol
                  )
    gc()
    return result
end

function calculate_traj_unbalanced()
    info("calculate trajectories...")
#    parallel_set_iter(i)
    prepare_U_prob()
    @time @sync begin
        for p = 2:nprocs()
            @async remotecall_wait(solve_traj,p)
        end
    end
    temp = cell(nworkers())
    for p = 2:nprocs()
        temp[p-1] = remotecall_fetch(get_result,p)
    end
    traj = cat(3,temp...)
    result = Dict(
                  "traj"=>traj,
                  "tspan"=>tspan,
                  "pos"=>Fields.fields.position,
                  "siz"=>Fields.fields.size,
                  "radial_temperature"=>radial_temperature,
                  "axial_temperature"=>axial_temperature,
                  "init_speed"=>init_speed,
                  "reltol"=>reltol,
                  "abstol"=>abstol
                  )
    return result
end

function solve_traj_one_shot()
    init_xv = distribute_atoms_one_shot()
    yout = Array(Float64,4,length(tspan))
    fill!(yout,NaN)
    mycvode(Fields.gradient!,init_xv,tspan,yout;reltol=reltol, abstol=abstol)
    return yout
end

function distribute_atoms_one_shot()
    pancake_id = rand(1:length(U_prob))
    init_xv = squeeze(distribute_atoms_inner(U_prob[pancake_id],1),2)
    return init_xv
end

#=
function solve_traj()
    fill!(result,NaN)
    init_xv = distribute_atoms()
    #initialize sundials
    if solver == "ADAMS"
        debug("Using solver ADAMS")
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_ADAMS, Sundials.CV_FUNCTIONAL))
    elseif solver == "BDF"
        debug("Using solver BDF")
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_BDF, Sundials.CV_NEWTON))
    else
        err("Unknown ODE solver $solver.")
    end
    init = zeros(Float64,4)
    for i = 1:my_trajnum::Int64
        init[:] = init_xv[:,i]
        yout = slice(result::Array{Float64,3},:,:,i)
        mycvode(mem,Fields.gradient!,init,tspan,yout; reltol = reltol, abstol =abstol)
    end
    #TODO: figure out how to delete mem
    gc()
end
=#
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

@inbounds function mycvode{T<:AbstractArray}(f::Function, y0::Vector{Float64}, t::Vector{Float64} , yout::T; reltol::Float64=1e-8, abstol::Float64=1e-7, mxstep::Int64=Integer(1e6), userdata::Any=nothing)
    # f, Function to be optimized of the form f(y::Vector{Float64}, fy::Vector{Float64}, t::Float64)
    #    where `y` is the input vector, and `fy` is the
    # y0, Vector of initial values
    # t, Vector of time values at which to record integration results
    # reltol, Relative Tolerance to be used (default=1e-4)
    # abstol, Absolute Tolerance to be used (default=1e-6)
    if solver == "BDF"
        mem = Sundials.CVodeCreate(Sundials.CV_BDF, Sundials.CV_NEWTON)
    elseif solver == "ADAMS"
        mem = Sundials.CVodeCreate(Sundials.CV_ADAMS, Sundials.CV_FUNCTIONAL)
    end
    if mem == C_NULL
        err("Failed to allocate CVODE solver object")
    end

    try
        userfun = Sundials.UserFunctionAndData(f, userdata)
        y0nv = Sundials.NVector(y0)
        #flag = Sundials.CVodeInit(mem, cfunction(Sundials.cvodefun, Int32, (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Ref{Function})), t[1], Sundials.nvector(y0).ptr[1])
        flag = Sundials.CVodeInit(mem, cfunction(Sundials.cvodefun, Cint, (Sundials.realtype, Sundials.N_Vector, Sundials.N_Vector, Ref{typeof(userfun)})), t[1], convert(Sundials.N_Vector, y0nv))
        flag = Sundials.CVodeSetUserData(mem, userfun)
        flag = Sundials.CVodeSStolerances(mem, reltol, abstol)
        flag = Sundials.CVDense(mem, length(y0))
#        flag = Sundials.CVodeSetUserData(mem, f)
#        flag = Sundials.CVodeSStolerances(mem, reltol, abstol)
#        flag = Sundials.CVDense(mem, length(y0))
#        flag = Sundials.CVodeSetMaxNumSteps(mem,mxstep)
        yout[1:2,1] = y0[1:2]
        #yout[3,1] = Fields.value(y0[1:2],t[1])
        yout[3:4,1] = y0[3:4]
        ynv = Sundials.NVector(copy(y0))
        tout = [t[1]]
        for k in 2:length(t)
            flag = Sundials.CVode(mem, t[k], ynv, tout, Sundials.CV_NORMAL)
            #yout[3:4] are free spaces.
            #yout[1:2,k] = y[1:2]
            #yout[3:4,k] = y[3:4]
            yout[:,k] = convert(Vector,ynv)
            #        yout[3,k] = Fields.value(y[1:2],t[k])
            #        yout[4,k] = t[k]

            if boundary(yout[1:2,k]) == false
                break
            end
        end
    finally
        Sundials.CVodeFree(Ref{Sundials.CVODEMemPtr}(mem))
    end
end

end
