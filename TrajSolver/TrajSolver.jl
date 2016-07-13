module TrajSolver
using Sundials
using Fields
using Lumberjack
using Optim
include("../constant.jl")
include("../fileio.jl")
include("polygon.jl")
include("TrajSolver_init.jl")
include("TrajSolver_get.jl")
include("../TrajAnalyzer/TrajAnalyzer_output.jl")
include("fit_trap.jl")
export calculate_traj
global iter
global my_trajnum
global solver, reltol, abstol
global trajnum, tspan, tdiv
global radial_temperature, axial_temperature, init_speed, init_range
global in_boundaries
global out_boundaries
global result
global U_prob #probability for atom distribution

function parallel_set_iter(i::Int64)
    @sync begin
        for p = 1:nprocs()
            @async remotecall_fetch(p,set_iter,i)
        end
    end
end

function set_iter(i::Int64)
    global iter
    iter = i
end

function calculate_traj(i::Int64)
    Lumberjack.info("calculate trajectories...")
    parallel_set_iter(i)
    prepare_U_prob()
    @time @sync begin
        for p = 2:nprocs()
            @async remotecall_wait(p,solve_traj)
        end
    end
    temp = cell(nworkers())
    for p = 2:nprocs()
        temp[p-1] = remotecall_fetch(p,get_result)
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

function solve_traj()
    fill!(result,NaN)
    init_xv = distribute_atoms()
    #initialize sundials
    if solver == "ADAMS"
        Lumberjack.debug("Using solver ADAMS")
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_ADAMS, Sundials.CV_FUNCTIONAL))
    elseif solver == "BDF"
        Lumberjack.debug("Using solver BDF")
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_BDF, Sundials.CV_NEWTON))
    else
        Lumberjack.error("Unknown ODE solver $solver.")
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

@inbounds function mycvode(mem, f::Function, y0::Vector{Float64}, t::Vector{Float64} , yout::SubArray; reltol::Float64=1e-8, abstol::Float64=1e-7, mxstep::Int64=Integer(1e6))
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
    flag = Sundials.CVodeSetMaxNumSteps(mem,mxstep)
    yout[1:2,1] = y0[1:2]
    #yout[3,1] = Fields.value3(y0[1:2],t[1])
    yout[3:4,1] = y0[3:4]
    y = copy(y0)
    tout = [t[1]]
    for k in 2:length(t)
        flag = Sundials.CVode(mem, t[k], y, tout, Sundials.CV_NORMAL)
        #yout[3:4] are free spaces.
        yout[1:2,k] = y[1:2]
        yout[3:4,k] = y[3:4]
#        yout[3,k] = Fields.value3(y[1:2],t[k])
        #        yout[4,k] = t[k]
        ###### ugly periodic hack
        if yout[2,k]>=69800
            yout[4,k] = -yout[4,k]
        elseif yout[2,k]<=200
            yout[4,k] = -yout[4,k]            
        end
            
        ######
        if boundary(yout[1:2,k]) == false
             break
        end
    end
end

end
