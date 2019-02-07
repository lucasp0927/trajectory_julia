module TrajSolver
using Distributed
using SharedArrays
using Statistics
using Random
using Printf
using DifferentialEquations
using Sundials
using Fields
using ProgressMeter

export calculate_traj
global sim_type
global my_trajnum
global solver, reltol, abstol
global trajnum, tspan, tdiv
global radial_temperature, axial_temperature, init_speed, init_range
global in_boundaries
global out_boundaries
#global result
global U_prob #probability for atom distribution
global trajsolver_config
global periodic_condition
global boundary_condition
include("../constant.jl")
include("../fileio.jl")
include("polygon.jl")
include("TrajSolver_init.jl")
include("TrajSolver_get.jl")
include("TrajSolver_conditions.jl")
include("fit_trap.jl")
include("../TrajAnalyzer/TrajAnalyzer_output.jl")
include("../TrajAnalyzer/TrajAnalyzer_trajectories.jl")
function calculate_traj()
    @info "distribute atoms..."
    #init trajnum initial xv pairs
    if trajsolver_config["atom-config"]["init-type"] == "fit-trap"
        @info "fitting trap..."
        prepare_U_prob()
        init_xv = distribute_atoms()
    elseif trajsolver_config["atom-config"]["init-type"] == "from-file"
        @info "reading initial condition from "*trajsolver_config["atom-config"]["filename"]*"..."
        result = h5todict(trajsolver_config["atom-config"]["filename"])
        result["tspan"] = vec(result["tspan"])
        #result = matread(trajsolver_config["atom-config"]["filename"])
        traj_s = copy_to_sharedarray!(result["traj"])
        trajectories = Trajectories(result,traj_s)
        t = float(trajsolver_config["atom-config"]["time"])
        @assert t==tspan[1]
        init_xv = trajectories[t,:]
        
        #cleanup shared array
        foreach(traj_s.refs) do r
            @spawnat r.where finalize(fetch(r))
        end
        finalize(traj_s.s)
        finalize(traj_s)        
        foreach(trajectories.traj.refs) do r
            @spawnat r.where finalize(fetch(r))
        end
        finalize(trajectories.traj.s)
        finalize(trajectories.traj)
        @everywhere GC.gc()        
    else
        err("Unknown init-type in atom-config.")
    end
    #preallocate result according to dim
    if sim_type == "2D"
        traj = Array{Float64}(undef,4,length(tspan),trajnum)
        @info "calculate trajectories..."
        @time calculate_traj_inner_parallel_2d(init_xv,traj)
    else
        traj = Array{Float64}(undef,6,length(tspan),trajnum)
        @info "calculate trajectories..."
        @time calculate_traj_inner_parallel_3d(init_xv,traj)
    end
    #truncate save time range
    t_idx_start = searchsortedlast(tspan,trajsolver_config["save-range"]["tstart"])
    t_idx_end = searchsortedfirst(tspan,trajsolver_config["save-range"]["tend"])
    traj_save = traj[:,t_idx_start:t_idx_end,:]
    tspan_save = tspan[t_idx_start:t_idx_end]
    traj=nothing
    result = Dict(
                  "traj"=>traj_save,
                  "tspan"=>tspan_save,
                  "pos"=>Fields.fields.position,
                  "siz"=>Fields.fields.size,
                  "reltol"=>reltol,
                  "abstol"=>abstol
    )
    @everywhere GC.gc()
    return result
end

function calculate_traj_inner_parallel_2d(init_xv,traj)
    i = 1
    pm = Progress(trajnum,1)
    @sync begin
        for p = 2:nprocs()
            @async begin
                while true
                    next!(pm)
                    idx=i
                    i+=1
                    if idx>trajnum
                        break
                    end
                    traj[:,:,idx] = remotecall_fetch(solve_traj_one_shot_2d,p,init_xv[:,idx])
                end
            end
        end
    end
end

function calculate_traj_inner_parallel_3d(init_xv,traj)
    i = 1
    pm = Progress(trajnum,1)
    @sync begin
        for p = 2:nprocs()
            @async begin
                while true
                    next!(pm)
                    idx=i
                    i+=1
                    if idx>trajnum
                        break
                    end
                    traj[:,:,idx] = remotecall_fetch(solve_traj_one_shot_3d,p,init_xv[:,idx])
                end
            end
        end
    end
end

# function calculate_traj_inner(init_xv,traj)
#     i = 1
#     pm = Progress(trajnum,1)
#     for i = 1:trajnum
#         next!(pm)
#         if i>trajnum
#             break
#         end
#         traj[:,:,i] = solve_traj_one_shot(init_xv[:,i])
#     end
# end

function solve_traj_one_shot_2d(init_xv::Vector{Float64})
    yout = Array{Float64}(undef,4,length(tspan))
    fill!(yout,NaN)
    if any(isnan.(init_xv)) == false
        solve_eq_of_motion_2d(Fields.gradient_odejl!,init_xv,tspan,yout;reltol=reltol, abstol=abstol)
    end
    return yout
end

function solve_traj_one_shot_3d(init_xv::Vector{Float64})
    yout = Array{Float64}(undef,6,length(tspan))
    fill!(yout,NaN)
    if any(isnan.(init_xv)) == false
        solve_eq_of_motion_3d(Fields.gradient_odejl!,init_xv,tspan,yout;reltol=reltol, abstol=abstol)
    end
    return yout
end

function solve_eq_of_motion_2d(f::Function, y0::Vector{Float64}, t::Vector{Float64} , yout::T; reltol::Float64=1e-8, abstol::Float64=1e-7, mxstep::Int64=Integer(1e6)) where T<:AbstractArray
    #TODO: add more solver options
    if solver == "ADAMS"
        solver_alg = CVODE_Adams()
    elseif solver == "BDF"
        solver_alg = CVODE_BDF()
    end
    prob = ODEProblem(f,y0,(t[1],t[end]))
    integrator = init(prob, solver_alg; abstol=abstol,reltol=reltol)
    yout[:,1] = y0
    for i = 2:length(t)
        dt = t[i] - t[i-1]
        step!(integrator,dt,true)
        yout[:,i] = integrator.u
        if boundary_2d(yout[1:2,i]) == false
            break
        end
    end
end

function solve_eq_of_motion_3d(f::Function, y0::Vector{Float64}, t::Vector{Float64} , yout::T; reltol::Float64=1e-8, abstol::Float64=1e-7, mxstep::Int64=Integer(1e6)) where T<:AbstractArray
    #TODO: add more solver options
    if solver == "ADAMS"
        solver_alg = CVODE_Adams()
    elseif solver == "BDF"
        solver_alg = CVODE_BDF()
    end
    #remove atom when atom is in dielectric
    #t_diff = mean(diff(t));
    #cb_material = PeriodicCallback(check_material_boundary!,t_diff) #not working
    cb_material = DiscreteCallback(condition_material,affect_material!)
    cb = cb_material    
    #periodic boundary condition with callback
    prob = ODEProblem(f,y0,(t[1],t[end]))
    if (periodic_condition::PeriodicCondition).periodic_condition::Bool
        #cb1 = ContinuousCallback(condition1,affect1!)
        #cb2 = ContinuousCallback(condition2,affect2!)
        cb1 = DiscreteCallback(condition1_fast,affect1_fast!)
        cb2 = DiscreteCallback(condition2_fast,affect2_fast!)
        cb = CallbackSet(cb1,cb2,cb)
        #integrator = init(prob, solver_alg; abstol=abstol,reltol=reltol,callback=cb)
    end

    if (boundary_condition::BoundaryCondition).boundary_condition::Bool
        cb_boundary = DiscreteCallback(condition_boundary,affect_boundary!)
        cb = CallbackSet(cb_boundary,cb)
    end
    
    sol = solve(prob, solver_alg; abstol=abstol,reltol=reltol,callback=cb,saveat=t)
    #copy results to yout
    j = 1
    for i = 1:length(t)
        while j<= length(sol.t) && sol.t[j] < t[i]
            j = j+1
        end
        if j>length(sol.t)
            break
        else
            @assert isapprox(sol.t[j],t[i])
            yout[:,i] = sol.u[j]
        end
    end
end

@inbounds function boundary_2d(pos::Vector{Float64})
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

@inbounds function boundary_3d(pos::Vector{Float64},t::Float64)
    if Fields.in_field(Fields.material,pos)
        if Fields.value(pos,t,Fields.material) > 0.5
            return false
        else
            return true
        end
    else
        return true
    end
end

end
