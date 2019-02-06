
function condition1(u,t,integrator) # Event when event_f(u,t) == 0
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_start = (periodic_condition::PeriodicCondition).periodic_start::Float64
    u[dim] - p_start
end

function affect1!(integrator)
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_end = (periodic_condition::PeriodicCondition).periodic_end::Float64
    integrator.u[dim] = p_end
end

function condition1_fast(u,t,integrator) # Event when event_f(u,t) == 0
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_start = (periodic_condition::PeriodicCondition).periodic_start::Float64
    u[dim] < p_start
end

function affect1_fast!(integrator)
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_start = (periodic_condition::PeriodicCondition).periodic_start::Float64
    p_end = (periodic_condition::PeriodicCondition).periodic_end::Float64
    dist = p_end - p_start
    integrator.u[dim] += dist
end

function condition2(u,t,integrator) # Event when event_f(u,t) == 0
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_end = (periodic_condition::PeriodicCondition).periodic_end::Float64
    u[dim] - p_end
end

function affect2!(integrator)
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_start = (periodic_condition::PeriodicCondition).periodic_start::Float64
    integrator.u[dim] = p_start
end

function condition2_fast(u,t,integrator) # Event when event_f(u,t) == 0
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_end = (periodic_condition::PeriodicCondition).periodic_end::Float64
    u[dim] > p_end
end

function affect2_fast!(integrator)
    dim = (periodic_condition::PeriodicCondition).dim::Int64
    p_start = (periodic_condition::PeriodicCondition).periodic_start::Float64
    p_end = (periodic_condition::PeriodicCondition).periodic_end::Float64
    dist = p_end - p_start    
    integrator.u[dim] -= dist
end


function condition_material(u,t,integrator)
    !boundary_3d(u[1:3],t)
end
affect_material!(integrator) = terminate!(integrator)

function condition_boundary(u,t,integrator)
    return (
        (u[2] < (boundary_condition::BoundaryCondition).ymin::Float64) ||
        (u[2] > (boundary_condition::BoundaryCondition).ymax::Float64) ||
        (u[1] < (boundary_condition::BoundaryCondition).xmin::Float64) ||
        (u[1] > (boundary_condition::BoundaryCondition).xmax::Float64) ||
        (u[3] < (boundary_condition::BoundaryCondition).zmin::Float64) ||
        (u[3] > (boundary_condition::BoundaryCondition).zmax::Float64)
    )
end
affect_boundary!(integrator) = terminate!(integrator)

function check_material_boundary!(integrator)
    if condition_boundary(integrator.u,integrator.t,integrator) || condition_material(integrator.u,integrator.t,integrator)
        terminate!(integrator)
    end
end
