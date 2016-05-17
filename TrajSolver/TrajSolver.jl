module TrajSolver
using Sundials
using Fields
include("constant.jl")

function pcalc_trajectory!(result::SharedArray{Float64}, traj_num::Int64,t_span::Vector{Float64}, init_xv::SharedArray{Float64,2}, args::Tuple, progress::SharedArray{Bool}, solver::Vector{ASCIIString}, tol::Vector{Float64})
    idx = indexpids(result)
    if idx == 0
        return
    end
    set_parameters(args...)
    nchunks = length(procs(result))
    jobs = fill(fld(traj_num,nchunks),nchunks)
    jobs[1:traj_num%nchunks] += 1
    start = Integer(sum(jobs[1:idx-1])+1)
    stop = Integer(sum(jobs[1:idx]))
        #mem = Sundials.CVodeCreate(Sundials.CV_ADAMS, Sundials.CV_FUNCTIONAL)
    if solver[2] == "ADAMS"
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_ADAMS, Sundials.CV_FUNCTIONAL))
    elseif solver[2] == "BDF"
        mem = convert(Sundials.CVODE_ptr,Sundials.CVodeHandle(Sundials.CV_BDF, Sundials.CV_NEWTON))
    else
        error("Unknown ODE solver $solver.")
    end
    for i in collect(start:stop)
        init = zeros(Float64,4)
        for j = 1:4
            init[j] = init_xv[i,j]
        end
        yout = slice(result,:,:,i)
        mycvode(mem,paccfunc2!,init,t_span,yout; reltol = tol[1], abstol = tol[2])
        progress[i] = true
    end
    # if solver =="sundials"
    #     Sundials.CVodeFree([mem])
    # end
    gc()
end

function pcalc_trajectory_wrap!(result::SharedArray{Float64,3}, traj_num::Int64, t_span::Vector{Float64}, init_xv::SharedArray{Float64,2}, args::Tuple, solver::Vector{ASCIIString}, tol::Vector{Float64})
    progress =  SharedArray(Bool,(traj_num),pids=procs())
    fill!(progress,false)
    @sync begin
        for p in procs(result)
            @async begin
                remotecall_wait(p, pcalc_trajectory!, result, traj_num, t_span, init_xv, args, progress, solver, tol)
            end
        end
#       prog = Progress(traj_num, 1)
       while sum(progress) < traj_num
#           update!(prog, sum(progress))
           sleep(0.1)
       end
#       println("")
    end
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
    yout[:,1] = y0[1:2]
    y = copy(y0)
    tout = [t[1]]
    for k in 2:length(t)
        flag = Sundials.CVode(mem, t[k], y, tout, Sundials.CV_NORMAL)
        yout[:,k] = y[1:2]
        # if outside(yout[:,k]) == true
        #     break
        # end
    end
end

function pinitialize!(init_xv::SharedArray{Float64,2}, x_grid, x_div, y_grid, y_div, U_t::SharedArray{Float64})
    idx = indexpids(init_xv)
    if idx == 0
        return
    end
    #srand()
    #srand(idx)#TODO remove this!!
    x_center = 31
    x_width = 10
    traj_num = size(init_xv,1)
    nchunks = length(procs(init_xv))
    jobs = fill(fld(traj_num,nchunks),nchunks)
    jobs[1:traj_num%nchunks] += 1
    start = Integer(sum(jobs[1:idx-1])+1)
    stop = Integer(sum(jobs[1:idx]))
    t_min = minimum(U_t[:,x_center-x_width:x_center+x_width,1])
#    U_t_itp = Interpolations.interpolate(U_t[:,:,1], BSpline(Quadratic(Flat())),OnGrid())
    for i in collect(start:stop)
        flag = false
        while flag == false
            # randomize position
            #x = round(Int64,x_center-x_width+rand()*(2*x_width))
            #y = round(Int64,110/y_div+ rand()*(y_grid-220/y_div))
            x = x_center-x_width+rand()*(2*x_width)
            y = 110/y_div+ rand()*(y_grid-220/y_div)
            x_floor = floor(Int64,x)
            y_floor = floor(Int64,y)
            #evaluate U for linear bilinear interpolation
            U11 = U_t[y_floor,x_floor,1]
            U12 = U_t[y_floor+1,x_floor,1]
            U21 = U_t[y_floor,x_floor+1,1]
            U22 = U_t[y_floor+1,x_floor+1,1]
            # randomize velocity
            vp = sqrt(2.0*KB*Config.temp/M_CS)
            vx = -3.5*vp+6.0*vp*rand()
            vy = -3.5*vp+6.0*vp*rand()
            vz = -3.5*vp+6.0*vp*rand()
            ek = 0.5*M_CS*(vx*vx+vy*vy+vz*vz)/KB
            #etot = ek+U_t[y,x,1]-t_min
            eu = U11*(x_floor+1-x)*(y_floor+1-y)+U21*(x-x_floor)*(y_floor+1-y)+U12*(x_floor+1-x)*(y-y_floor)+U22*(x-x_floor)*(y-y_floor)
            etot = ek+eu-t_min #in K
            p = exp(-1.0*etot/Config.temp)
            if rand() <= p
                flag = true
                init_xv[i,1] = (x-1)*x_div
                init_xv[i,2] = (y-1)*y_div
                init_xv[i,3] = vx+Config.lattice_speed
                init_xv[i,4] = vy
            end
        end
    end
end

function pinitialize_wrap!(init_xv::SharedArray{Float64,2},x_grid,x_div,y_grid,y_div,U_t::SharedArray{Float64})
    @sync begin
        for p in procs(init_xv)
            @async begin
                remotecall_wait(p, pinitialize!, init_xv, x_grid,x_div,y_grid,y_div, U_t)
            end
        end
    end
end
#=
@everywhere @inbounds @fastmath function paccfunc2(t::Float64, x::Vector{Float64})
    arr = zeros(Float64,4)
    if check_wall(x[1:2]) == true
        fill!(arr,0.0)
    else
        #prepare index
        tindex = rem(t/acc_dt::Float64+1.0,Nf)
        t_b = floor(Integer,tindex)

        for i = 1:5
            tarr_orig[i] = t_b + i - 3
            tarr_orig_f[i] = Float64(t_b + i - 3)
            tarr[i] = t_b + i - 3
        end

        for i = 1:length(tarr)
            if tarr[i] < 1
                tarr[i] = N+tarr[i]
            elseif tarr[i]>N
                tarr[i] = tarr[i]-N
            end
        end
        interpolate!(acc_U::SharedArray{Float64,3},grad_arr,tarr,acc_dx::Float64,acc_dy::Float64,x)

        #interpolate in time
        for i = 1:5
            L[i] = 1
            for j = 1:5
                if i!=j
                    L[i] = L[i]*(tindex-tarr_orig_f[j])/(i-j)
                end
            end
        end
        arr[1] = x[3]
        arr[2] = x[4]
        for i = 1:5
            arr[3] = arr[3]+L[i]*grad_arr[i,1]
            arr[4] = arr[4]+L[i]*grad_arr[i,2]
        end
    end
    return arr
end

@everywhere @inbounds @fastmath function paccfunc2!(t::Float64, x::Vector{Float64},arr::Vector{Float64})
    if check_wall(x[1:2]) == true
        fill!(arr,0.0)
    else
        #prepare index
        tindex = rem(t/acc_dt::Float64+1.0,Nf)
        t_b = floor(Integer,tindex)

        for i = 1:5
            tarr_orig[i] = t_b + i - 3
            tarr_orig_f[i] = Float64(t_b + i - 3)
            tarr[i] = t_b + i - 3
        end

        for i = 1:length(tarr)
            if tarr[i] < 1
                tarr[i] = N+tarr[i]
            elseif tarr[i]>N
                tarr[i] = tarr[i]-N
            end
        end
        interpolate!(acc_U::SharedArray{Float64,3},grad_arr,tarr,acc_dx::Float64,acc_dy::Float64,x)

        #interpolate in time
        for i = 1:5
            L[i] = 1
            for j = 1:5
                if i!=j
                    L[i] = L[i]*(tindex-tarr_orig_f[j])/(i-j)
                end
            end
        end
        arr[1] = x[3]
        arr[2] = x[4]
        arr[3] = 0.0
        arr[4] = 0.0
        for i = 1:5
            arr[3] = arr[3]+L[i]*grad_arr[i,1]
            arr[4] = arr[4]+L[i]*grad_arr[i,2]
        end
    end
end

=#
end
