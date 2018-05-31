using Base.Test
@generated function ngamma1d(idx::Vector,filename)
    quote
        if length(idx) == 1
            a = $(zeros(Float64,length(range_i)))
            idx_a = [searchsorted(range_i,idx[1])]
            a[idx_a...] = calc_avg_ngamma1d_parallel()
            if idx[1] == range_i[end]
                @info "saving ngamma1d..."
                matwrite(filename, Dict("ng1d" => a))
            end
        elseif length(idx) == 2
            a = $(zeros(Float64,length(range_i),length(range_j)))
            idx_a = [searchsorted(range_i,idx[1]),searchsorted(range_j,idx[2])]
            a[idx_a...] = calc_avg_ngamma1d_parallel()
            if idx[1] == range_i[end] && idx[2] == range_j[end]
                @info "saving ngamma1d..."
                matwrite(filename, Dict("ng1d" => a))
            end
        end
    end
end

function calc_avg_ngamma1d()
    ng1d = sum(calc_ngamma1d,Trajs.tspan)
    ng1d = ng1d/length(Trajs.tspan)
end

function calc_avg_ngamma1d_parallel()
    ng1d = @parallel (+) for t in Trajs.tspan
        calc_ngamma1d(t)
    end
    ng1d = ng1d/length(Trajs.tspan)
end

function calc_ngamma1d(t::Float64)
    traj_snapshot = Trajs[t,:]
    #remove NaN coordinate for vanished atom
    traj_snapshot = removenan(traj_snapshot)
    @assert any(isnan(traj_snapshot)) == false
    traj_snapshot = traj_snapshot[1:2,:]
    sum(calc_gamma1d(traj_snapshot[:,i],t) for i in 1:size(traj_snapshot,2))
end
