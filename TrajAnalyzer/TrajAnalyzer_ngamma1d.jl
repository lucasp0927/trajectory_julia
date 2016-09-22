using Base.Test
@generated function ngamma1d(idx::Vector,filename)
    quote
        if length(idx) == 1
            a = $(zeros(Float64,length(range_i)))
        elseif length(idx) == 2
            a = $(zeros(Float64,length(range_i),length(range_j)))
        end
        #        @test_approx_eq calc_avg_ngamma1d_parallel() calc_avg_ngamma1d()
        a[idx...] = calc_avg_ngamma1d_parallel()
        matwrite(filename, Dict("ng1d" => a))
    end
end

function calc_avg_ngamma1d()
    ng1d = sum(calc_ngamma1d,Trajs.tspan)
    # ng1d = @parallel (+) for t in Trajs.tspan
    #     calc_ngamma1d(t)
    # end
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
    sum(i->calc_gamma1d(traj_snapshot[:,i],t),1:size(traj_snapshot,2))
end
