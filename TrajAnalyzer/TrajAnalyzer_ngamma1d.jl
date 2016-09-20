@generated function ngamma1d(idx::Vector,filename)
    quote
        if length(idx) == 1
            a = $(zeros(Float64,length(range_i)))
        elseif length(idx) == 2
            a = $(zeros(Float64,length(range_i),length(range_j)))
        end
        a[idx...] = calc_ngamma1d()
        matwrite(filename, Dict("ng1d" => a))
    end
end

function calc_ngamma1d()
    ng1d = 0.0
    for t in Trajs.tspan
        traj_snapshot = Trajs[171.0,:]
        traj_snapshot = traj_snapshot[1:2,:]
        for i in 1:Trajs.atom_num
            ng1d += calc_gamma1d(traj_snapshot[:,i],t)
        end
    end
    println(ng1d)
    ng1d = ng1d/length(Trajs.tspan)
end
