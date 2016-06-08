function set_trajs_parallel!(result::Dict)
    traj_s = copy_to_sharedarray!(result["traj"])
    @sync begin
        for p = 1:nprocs()
            @async remotecall_wait(p,set_trajs!,result,traj_s)
        end
    end
end

function set_trajs!(result::Dict,traj_s::SharedArray{Float64})
    global Trajs
    Trajs = Trajectories(result,traj_s)
end
