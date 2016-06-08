function set_trajs_parallel!(result::Dict,probe_sfn::Fields.ScalarFieldNode)
    traj_s = copy_to_sharedarray!(result["traj"])
    @sync begin
        for p = 1:nprocs()
            @async remotecall_wait(p,set_trajs!,result,traj_s,probe_sfn)
        end
    end
end

function set_trajs!(result::Dict,traj_s::SharedArray{Float64},probe_sfn::Fields.ScalarFieldNode)
    global Trajs, Probe
    Trajs = Trajectories(result,traj_s)
    Probe = Fields.copyfield(probe_sfn)
end
