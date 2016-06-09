function init_parallel!(result::Dict,probe_sfn::ScalarFieldNode,ForceFields_sfn::ScalarFieldNode,config::Dict)
    traj_s = copy_to_sharedarray!(result["traj"])
    @sync begin
        for p = 1:nprocs()
            @async remotecall_wait(p,init!,result,traj_s,probe_sfn,ForceFields_sfn,config)
        end
    end
end

function init!(result::Dict,traj_s::SharedArray{Float64},probe_sfn::ScalarFieldNode,ForceFields_sfn::ScalarFieldNode,config::Dict)
    global Trajs, Probe, ForceFields, TA_Config
    Trajs = Trajectories(result,traj_s)
    Probe = Fields.copyfield(probe_sfn)
    ForceFields = Fields.copyfield(ForceFields_sfn)
    TA_Config = config
end
