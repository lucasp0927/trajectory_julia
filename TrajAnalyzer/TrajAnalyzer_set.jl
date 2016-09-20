function init_parallel!(result::Dict,probe_sfn::ScalarFieldNode,ForceFields_sfn::ScalarFieldNode,config::Dict)
    traj_s = copy_to_sharedarray!(result["traj"])
    result_wo_traj = result
    delete!(result_wo_traj,"traj")
    @sync begin
        for p = 1:nprocs()
            Fields.clean_scaling!(probe_sfn)
            Fields.clean_scaling!(ForceFields_sfn)
            @async remotecall_wait(init!,p,result_wo_traj,traj_s,probe_sfn,ForceFields_sfn,config)
        end
    end
    Fields.eval_scaling!(probe_sfn)
    Fields.eval_scaling!(ForceFields_sfn)
end

function init!(result::Dict,traj_s::SharedArray{Float64},probe_sfn::ScalarFieldNode,ForceFields_sfn::ScalarFieldNode,config::Dict)
    global Trajs, Probe, ForceFields, TA_Config
    global avg_atom_num,lattice_width,lattice_unit,k_ratio,gamma_1d,gamma_prime
    Trajs = Trajectories(result,traj_s)
    Probe = Fields.copyfield(probe_sfn)
    ForceFields = Fields.copyfield(ForceFields_sfn)
    TA_Config = config
    avg_atom_num = Int64(TA_Config["spectrum"]["avg-atom-number"])
    lattice_width = Float64(TA_Config["spectrum"]["lattice-width"])
    lattice_unit =  Float64(TA_Config["spectrum"]["lattice-unit"])
    k_ratio = Float64(TA_Config["spectrum"]["k-ratio"])
    gamma_1d = Float64(TA_Config["spectrum"]["gamma-1d"])
    gamma_prime = Float64(TA_Config["spectrum"]["gamma-prime"])
end
