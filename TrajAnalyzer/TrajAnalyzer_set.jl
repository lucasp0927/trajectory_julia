function init_parallel!(result::Dict,probe_sfn::ScalarFieldNode,ForceFields_sfn::ScalarFieldNode,config::Dict)
    #filter trajectories
    filter_traj(result,config["filter"])
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
    global range_i, range_j
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
    if TA_Config["type"] == "single-scan-scaling"
        range_i = collect(Int(TA_Config["range_i_start"]):Int(TA_Config["range_i_end"]))
    elseif TA_Config["type"] == "double-scan-scaling"
        range_i = collect(Int(TA_Config["range_i_start"]):Int(TA_Config["range_i_end"]))
        range_j = collect(Int(TA_Config["range_j_start"]):Int(TA_Config["range_j_end"]))
    end
end
"""
    filter:
        type: crashed
        tstart: 1350.0
        tend: 1420.0
        gap: [9890.0, 10110.0, 10110.0, 9890.0, 25120.0, 25120.0, 24880.0, 24880.0]
"""
function filter_traj(result::Dict,config::Dict)
    if config["type"] == "none"
        info("Not filtering trajectories.")
    elseif config["type"] == "crashed"
        info("Select crashed trajectories.")
        start_idx = searchsortedlast(result["tspan"],config["tstart"])
        end_idx = searchsortedlast(result["tspan"],config["tend"])
        traj = result["traj"]
        pp = Polygon([promote(config["gap"]...)...])
        selected = filter(i->any(isnan(traj[:,start_idx:end_idx,i])),1:size(traj,3))
        selected = filter(i->anyPointInPolygon(pp,traj[1:2,start_idx:end_idx,i])==false,selected)
        info("selected $(length(selected)) trajectories from $(size(traj,3)) trajectories.")
        result["traj"] = cat(3,map(i->traj[:,:,i],selected)...)
    elseif config["type"] == "gap"
        info("Select trajectories pass through the gap.")
        start_idx = searchsortedlast(result["tspan"],config["tstart"])
        end_idx = searchsortedlast(result["tspan"],config["tend"])
        traj = result["traj"]
        pp = Polygon([promote(config["gap"]...)...])
        #first remove crashed atom
        selected = filter(i->any(isnan(traj[:,start_idx:end_idx,i]))==false,1:size(traj,3))
        #select trajectories pass through the gap
        selected = filter(i->anyPointInPolygon(pp,traj[1:2,start_idx:end_idx,i]),selected)
        info("selected $(length(selected)) trajectories from $(size(traj,3)) trajectories.")
        result["traj"] = cat(3,map(i->traj[:,:,i],selected)...)
    end
end
