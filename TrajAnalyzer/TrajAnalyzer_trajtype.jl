function traj_iscrashed()
    config = TA_Config["filter"]
    start_idx = searchsortedlast(Trajs.tspan,config["tstart"])
    end_idx = searchsortedlast(Trajs.tspan,config["tend"])
    gap_p = Polygon([promote(config["gap"]...)...])
    roi_p = Polygon([promote(config["roi"]...)...])
    # crashed (may pass through the gap)
    selected = collect(1:size(Trajs.traj,3))
    selected = filter(i->anyPointInPolygon(roi_p,Trajs.traj[1:2,start_idx:end_idx,i]),selected)
    selected = filter(i->any(isnan(Trajs.traj[:,start_idx:end_idx,i])),selected)
    selected = filter(i->all(isnan(Trajs.traj[:,start_idx:end_idx,i]))==false,selected)
    selected = filter(i->anyPointInPolygon(gap_p,Trajs.traj[1:2,start_idx:end_idx,i])==false,selected)
    return length(selected)
end

function traj_ingap(exclude_crashed::Bool)
    config = TA_Config["filter"]
    start_idx = searchsortedlast(Trajs.tspan,config["tstart"])
    end_idx = searchsortedlast(Trajs.tspan,config["tend"])
    gap_p = Polygon([promote(config["gap"]...)...])
    roi_p = Polygon([promote(config["roi"]...)...])
    # crashed (may pass through the gap)
    selected = collect(1:size(Trajs.traj,3))
    selected = filter(i->anyPointInPolygon(roi_p,Trajs.traj[1:2,start_idx:end_idx,i]),selected)
    if exclude_crashed
        selected = filter(i->any(isnan(Trajs.traj[:,start_idx:end_idx,i]))==false,selected)
    end
    selected = filter(i->anyPointInPolygon(gap_p,Trajs.traj[1:2,start_idx:end_idx,i]), selected)
    return length(selected)
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