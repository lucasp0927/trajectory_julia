function traj_iscrashed()
    #need optimization
    config = TA_Config["filter"]
    start_idx = searchsortedlast(Trajs.tspan,config["tstart"])
    end_idx = searchsortedlast(Trajs.tspan,config["tend"])
    gap_p = Polygon([promote(config["gap"]...)...])
    roi_p = Polygon([promote(config["roi"]...)...])
    selected = collect(1:size(Trajs.traj,3))
    selected = filter(i->anyPointInPolygon(roi_p,Trajs.traj[1:2,start_idx:end_idx,i]),selected)
    selected = filter(i->any(isnan.(Trajs.traj[:,start_idx:end_idx,i])),selected)
    selected = filter(i->all(isnan.(Trajs.traj[:,start_idx:end_idx,i]))==false,selected)
    selected = filter(i->anyPointInPolygon(gap_p,Trajs.traj[1:2,start_idx:end_idx,i])==false,selected)
    return selected
end

function traj_ingap(exclude_crashed::Bool)
    #need optimization
    config = TA_Config["filter"]
    start_idx = searchsortedlast(Trajs.tspan,config["tstart"])
    end_idx = searchsortedlast(Trajs.tspan,config["tend"])
    gap_p = Polygon([promote(config["gap"]...)...])
    roi_p = Polygon([promote(config["roi"]...)...])
    selected = collect(1:size(Trajs.traj,3))
    selected = filter(i->anyPointInPolygon(roi_p,Trajs.traj[1:2,start_idx:end_idx,i]),selected)
    if exclude_crashed
        selected = filter(i->any(isnan(Trajs.traj[:,start_idx:end_idx,i]))==false,selected)
    end
    selected = filter(i->anyPointInPolygon(gap_p,Trajs.traj[1:2,start_idx:end_idx,i]), selected)
    return selected
end

"""
    filter:
        type: crashed
        tstart: 1350.0
        tend: 1420.0
        gap: [9890.0, 10110.0, 10110.0, 9890.0, 25120.0, 25120.0, 24880.0, 24880.0]
"""
function filter_traj(result::Dict,config::Dict)
    start_idx = searchsortedlast(vec(result["tspan"]),config["tstart"])
    end_idx = searchsortedlast(vec(result["tspan"]),config["tend"])
    traj = result["traj"]
    gap_p = Polygon([promote(config["gap"]...)...])
    side1_p = Polygon([promote(config["side1"]...)...])
    side2_p = Polygon([promote(config["side2"]...)...])
    if config["type"] == "none"
        @info "Not filtering trajectories."
    elseif config["type"] == "gap"
        @info "Select trajectories pass through the gap."
        #all trajectories ID
        selected = collect(1:size(traj,3))
        selected = filter(i->anyPointInPolygon(gap_p,traj[1:2,start_idx:end_idx,i]), selected)
        @info "selected $(length(selected)) trajectories from $(size(traj,3)) trajectories."
        result["traj"] = cat(map(i->traj[:,:,i],selected)...,dims=3)
    elseif config["type"] == "sides"
        @info "Select side trajectories."
        selected = collect(1:size(traj,3))
        selected1 = filter(i->anyPointInPolygon(side1_p,traj[1:2,start_idx:end_idx,i]), selected)
        selected2 = filter(i->anyPointInPolygon(side2_p,traj[1:2,start_idx:end_idx,i]), selected)
        selected = cat(selected1,selected2,dims=1)
        @info "selected $(length(selected)) trajectories from $(size(traj,3)) trajectories."
        result["traj"] = cat(map(i->traj[:,:,i],selected)...,dims=3)
    elseif config["type"] == "other"
        @info "Select top-bottom trajectories."
        #not in "gap" and not in "sides"
        selected = collect(1:size(traj,3))
        selected = filter(i->anyPointInPolygon(gap_p,traj[1:2,start_idx:end_idx,i])==false, selected)
        selected = filter(i->anyPointInPolygon(side1_p,traj[1:2,start_idx:end_idx,i])==false, selected)
        selected = filter(i->anyPointInPolygon(side2_p,traj[1:2,start_idx:end_idx,i])==false, selected)
        @info "selected $(length(selected)) trajectories from $(size(traj,3)) trajectories." 
        result["traj"] = cat(map(i->traj[:,:,i],selected)...,dims=3)
    end

#     if config["type"] == "none"
#         @info "Not filtering trajectories."
#     elseif config["type"] == "crashed"
#         @info "Select crashed trajectories."
#         start_idx = searchsortedlast(result["tspan"],config["tstart"])
#         end_idx = searchsortedlast(result["tspan"],config["tend"])
#         traj = result["traj"]
#         gap_p = Polygon([promote(config["gap"]...)...])
#         roi_p = Polygon([promote(config["roi"]...)...])
#         selected = collect(1:size(traj,3))
#         selected = filter(i->anyPointInPolygon(roi_p,traj[1:2,start_idx:end_idx,i]),selected)
#         selected = filter(i->any(isnan(traj[:,start_idx:end_idx,i])),selected)
#         selected = filter(i->all(isnan(traj[:,start_idx:end_idx,i]))==false,selected)
#         selected = filter(i->anyPointInPolygon(gap_p,traj[1:2,start_idx:end_idx,i])==false,selected)
#         @info "selected $(length(selected)) trajectories from $(size(traj,3)) trajectories."
#         result["traj"] = cat(3,map(i->traj[:,:,i],selected)...)
#     elseif config["type"] == "gap"
#         @info "Select trajectories pass through the gap."
#         start_idx = searchsortedlast(result["tspan"],config["tstart"])
#         end_idx = searchsortedlast(result["tspan"],config["tend"])
#         traj = result["traj"]
#         gap_p = Polygon([promote(config["gap"]...)...])
#         roi_p = Polygon([promote(config["roi"]...)...])
#         selected = collect(1:size(traj,3))
#         selected = filter(i->anyPointInPolygon(roi_p,traj[1:2,start_idx:end_idx,i]),selected)
#         #if exclude_crashed
#         #    selected = filter(i->any(isnan(traj[:,start_idx:end_idx,i]))==false,selected)
#         #end
#         selected = filter(i->anyPointInPolygon(gap_p,traj[1:2,start_idx:end_idx,i]), selected)
#         @info "selected $(length(selected)) trajectories from $(size(traj,3)) trajectories."
#         result["traj"] = cat(3,map(i->traj[:,:,i],selected)...)
#     elseif config["type"] == "non-gap"
#         @info "Select trajectories dont pass through the gap."
#         start_idx = searchsortedlast(result["tspan"],config["tstart"])
#         end_idx = searchsortedlast(result["tspan"],config["tend"])
#         traj = result["traj"]
#         gap_p = Polygon([promote(config["gap"]...)...])
#         roi_p = Polygon([promote(config["roi"]...)...])
#         selected = collect(1:size(traj,3))
# #        selected = filter(i->anyPointInPolygon(roi_p,traj[1:2,start_idx:end_idx,i]),selected)
#         #if exclude_crashed
#         #    selected = filter(i->any(isnan(traj[:,start_idx:end_idx,i]))==false,selected)
#         #end
#         selected = filter(i->anyPointInPolygon(gap_p,traj[1:2,start_idx:end_idx,i])==false, selected)
#         @info "selected $(length(selected)) trajectories from $(size(traj,3)) trajectories."
#         result["traj"] = cat(3,map(i->traj[:,:,i],selected)...)
#     end
end
