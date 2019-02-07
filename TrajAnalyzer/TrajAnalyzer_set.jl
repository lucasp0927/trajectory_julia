function init_probe_parallel!(probe_sfn::ScalarFieldNode)
    @sync begin
        for p = 1:nprocs()
            Fields.clean_scaling!(probe_sfn)
            @async remotecall_wait(init_probe!,p,probe_sfn)
        end
    end
    Fields.eval_scaling!(probe_sfn)
end

function init_probe!(probe_sfn::ScalarFieldNode)
    global Probe
    Probe = Fields.copyfield(probe_sfn)
end


function init_parallel!(result::Dict,probe_sfn::ScalarFieldNode,ForceFields_sfn::ScalarFieldNode,config::Dict)
    #filter trajectories
    filter_traj(result,config["filter"])
    traj_s = copy_to_sharedarray!(result["traj"])
    result_wo_traj = result
    result = nothing
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
    #ploting backend
#    @info "set gr() as Plots.jl backend."
#    gr()
    global spectrum_mode, vector_shift
    global Trajs, Probe, ForceFields, TA_Config
    global avg_atom_num,lattice_width,lattice_unit,k_ratio,gamma_1d,gamma_prime,pos_variance,atom_beam_waist, probe_contrast
    global range_i, range_j
    Trajs = Trajectories(result,traj_s)
    Probe = Fields.copyfield(probe_sfn)
    ForceFields = Fields.copyfield(ForceFields_sfn)
    TA_Config = config
    spectrum_mode = Int64(TA_Config["spectrum"]["mode"])
    vector_shift = Int64(TA_Config["spectrum"]["vector-shift"])
    avg_atom_num = Int64(TA_Config["spectrum"]["avg-atom-number"])
    lattice_width = Float64(TA_Config["spectrum"]["lattice-width"])
    lattice_unit =  Float64(TA_Config["spectrum"]["lattice-unit"])
    k_ratio = Float64(TA_Config["spectrum"]["k-ratio"])
    gamma_1d = Float64(TA_Config["spectrum"]["gamma-1d"])
    gamma_prime = Float64(TA_Config["spectrum"]["gamma-prime"])
    pos_variance = Float64(TA_Config["spectrum"]["pos-variance"]) #position variance in lattice unit
    atom_beam_waist = Float64(TA_Config["spectrum"]["atom-beam-waist"])
    probe_contrast = Float64(TA_Config["spectrum"]["probe-contrast"])
    if TA_Config["type"] == "single-scan-scaling"
        range_i = collect(Int(TA_Config["range_i_start"]):Int(TA_Config["range_i_end"]))
    elseif TA_Config["type"] == "double-scan-scaling"
        range_i = collect(Int(TA_Config["range_i_start"]):Int(TA_Config["range_i_end"]))
        range_j = collect(Int(TA_Config["range_j_start"]):Int(TA_Config["range_j_end"]))
    end
end

function cleanup()
    global Trajs
    finalize_shared_array!(Trajs.traj)
    @everywhere GC.gc()
end
