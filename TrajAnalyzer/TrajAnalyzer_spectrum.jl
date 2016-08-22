include("TrajAnalyzer_transfermatrix.jl")
using PyCall
function spectrum(filename)
    output,output_matrix = calculate_transmission()
    spectrum_data = Dict(
                    "output"=>output,
                    "output_matrix"=>output_matrix
                         )
    average_spectrum = squeeze(mean(abs2(output),3),3)
    h5_filename = filename*"_avg_spectrum_tm.h5"
    if isfile(h5_filename) == true
        rm(h5_filename)
    end
    h5write(h5_filename, "/spectrum", average_spectrum)
    freq_config = TA_Config["spectrum"]["frequency"]
    time_config = TA_Config["spectrum"]["time"]
    fstart = Float64(freq_config["start"])
    fend   = Float64(freq_config["end"])
    tstart = Float64(time_config["start"])
    tend = Float64(time_config["end"])
    #plot using matplotlib
    pyimport("matplotlib")[:use]("Agg")
    @pyimport matplotlib.pyplot as plt
    @pyimport numpy as np
    freq_range = Float64(freq_config["start"]):Float64(freq_config["step"]):Float64(freq_config["end"])
    time_range = Float64(time_config["start"]):Float64(time_config["step"]):Float64(time_config["end"])
    plt.pcolormesh(np.asarray(time_range),np.asarray(freq_range),average_spectrum,cmap="RdBu")
    plt.axis([tstart, tend, fstart, fend])
    plt.colorbar()
    plt.xlabel("time(us)")
    plt.ylabel("detuning (MHz)")
    plt.savefig(filename*"_spectrum.png")
    plt.clf()
end

function calculate_transmission()
    # traj selection
#    traj_selected = select_traj(paras.traj,paras.t_span,paras.select)
    # calculate transmission
#    init_parallel(paras.phy_para,paras.sim_para,paras.time_arr,paras.freq_arr,paras.t_span,traj_selected,paras.probe)
    #preallocate matrix
    freq_config = TA_Config["spectrum"]["frequency"]
    time_config = TA_Config["spectrum"]["time"]
    freq_range = Float64(freq_config["start"]):Float64(freq_config["step"]):Float64(freq_config["end"])
    time_range = Float64(time_config["start"]):Float64(time_config["step"]):Float64(time_config["end"])
    @assert Trajs.tspan[1]<time_range[1] && Trajs.tspan[end]>time_range[end] "spectrum time range out of range!"
    iter = TA_Config["spectrum"]["iteration"]
    output = SharedArray(Complex{Float64},(length(freq_range),length(time_range),iter))
    output_matrix = SharedArray(Complex{Float64},(2,2,length(freq_range),length(time_range),iter))
    for i in 1:iter
        Lumberjack.info("iteration: $i")
        lattice_sites::Float64 = lattice_width/lattice_unit
        #TODO: other distribution of atom_num
        atom_num = avg_atom_num::Int64
        Lumberjack.info("atom number: $atom_num")
        atom_arr::Array{Int64,2} = generate_atom_array(atom_num,Trajs.atom_num,lattice_sites)
@time        @sync @parallel for fidx in collect(eachindex(freq_range))
            for tidx in eachindex(time_range)
                output_matrix[:,:,fidx,tidx,i],output[fidx,tidx,i] = transmission(time_range[tidx],freq_range[fidx],atom_arr)
            end
        end
    end
    return sdata(output), sdata(output_matrix)
end

function generate_atom_array(atom_num,total_atom_num,lattice_sites)
    #lattice_scale = lattice_width/lattice_unit
    #(traj id,lattice id)*atom_num
    atom_array = zeros(Int64,(2,atom_num))
    atom_array[1,:] = shuffle(collect(1:total_atom_num))[1:atom_num] #atom id
    atom_array[2,:] = sort(round(Int64,randn(atom_num)*lattice_sites)) #atom's position in unit of lattice width
    return atom_array
end

@inbounds function transmission(t::Float64,detune::Float64,atom_arr::Array{Int64,2})
    # calculate transmission at time t and detuning detune.
    atom_num = size(atom_arr,2)
    x_point_k::Float64 = pi/lattice_unit::Float64
    k::Float64 = x_point_k*k_ratio::Float64
    #generate waveguide transfer matrix
    ldiff::Vector{Float64} = diff(squeeze(atom_arr[2,:],1))*lattice_unit::Float64
    M_wg::Array{Complex{Float64},3} = reduce((x,y)->cat(3,x,y),map(x->wg_transfer_matrix(k,x),ldiff))
    #generate atomic transfer matrix
    M_atom::Array{Complex{Float64},3} = zeros(Complex{Float64},(2,2,atom_num))
    for i = 1:atom_num
        pos::Vector{Float64} = Trajs[t,atom_arr[1,i]]
        if any(isnan(pos[1:2]))
            M_atom[:,:,i] = wg_transfer_matrix(0.0,0.0)
        else
            f_0 = Fields.value(pos[1:2],t,ForceFields::ScalarFieldNode)/(-1e-3)*20
            p_0 = Fields.value(pos[1:2],t,Probe::ScalarFieldNode)
#            @assert p_0 >= 0.0 "negative probe power!"
            M_atom[:,:,i] = atom_transfer_matrix(detune,f_0,p_0*gamma_1d::Float64,gamma_prime::Float64)
        end
    end
    M_tot::Array{Complex{Float64},2} = M_atom[:,:,1];
    @fastmath for i = 1:atom_num - 1
        M_tot *= M_atom[:,:,i+1]*M_wg[:,:,i]
    end
    return M_tot,one(Complex{Float64})/M_tot[2,2]
end

#=
function parallel_transmission!(tran,output_matrix)
    #norm_time_arr = calc_norm_tarr(time_arr,t_span) #normalize time array
    pm = Progress(sim_para["iteration"], 1)
    np = nprocs()
    i = 1
    nextidx() = (idx=i; i+=1; idx)
    @sync begin
        for p=1:np
            if p != myid() || np == 1
                @async begin
                    while true
                        idx = nextidx()
                        if idx > sim_para["iteration"]
                            break
                        end
                        (tran[:,:,idx],output_matrix[:,:,:,:,idx]) = remotecall_fetch(transmission,p)
                        next!(pm)
                    end
                end
            end
        end
    end
end

@inbounds @fastmath function transmission()
    global norm_time_arr, freq_arr, sim_para, phy_para, traj, probe_itp
    global U_itp, X_itp, Y_itp
    # calculate transmission at index idx and frequency f.
    total_atom_num::Int64 = sim_para["use-all-atom"]?size(traj,3):sim_para["total-atom-number"]
    atom_num::Int64 = sim_para["use-all-atom"]?total_atom_num:sim_para["avg-atom-number"]
    atom_arr::Array{Int64,2} = generate_atom_array(atom_num,total_atom_num,phy_para["lattice-width"]/phy_para["lattice-unit"])
    M_atom::Array{Complex{Float64},3} = zeros(Complex{Float64},(2,2,atom_num))
    M_wg::Array{Complex{Float64},3} = zeros(Complex{Float64},(2,2,atom_num-1))
    x_point_k::Float64 = pi/phy_para["lattice-unit"]
    k::Float64 = x_point_k*phy_para["k-ratio"]
    for i = 1:atom_num-1
        ldiff = diff(squeeze(atom_arr[2,:],1))
        M_wg[:,:,i] = wg_transfer_matrix(k,Float64(ldiff[i]*phy_para["lattice-unit"]))
    end
    U_itp_arr = Array(Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},Interpolations.BSpline{Interpolations.Linear},Interpolations.OnGrid,0},atom_num)
    X_itp_arr = Array(Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},Interpolations.BSpline{Interpolations.Linear},Interpolations.OnGrid,0},atom_num)
    Y_itp_arr = Array(Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},Interpolations.BSpline{Interpolations.Linear},Interpolations.OnGrid,0},atom_num)
    for i = 1:atom_num
        #TODO: use permute
        X_itp_arr[i] = X_itp[atom_arr[1,i]]
        Y_itp_arr[i] = Y_itp[atom_arr[1,i]]
        U_itp_arr[i] = U_itp[atom_arr[1,i]]
    end
    result = Array(Complex{Float64},(length(freq_arr),length(norm_time_arr)))
    M_total = Array(Complex{Float64},(2,2,length(freq_arr),length(norm_time_arr)))
    transmission_loop!(result,M_total,atom_num::Int64,U_itp_arr,X_itp_arr,Y_itp_arr,M_atom,M_wg)
    return (result,M_total)
end

@inbounds function transmission_loop!(result::Array{Complex{Float64},2},M_total::Array{Complex{Float64},4},atom_num::Int64,U_itp_arr,X_itp_arr,Y_itp_arr,M_atom::Array{Complex{Float64},3},M_wg::Array{Complex{Float64},3})
    global norm_time_arr, freq_arr
    global c_factor, linewidth
    global probe_itp
    for t = enumerate(norm_time_arr), f = enumerate(freq_arr)
        for i = 1:atom_num::Int64
            #f_0 AC stark shift frequency
            #p power of probe beam
            f_0::Float64 = (U_itp_arr[i]::Interpolations.BSplineInterpolation{Float64,1,Array{Float64,1},Interpolations.BSpline{Interpolations.Linear},Interpolations.OnGrid,0})[t[2]::Float64]/(-1e-3)*20;
            if isnan(f_0)
                M_atom[:,:,i]::Array{Complex{Float64},2} = eye(Complex{Float64},2)
            else
                gmxy::Tuple{Float64,Float64} = gm_xy(X_itp_arr[i][t[2]::Float64],Y_itp_arr[i][t[2]::Float64])
                p::Float64 = (probe_itp::Interpolations.BSplineInterpolation{Float64,2,Array{Float64,2},Interpolations.BSpline{Interpolations.Linear},Interpolations.OnGrid,0})[gmxy[1]::Float64,gmxy[2]::Float64]
                #TODO: minus p value!?
                M_atom[:,:,i]::Array{Complex{Float64},2} = atom_transfer_matrix(f[2],f_0,p*c_factor::Float64,linewidth::Float64)::Array{Complex{Float64},2}
            end
        end
        M_tot::Array{Complex{Float64},2} = M_atom[:,:,1];
        @fastmath for i = 1:atom_num - 1
            M_tot *= M_atom[:,:,i+1]*M_wg[:,:,i]
        end
        (M_total::Array{Complex{Float64},4})[:,:,f[1]::Int64,t[1]::Int64] = M_tot
        (result::Array{Complex{Float64},2})[f[1]::Int64,t[1]::Int64] = one(Complex{Float64})/M_tot[2,2]
    end
end
=#
