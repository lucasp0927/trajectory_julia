include("TrajAnalyzer_transfermatrix.jl")
using PyCall
using StatsBase
using Distributions

function spectrum(filename)
    output,output_matrix = calculate_transmission()
    spectrum_data = Dict(
                    "output"=>output,
                    "output_matrix"=>output_matrix
                         )
    average_spectrum = squeeze(mean(abs2(output),3),3)

    h5_filename = filename*"_avg_spectrum.h5"
    if isfile(h5_filename) == true
        rm(h5_filename)
    end
    # h5open(h5_filename,"w") do file
    #     write(file,"avg_spectrum",average_spectrum)
    #     write(file,"spectrum",output)
    #     write(file,"transfer_matrix",output_matrix)
    # end
    h5write(h5_filename, "/spectrum", output)
    h5write(h5_filename, "/avg_spectrum", average_spectrum)
    h5write(h5_filename, "/transfer_matrix",output_matrix)
    #plot using matplotlib
    freq_config = TA_Config["spectrum"]["frequency"]
    time_config = TA_Config["spectrum"]["time"]
    fstart = Float64(freq_config["start"])
    fend   = Float64(freq_config["end"])
    tstart = Float64(time_config["start"])
    tend = Float64(time_config["end"])
    freq_range = fstart:Float64(freq_config["step"]):fend
    time_range = tstart:Float64(time_config["step"]):tend
    pyimport("matplotlib")[:use]("Agg")
    @pyimport matplotlib.pyplot as plt
    @pyimport numpy as np
    plt.pcolormesh(np.asarray(time_range),np.asarray(freq_range),average_spectrum,cmap="RdBu")
    plt.axis([tstart, tend, fstart, fend])
    plt.colorbar()
    plt.xlabel("time(us)")
    plt.ylabel("detuning (MHz)")
    plt.savefig(filename*"_spectrum.png")
    plt.clf()
end

@inbounds function calculate_transmission()
    # traj selection
#    traj_selected = select_traj(paras.traj,paras.t_span,paras.select)
    # calculate transmission
#    init_parallel(paras.phy_para,paras.sim_para,paras.time_arr,paras.freq_arr,paras.t_span,traj_selected,paras.probe)
    #preallocate matrix
    freq_config = TA_Config["spectrum"]["frequency"]
    time_config = TA_Config["spectrum"]["time"]
    freq_range = Float64(freq_config["start"]):Float64(freq_config["step"]):Float64(freq_config["end"])
    time_range = Float64(time_config["start"]):Float64(time_config["step"]):Float64(time_config["end"])
    @assert Trajs.tspan[1]<=time_range[1] && Trajs.tspan[end]>=time_range[end] "spectrum time range out of range!"
    iter = TA_Config["spectrum"]["iteration"]
    output = SharedArray(Complex{Float64},(length(freq_range),length(time_range),iter))
    output_matrix = SharedArray(Complex{Float64},(2,2,length(freq_range),length(time_range),iter))
    for i in 1:iter
        info("iteration: $i")
        lattice_scale::Float64 = lattice_width/lattice_unit
        #TODO: other distribution of atom_num
        atom_num = round(Int,avg_atom_num*(Trajs.atom_num/TA_Config["spectrum"]["total-atom-number"]))
        @assert atom_num <= Trajs.atom_num
        info("atom number: $atom_num")
        atom_arr::Array{Int64,2} = generate_atom_array(atom_num,Trajs.atom_num,lattice_scale)
        @time @sync @parallel for fidx in collect(eachindex(freq_range))
            for tidx in eachindex(time_range)
                output_matrix[:,:,fidx,tidx,i],output[fidx,tidx,i] = transmission(time_range[tidx],freq_range[fidx],atom_arr)
            end
        end
    end
    return sdata(output), sdata(output_matrix)
end

function generate_atom_array(atom_num,total_atom_num,lattice_scale)
    #lattice_scale = lattice_width/lattice_unit
    #(traj id,lattice id)*atom_num
    atom_array = zeros(Int64,(2,atom_num))
    atom_array[1,:] = sample(collect(1:total_atom_num),atom_num)
    #    atom_array[1,:] = shuffle(collect(1:total_atom_num))[1:atom_num] #atom id
    d = Truncated(Normal(lattice_scale/2, atom_beam_waist/lattice_unit), 0, lattice_scale)
    atom_pos = round(Int64,rand(d,atom_num))
    @assert all(lattice_scale .>= atom_pos .>= 0)
    atom_array[2,:] = sort(atom_pos) #atom's position in unit of lattice width
#    atom_array[2,:] = sort(round(Int64,randn(atom_num)*lattice_scale)) #atom's position in unit of lattice width
    return atom_array
end

function calc_gamma1d(pos,t)
    @assert any(isnan(pos)) == false
    g1d = Float64(Fields.value(pos[1:2],t,Probe::ScalarFieldNode)*gamma_1d)
    @assert isnan(g1d) == false
    return g1d
end

@inbounds function transmission(t::Float64,detune::Float64,atom_arr::Array{Int64,2})
    # calculate transmission at time t and detuning detune.
    # put wg from -2*Lattice_unit to first atom, and after last atom to lattice_width+2*lattice_unit
    atom_num = size(atom_arr,2)
    x_point_k::Float64 = pi/lattice_unit::Float64
    k::Float64 = x_point_k*k_ratio::Float64
    #generate waveguide transfer matrix
    d = Truncated(Normal(0, sqrt(pos_variance)), -1, 1)
    atom_pos_div = rand(d,atom_num)
    atom_pos = (atom_arr[2,:] + atom_pos_div)*lattice_unit
    @assert all(lattice_width+2*lattice_unit .> atom_pos .> -2*lattice_unit)
    ldiff::Vector{Float64} = diff(atom_pos)

#    println(size(ldiff))
#    M_wg::Array{Complex{Float64},3} = reduce((x,y)->cat(3,x,y),map(x->wg_transfer_matrix(k,x),ldiff))
#    println(size(M_wg))
    M_wg = zeros(Complex{Float64},2,2,length(ldiff))
    for i = 1:length(ldiff)
        M_wg[:,:,i] = wg_transfer_matrix(k,ldiff[i])
    end
    #generate atomic transfer matrix
    M_atom::Array{Complex{Float64},3} = zeros(Complex{Float64},(2,2,atom_num))
    for i = 1:atom_num
        pos::Vector{Float64} = Trajs[t,atom_arr[1,i]]
        if any(isnan(pos[1:2]))
            M_atom[:,:,i] = wg_transfer_matrix(0.0,0.0)
        else
            f_0 = Fields.value(pos[1:2],t,ForceFields::ScalarFieldNode)*(-2.08e4) #*20.8/(-1e-3)
            ### vector shift?
            g1d = calc_gamma1d(pos[1:2],t)*cospi(atom_pos_div[i])^2
#            @assert p_0 >= 0.0 "negative probe power!"
            M_atom[:,:,i] = atom_transfer_matrix(detune,f_0,g1d,gamma_prime::Float64)
        end
    end

    first_wg = wg_transfer_matrix(k,atom_pos[1]-(-2lattice_unit))
    last_wg = wg_transfer_matrix(k,lattice_width+2lattice_unit-atom_pos[end])
    M_tot::Array{Complex{Float64},2} = first_wg*M_atom[:,:,1];
    @fastmath for i = 1:atom_num - 1
        M_tot *= M_atom[:,:,i+1]*M_wg[:,:,i]
    end
    M_tot *= last_wg
    return M_tot,one(Complex{Float64})/M_tot[2,2]
end
