#using PyCall
function spectrum3d(filename,gm_name)
    output,output_matrix = calculate_transmission3d()
    average_spectrum = dropdims(mean(abs2.(output),dims=3),dims=3)
    spectrum_data = Dict(
                         "spectrum"=>output,
                         "avg_spectrum"=>average_spectrum,
                         "transfer_matrix"=>output_matrix
                         )
    dicttoh5(filename*"_spectrum_data_"*gm_name*".h5",spectrum_data)
    #plot using matplotlib
    # freq_config = TA_Config["spectrum"]["frequency"]
    # time_config = TA_Config["spectrum"]["time"]
    # fstart = Float64(freq_config["start"])
    # fend   = Float64(freq_config["end"])
    # tstart = Float64(time_config["start"])
    # tend = Float64(time_config["end"])
    # freq_range = fstart:Float64(freq_config["step"]):fend
    # time_range = tstart:Float64(time_config["step"]):tend
    # TODO: fix pyimport for julia 1.0
    # pyimport("matplotlib")[:use]("Agg")
    # @pyimport matplotlib.pyplot as plt
    # @pyimport numpy as np
    # plt.pcolormesh(np.asarray(time_range),np.asarray(freq_range),average_spectrum,cmap="RdBu")
    # plt.axis([tstart, tend, fstart, fend])
    # plt.colorbar()
    # plt.xlabel("time(us)")
    # plt.ylabel("detuning (MHz)")
    # plt.savefig(filename*"_spectrum_"*gm_name*".png")
    # plt.clf()
end

@inbounds function calculate_transmission3d()
    # traj selection
#    traj_selected = select_traj(paras.traj,paras.t_span,paras.select)
    # calculate transmission
#    init_parallel(paras.phy_para,paras.sim_para,paras.time_arr,paras.freq_arr,paras.t_span,traj_selected,paras.probe)
    #preallocate matrix
    @debug "in calculate_transmission3d()"
    freq_config = TA_Config["spectrum"]["frequency"]
    time_config = TA_Config["spectrum"]["time"]
    freq_range = Float64(freq_config["start"]):Float64(freq_config["step"]):Float64(freq_config["end"])
    time_range = Float64(time_config["start"]):Float64(time_config["step"]):Float64(time_config["end"])
    @assert Trajs.tspan[1]<=time_range[1] && Trajs.tspan[end]>=time_range[end] "spectrum time range out of range!"
    iter = TA_Config["spectrum"]["iteration"]
    output = SharedArray{Complex{Float64}}(length(freq_range),length(time_range),iter)
    output_matrix = SharedArray{Complex{Float64}}((2,2,length(freq_range),length(time_range),iter))
    # if spectrum_mode == 1
    #     d = Uniform(-1,1)
    # elseif spectrum_mode == 2
    #     d = Truncated(Normal(0, sqrt(pos_variance)), -1, 1)
    # end
    @time @sync @distributed for i in 1:iter
        @debug "iteration: $i"
        lattice_scale::Float64 = lattice_width/lattice_unit
        #TODO: other distribution of atom_num
        atom_num = round(Int,avg_atom_num*(Trajs.atom_num/TA_Config["spectrum"]["total-atom-number"]))
        @debug "atom_num:", atom_num
        @assert atom_num <= Trajs.atom_num
        @debug "generating_atom_array..."
        atom_arr::Array{Int64,2} = generate_atom_array3d(atom_num,Trajs.atom_num,lattice_scale)
        #@sync @parallel for fidx in collect(eachindex(freq_range))
        #preallocate array
        M_wg = Array{Complex{Float64}}(undef,2,2,atom_num-1)
        M_atom::Array{Complex{Float64},3} = Array{Complex{Float64}}(undef,2,2,atom_num)
        @debug "calculating output_matrix..."
        for fidx in collect(eachindex(freq_range))
            for tidx in eachindex(time_range)
                @debug "fidx, tix: ", fidx, tidx
                output_matrix[:,:,fidx,tidx,i],output[fidx,tidx,i] = transmission3d(time_range[tidx],freq_range[fidx],atom_arr,M_wg,M_atom)
            end
        end
    end
    return sdata(output), sdata(output_matrix)
end

function generate_atom_array3d(atom_num,total_atom_num,lattice_scale)
    #lattice_scale = lattice_width/lattice_unit
    #(traj id,lattice id)*atom_num
    atom_array = zeros(Int64,(2,atom_num))
    atom_array[1,:] = sample(collect(1:total_atom_num),atom_num)
    #    atom_array[1,:] = shuffle(collect(1:total_atom_num))[1:atom_num] #atom id
    d = Truncated(Normal(lattice_scale/2, atom_beam_waist/lattice_unit), 0, lattice_scale)
    atom_pos = round.(Int64,rand(d,atom_num))
    @assert all(lattice_scale .>= atom_pos .>= 0)
    atom_array[2,:] = sort(atom_pos) #atom's position in unit of lattice width
#    atom_array[2,:] = sort(round(Int64,randn(atom_num)*lattice_scale)) #atom's position in unit of lattice width
    return atom_array
end

function calc_gamma1d_3d(pos,t)
    @assert any(isnan.(pos)) == false
    g1d = Float64(Fields.value(pos[1:3],t,Probe::ScalarFieldNode)*gamma_1d)
    @assert isnan(g1d) == false
    return g1d
end

function transmission3d(t::Float64,detune::Float64,atom_arr::Array{Int64,2},M_wg::Array{Complex{Float64},3},M_atom::Array{Complex{Float64},3})
    # calculate transmission at time t and detuning detune.
    # put wg from -2*Lattice_unit to first atom, and after last atom to lattice_width+2*lattice_unit
    atom_num::Int64 = size(atom_arr,2)
    x_point_k::Float64 = pi/lattice_unit::Float64
    k::Float64 = x_point_k*k_ratio::Float64
    #generate waveguide transfer matrix

    atom_pos = atom_arr[2,:]*lattice_unit
    for i = 1:atom_num
        pos::Vector{Float64} = Trajs[t,atom_arr[1,i]]
        if any(isnan.(pos[1:3]))
            wg_transfer_matrix(M_atom,i,0.0,0.0)
        else
            atom_pos[i] += pos[3]            
            g1d = calc_gamma1d_3d(pos[1:3],t)
            f_0 = Fields.value(pos[1:3],t,ForceFields::ScalarFieldNode)*(-2.08e4) #*20.8/(-1e-3)
            #vector shift
            if vector_shift == 1
                mf = sample(-3:3)
                gf = -0.25
                f_0 *= 1+0.5*gf*mf
            end
            atom_transfer_matrix(M_atom,i,detune,f_0,g1d,gamma_prime::Float64)
        end
    end
    
    @assert all(lattice_width+3*lattice_unit .> atom_pos .> -3*lattice_unit)
    ldiff::Vector{Float64} = diff(atom_pos)
    
    for i = 1:length(ldiff)
        wg_transfer_matrix(M_wg,i,k,ldiff[i])
    end
    first_wg = wg_transfer_matrix(k,atom_pos[1]-(-2lattice_unit))
    last_wg = wg_transfer_matrix(k,lattice_width+2lattice_unit-atom_pos[end])
    M_tot::Array{Complex{Float64},2} = first_wg*M_atom[:,:,1];
    @fastmath for i = 1:atom_num - 1
        @inbounds M_tot *= @view M_atom[:,:,i+1]
        @inbounds M_tot *= @view M_wg[:,:,i]
    end
    M_tot *= last_wg
    return M_tot,one(Complex{Float64})/M_tot[2,2]
end
