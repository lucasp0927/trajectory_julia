include("TrajAnalyzer_transfermatrix.jl")

function calculate_transmission(paras::Para)
    # traj selection
    TrajAnalysis.initialize(paras.phy_para["beams-geo"])
    traj_selected = select_traj(paras.traj,paras.t_span,paras.select)
    # calculate transmission
    init_parallel(paras.phy_para,paras.sim_para,paras.time_arr,paras.freq_arr,paras.t_span,traj_selected,paras.probe)
    #preallocate matrix
    output = Array(Complex{Float64},(length(paras.freq_arr),length(paras.time_arr),paras.sim_para["iteration"]))
    output_matrix = Array(Complex{Float64},(2,2,length(paras.freq_arr),length(paras.time_arr),paras.sim_para["iteration"]))
    parallel_transmission!(output,output_matrix)
    return output, output_matrix
    #output_data(output,output_matrix,paras.outputfile)
end

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
