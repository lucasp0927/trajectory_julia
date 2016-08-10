using PyCall
include("output.jl")
include("flux.jl")
function single_scan_scaling(trajsolver_config::Dict,config::Dict,sfn::ScalarFieldNode,output_file,calc_traj_flag::Bool,spectrum_flag::Bool,movie_flag::Bool)
    range = config["range"]
    field_name = config["field"]
    scaling = config["scaling"]
    score = [ascii(k)=>zeros(Int64,range) for k in keys(config["score"])]
    for i = 1:range
        s = replace(scaling,"@i",float(i))
        Lumberjack.info("change scaling of field $field_name to ",s)
        s_exp = eval(parse(s))
        Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        Fields.init_parallel!(sfn)
        # build probe beam
        if calc_traj_flag
            result = calculate_traj(i)
            Lumberjack.info("save results...")
            matwrite(output_file*string(i)*".mat",result)
            traj = result["traj"]
            tspan = result["tspan"]
        else
            Lumberjack.info("read results...")
            result = matread(output_file*string(i)*".mat")
            traj = result["traj"]
            tspan = result["tspan"]
        end
        if spectrum_flag || movie_flag
            Lumberjack.info("Initialize TrajAnalyzer...")
            probe_sfn = Fields.buildAndAlign(config["probe"]["field"],0,name=ascii([k for k in keys(config["probe"])][1]))
            TrajAnalyzer.init_parallel!(result,probe_sfn,sfn,config)
        end
        if spectrum_flag
            TrajAnalyzer.spectrum(output_file*string(i)*"_tm")
        end
        if movie_flag
            movie_range = [promote(config["movie-output"]["range"]...)...]
            TrajAnalyzer.output_movie_traj(config["movie-output"],output_file*string(i)*"_traj.mp4")
        end
        #score and flux
        # for (k,v) in config["score"]
        #     Lumberjack.info("calculating score for area $k...")
        #     @time (score[ascii(k)])[i] = TrajAnalyzer.calc_score(v)
        # end
        #@time flux = calc_flux(traj,tspan,config["flux"],output_file*string(i)*"_flux.mat")
        #output intial range potential
        #init_range = get_large_init_range(values(trajsolver_config["atom-config"]["init-range"]))
        #t0 = trajsolver_config["simulation-config"]["tstart"]
        #        TrajAnalyzer.output_image_gp(t0,init_range,output_file*string(i)*"_init_range.png",save_data = true, data_filename=output_file*string(i)*"_init_range.h5")
        # println("output 0")
        # TrajAnalyzer.output_image_gp(0.0,[10000.0-750.0,10750.0,25000.0-750.0,25000.0+750.0],output_file*string(i)*"_phase_potential_0.png",save_data = true, data_filename=output_file*string(i)*"_phase_potential_0.h5")
#        TrajAnalyzer.output_image_gp_traj(t0,init_range,10.0,10.0,output_file*string(i)*"_init_range_traj.png")
#        TrajAnalyzer.output_image_gp(0.0,movie_range,output_file*string(i)*"_probe.png",TrajAnalyzer.Probe)
            # output individual frames
            # config = config["movie-output"]
            # frame_range = [promote(config["range"]...)...]
            # res = config["res"]
            # res_x = res[1]
            # res_y = res[2]
            # TrajAnalyzer.output_image_gp_traj(171.5,frame_range,res_x,res_y,output_file*"frame_1.png")
            # TrajAnalyzer.output_image_gp_traj(172.4,frame_range,res_x,res_y,output_file*"frame_2.png")
            # TrajAnalyzer.output_image_gp_traj(173.1,frame_range,res_x,res_y,output_file*"frame_3.png")
            # TrajAnalyzer.output_image_gp_traj(173.9,frame_range,res_x,res_y,output_file*"frame_4.png")
            # TrajAnalyzer.output_image_gp_traj(174.3,frame_range,res_x,res_y,output_file*"frame_5.png")
    end
 #   matwrite(output_file*"score.mat",score)
end

function get_large_init_range(init_range)
    #find the largest area that contain the whole init_range
    init_range = convert(Array{Float64},cat(2,init_range...))
    init_range_large = zeros(Float64,4)
    init_range_large[1] = minimum(init_range[1,:])
    init_range_large[2] = maximum(init_range[2,:])
    init_range_large[3] = minimum(init_range[3,:])
    init_range_large[4] = maximum(init_range[4,:])
    return init_range_large
end
