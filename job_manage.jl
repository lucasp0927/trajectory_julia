using PyCall
include("output.jl")
include("flux.jl")
function single_scan_scaling(trajsolver_config::Dict,config::Dict,sfn::ScalarFieldNode,output_file,calc_traj_flag::Bool,movie_flag::Bool)
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
        probe_sfn = Fields.buildAndAlign(config["probe"]["field"],0,name=ascii([k for k in keys(config["probe"])][1]))
        if calc_traj_flag
            result = calculate_traj()
            TrajAnalyzer.init_parallel!(result,probe_sfn,sfn,config)
            Lumberjack.info("save results...")
            matwrite(output_file*string(i)*".mat",result)
            traj = result["traj"]
            tspan = result["tspan"]
        else
            Lumberjack.info("read results...")
            result = matread(output_file*string(i)*".mat")
            TrajAnalyzer.init_parallel!(result,probe_sfn,sfn,config)
            traj = result["traj"]
            tspan = result["tspan"]
        end
        # for (k,v) in config["score"]
        #     Lumberjack.info("calculating score for area $k...")
        #     @time (score[ascii(k)])[i] = TrajAnalyzer.calc_score(v)
        # end

        #output intial range potential
        init_range = [promote(trajsolver_config["atom-config"]["init-range"]...)...]
        t0 = trajsolver_config["simulation-config"]["tstart"]
        TrajAnalyzer.output_image_gp(t0,init_range,output_file*string(i)*"_init_range.png")

        Lumberjack.debug("testing calculate_transmission()")
        output,output_matrix = TrajAnalyzer.calculate_transmission()
        specname = output_file*string(i)*"_spectrum.mat"
        spec = Dict(
                    "output"=>output,
                    "output_matrix"=>output_matrix
                    )
        Lumberjack.debug("save output matrix go $specname")
        matwrite(specname,spec)

        movie_range = [promote(config["movie-output"]["range"]...)...]
        TrajAnalyzer.output_image_gp(0.0,movie_range,output_file*string(i)*"_probe.png",TrajAnalyzer.Probe)
#        @time flux = calc_flux(traj,tspan,config["flux"],output_file*string(i)*"_flux.mat")
        if movie_flag
            @time TrajAnalyzer.output_movie_traj(config["movie-output"],output_file*string(i)*"_traj.mp4")
        end
    end
    matwrite(output_file*"score.mat",score)
end