using PyCall
@everywhere using HDF5
function single_scan_scaling(config::Dict,sfn::ScalarFieldNode,output_file)
    range = config["range"]
    field_name = config["field"]
    scaling = config["scaling"]
    score = zeros(Int64,range)
    for i = 1:range
        s = replace(scaling,"@i",float(i))
        println("change scaling of field $field_name to ",s)
        s_exp = eval(parse(s))
        Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        Fields.init_parallel!(sfn)
        println("start calculation...")
        result = calculate()
        println("save results...")
        tspan = TrajSolver.get_tspan()
        matwrite(output_file*string(i)*".mat",Dict(
                                                   "result"=>result,
                                                   "tspan" =>tspan,
                                                   "pos"=>sfn.position,
                                                   "siz"=>sfn.size
                                                   ))
        score[i] = calc_score(result)
        @time output_movie_traj(config["movie-output"],output_file*string(i)*"_traj.mp4",result,tspan)
    end
    #save score and plot score
    savemat(output_file*"score.mat",score,"score")
    pyimport("matplotlib")[:use]("Agg")
    @pyimport matplotlib.pyplot as plt
    prange = config["plot-range"]
    xx = collect(linspace(prange...,range))
    plt.plot(xx,score)
    plt.savefig(output_file*"score.png")
    plt.clf()
end

function output_movie_traj(config,filename,result,tspan)
        mov_tspan = collect(config["tstart"]:config["tdiv"]:config["tend"])
        mov_range = [promote(config["range"]...)...]
        mov_res = config["res"]
        output_movie(mov_tspan,mov_range,mov_res,filename,traj=true,result=result,tspan=tspan)
end

function output_movie(output_tspan,range,res,filename;traj=false,result=[],tspan=[])
    output_0 = Fields.composite_slow(range,0.0)
    v_min = minimum(output_0)
    v_max = maximum(output_0)
    #TODO: parallel image output
    hash_key = string(hash(rand()))
    current_folder = pwd()
    movie_folder = "/tmp/movie"*hash_key
    mkdir(movie_folder)
    if traj==false
        @sync @parallel for t in collect(enumerate(output_tspan))
            output_image_gp(t[2],range,movie_folder*"/img"*@sprintf("%04d",t[1])*".png",v_min=v_min,v_max=v_max)
        end
    else
        res_x = res[1]
        res_y = res[2]
        result_s=copy_to_sharedarray!(result)
        tspan_s=copy_to_sharedarray!(tspan)
        @sync @parallel for t in collect(enumerate(output_tspan))
            output_image_gp_traj(t[2],range,result_s,tspan_s,res_x,res_y,movie_folder*"/img"*@sprintf("%04d",t[1])*".png",v_min=v_min,v_max=v_max)
        end
    end
    cd(movie_folder)
    run(pipeline(`ffmpeg -framerate 10 -i img%04d.png -s:v 1300x1000 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4`,stderr=current_folder*"/ffmpeg.log"))
    cd(current_folder)
    cp(movie_folder*"/out.mp4",filename,remove_destination=true)
    rm(movie_folder,recursive=true)
end

function traj_plot(result,tspan,range,filename)
    pyimport("matplotlib")[:use]("Agg")
    @pyimport matplotlib.pyplot as plt

    for i = 1:size(result,3)
        println(i)
        plt.plot(result[1,1:5:end,i],result[2,1:5:end,i],"r,")
    end
    plt.xlim(range[1:2])
    plt.ylim(range[3:4])
    plt.savefig(filename)
    plt.clf()
end

@everywhere function output_image_gp_traj(t,range,result,tspan,res_x,res_y,filename;v_min=0.0,v_max=0.0)
    #find closest index
    t_idx=indmin(abs(tspan-t))
    tmp = squeeze(result[1:2,t_idx,:],2)
    dots = zeros(Float64,3,size(tmp,2))
    dots[1:2,:] = tmp[:,:]
    output_data = Fields.composite_slow_with_position(range,t,[res_x,res_y])
    current_folder = pwd()
    hash_key = string(hash(rand()))
    image_folder = "/tmp/image"*hash_key
    mkdir(image_folder)
    cd(image_folder)
    h5write(image_folder*"/data.h5", "output", output_data)
    h5write(image_folder*"/dots.h5", "output", dots)
    run(`h5totxt data.h5 -o data.txt`)
    run(`h5totxt dots.h5 -o dots.txt`)
    cp(current_folder*"/output_image_gp_traj.gp",image_folder*"/output_image_gp_traj.gp")
    if v_min==0.0 && v_max==0.0
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);time=$t" output_image_gp_traj.gp`)
    else
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);time=$t;set cbrange [$v_min:$v_max]" output_image_gp_traj.gp`)
    end
    cd(current_folder)
    cp(image_folder*"/data.png",filename,remove_destination=true)
    rm(image_folder,recursive=true)
end

@everywhere function output_image_gp(t,range,filename;v_min=0.0,v_max=0.0)
    output_data = Fields.composite_slow_with_position(range,t,[20.0,20.0])
    current_folder = pwd()
    hash_key = string(hash(rand()))
    image_folder = "/tmp/image"*hash_key
    mkdir(image_folder)
    cd(image_folder)
    h5write(image_folder*"/data.h5", "output", output_data)
    run(`h5totxt data.h5 -o data.txt`)
    cp(current_folder*"/output_image_gp.gp",image_folder*"/output_image_gp.gp")
    if v_min==0.0 && v_max==0.0
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4])" output_image_gp.gp`)
    else
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);set cbrange [$v_min:$v_max]" output_image_gp.gp`)
    end
    cd(current_folder)
    cp(image_folder*"/data.png",filename,remove_destination=true)
    rm(image_folder,recursive=true)
end

function output_image(t,range,filename;v_min=0.0,v_max=0.0)#range x1 x2 y1 y2
    pyimport("matplotlib")[:use]("Agg")
    @pyimport matplotlib.pyplot as plt
    output = Fields.composite_slow(range,t)
    if v_min==0.0 && v_max==0.0
        plt.imshow(output, extent=[range[3],range[4],range[2],range[1]])
    else
        plt.imshow(output, vmin=v_min, vmax=v_max, extent=[range[3],range[4],range[2],range[1]])
    end
    plt.colorbar()
    plt.savefig(filename)
    plt.clf()
end

function calculate()
    @time    @sync begin
        for p = 2:nprocs()
            @async remotecall_wait(p,TrajSolver.solve_traj)
        end
    end
    temp = cell(nworkers())
    for p = 2:nprocs()
        temp[p-1] = remotecall_fetch(p,TrajSolver.get_result)
    end
    result = cat(3,temp...)
    return result
end

function calc_score(result)
    println("calculate score...")
    include("./TrajSolver/polygon.jl")
    pp = Polygon([9890.0,10110.0,10110.0,9890.0,25338.0,25338.0,24662,24662])
    score = 0
    for i = 1:size(result,3)
        for j = 1:size(result,2)
            if pointInPolygon(pp,result[1:2,j,i])
                score+=1
            end
        end
    end
    println("score:",score)
    return score
end
