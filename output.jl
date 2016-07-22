function output_movie_traj(config,filename,result,tspan)
    mov_tspan = collect(config["tstart"]:config["tdiv"]:config["tend"])
    mov_range = [promote(config["range"]...)...]
    mov_res = config["res"]
    output_movie(mov_tspan,mov_range,mov_res,filename,traj=true,result=result,tspan=tspan)
end

function output_movie_traj_flux(config,filename,result,tspan,flux)
    flux_arr = reduce((x,y)->cat(2,x,y),values(flux))
    #TODO: sort by time
    mov_tspan = collect(config["tstart"]:config["tdiv"]:config["tend"])
    mov_range = [promote(config["range"]...)...]
    mov_res = config["res"]
    output_movie(mov_tspan,mov_range,mov_res,filename,traj=true,result=result,tspan=tspan)
end

function output_movie(output_tspan,range,res,filename;traj=false,result=[],tspan=[])
    output_0 = Fields.composite_slow(range,0.0)
    v_min = minimum(output_0)
    v_max = maximum(output_0)
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
        tdiv = mean(diff(output_tspan))
        @sync @parallel for t in collect(enumerate(output_tspan))
            output_image_gp_traj(t[2],range,result_s,tspan_s,res_x,res_y,movie_folder*"/img"*@sprintf("%04d",t[1])*".png",v_min=v_min,v_max=v_max,tdiv=tdiv)
        end
    end
    cd(movie_folder)
    run(pipeline(`ffmpeg -framerate 5 -i img%04d.png -s:v 1300x1000 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4`,stderr=current_folder*"/ffmpeg.log"))
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

@everywhere function output_image_gp_traj(t,range,result,tspan,res_x,res_y,filename;v_min=0.0,v_max=0.0,tdiv=0.0)
    #find closest index
    t_idx=indmin(abs(tspan-t))
    tmp = squeeze(result[1:2,t_idx,:],2)
    dots = zeros(Float64,3,size(tmp,2))
    dots[1:2,:] = tmp[:,:]
    #make atom red before vanishing.
    if tdiv>0.0
        t_idx_nx=indmin(abs(tspan-(t+tdiv)))
        tmp_next = squeeze(result[1:2,t_idx_nx,:],2)
        dots[3,:] = [(isnan(tmp[1,i])==false && isnan(tmp_next[1,i])==true)?1.0:0.0 for i in 1:size(tmp,2)]
    end
    output_data = Fields.composite_slow_with_position(range,t,[res_x,res_y])
    current_folder = pwd()
    hash_key = string(hash(rand()))
    image_folder = "/tmp/image"*hash_key
    mkdir(image_folder)
    cd(image_folder)
    h5write(image_folder*"/data.h5", "output", output_data)
    h5write(image_folder*"/dots.h5", "output", dots)
    run(`h5totxt -o data.txt data.h5`)
    run(`h5totxt -o dots.txt data.h5`)
    cp(current_folder*"/gnuplot/output_image_gp_traj.gp",image_folder*"/output_image_gp_traj.gp")
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
    run(`h5totxt -o data.txt data.h5`)
    cp(current_folder*"/gnuplot/output_image_gp.gp",image_folder*"/output_image_gp.gp")
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
