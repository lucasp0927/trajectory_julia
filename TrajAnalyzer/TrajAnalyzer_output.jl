function output_movie_traj(config,filename)
    mov_tspan = collect(config["tstart"]:config["tdiv"]:config["tend"])
    mov_range = [promote(config["range"]...)...]
    mov_res = config["res"]
    output_movie(mov_tspan,mov_range,mov_res,filename,traj=true)
end

function output_movie_traj_flux(config,filename,result,tspan,flux)
    flux_arr = reduce((x,y)->cat(2,x,y),values(flux))
    #TODO: sort by time
    mov_tspan = collect(config["tstart"]:config["tdiv"]:config["tend"])
    mov_range = [promote(config["range"]...)...]
    mov_res = config["res"]
    output_movie(mov_tspan,mov_range,mov_res,filename,traj=true,result=result,tspan=tspan)
end

function output_movie(mov_tspan,range,res,filename;traj=false)
    output_0 = Fields.composite_slow(range,0.0)
    v_min = minimum(output_0)
    v_max = maximum(output_0)
    hash_key = string(hash(rand()))
    current_folder = pwd()
    movie_folder = "/tmp/movie"*hash_key
    mkdir(movie_folder)
    if traj==false
        @sync @parallel for t in collect(enumerate(mov_tspan))
            output_image_gp(t[2],range,movie_folder*"/img"*@sprintf("%04d",t[1])*".png",v_min=v_min,v_max=v_max)
        end
    else
        res_x = res[1]
        res_y = res[2]
        tdiv = mean(diff(mov_tspan))
        @sync @parallel for t in collect(enumerate(mov_tspan))
            output_image_gp_traj(t[2],range,res_x,res_y,movie_folder*"/img"*@sprintf("%04d",t[1])*".png",v_min=v_min,v_max=v_max,tdiv=tdiv)
        end
    end
    cd(movie_folder)
    run(pipeline(`ffmpeg -framerate 5 -i img%04d.png -s:v 1300x1000 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4`,stderr=current_folder*"/ffmpeg.log"))
    cd(current_folder)
    cp(movie_folder*"/out.mp4",filename,remove_destination=true)
    rm(movie_folder,recursive=true)
end


function output_image_gp(t,range,filename;v_min=0.0,v_max=0.0)
    output_data = Fields.composite_slow_with_position(range,t,[20.0,20.0])
    current_folder = pwd()
    Lumberjack.debug("current_folder: $current_folder")
    hash_key = string(hash(rand()))
    image_folder = "/tmp/image"*hash_key
    mkdir(image_folder)
    cd(image_folder)
    h5write(image_folder*"/data.h5", "output", output_data)
    run(`h5totxt data.h5 -o data.txt`)
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


function output_image_gp_traj(t,range,res_x,res_y,filename;v_min=0.0,v_max=0.0,tdiv=0.0)
    #find closest index
    t_idx=indmin(abs(Trajs.tspan-t))
    tmp = squeeze(Trajs.traj[1:2,t_idx,:],2)
    dots = zeros(Float64,3,size(tmp,2))
    dots[1:2,:] = tmp[:,:]
    #make atom red before vanishing.
    if tdiv>0.0
        t_idx_nx=indmin(abs(Trajs.tspan-(t+tdiv)))
        tmp_next = squeeze(Trajs.traj[1:2,t_idx_nx,:],2)
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
    run(`h5totxt data.h5 -o data.txt`)
    run(`h5totxt dots.h5 -o dots.txt`)
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
