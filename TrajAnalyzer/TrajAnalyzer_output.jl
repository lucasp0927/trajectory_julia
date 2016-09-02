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
    output_0 = Fields.composite(range,0.0)
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
    if contains(platform.platform(),"Ubuntu-12.04")||contains(platform.platform(),"Ubuntu-14.04")
        run(pipeline(`avconv -framerate 5 -i img%04d.png -s:v 1300x1000 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4`,stderr=current_folder*"/ffmpeg.log"))
    else
        run(pipeline(`ffmpeg -framerate 5 -i img%04d.png -s:v 1300x1000 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4`,stderr=current_folder*"/ffmpeg.log"))
    end
    cd(current_folder)
    cp(movie_folder*"/out.mp4",filename,remove_destination=true)
    rm(movie_folder,recursive=true)
end

function output_image_gp(t,range,filename,sfn;v_min=0.0,v_max=0.0,xres=1300,yres=1000,save_data=false,data_filename = "")
    Lumberjack.debug("in TrajAnalyzer.output_image_gp()... with sfn")
    output_data = Fields.composite_with_position(range,t,[10.0,10.0],sfn)
    current_folder = pwd()
    Lumberjack.debug("current_folder: $current_folder")
    hash_key = string(hash(rand()))
    image_folder = "/tmp/image"*hash_key
    mkdir(image_folder)
    cd(image_folder)
    h5write(image_folder*"/data.h5", "output", output_data)
    run(`h5totxt -o data.txt data.h5`)
    cp(current_folder*"/gnuplot/output_image_gp.gp",image_folder*"/output_image_gp.gp")
    if v_min==0.0 && v_max==0.0
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);xres=$xres;yres=$yres" output_image_gp.gp`)
    else
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);xres=$xres;yres=$yres;set cbrange [$v_min:$v_max]" output_image_gp.gp`)
    end
    cd(current_folder)
    cp(image_folder*"/data.png",filename,remove_destination=true)
    if save_data
        cp(image_folder*"/data.h5",data_filename,remove_destination=true)
    end
    rm(image_folder,recursive=true)
end

function output_image_gp(t,range,filename;v_min=0.0,v_max=0.0,xres=1300,yres=1000,save_data=false,data_filename = "")
    Lumberjack.debug("in TrajAnalyzer.output_image_gp()")
    output_data = Fields.composite_with_position(range,t,[10.0,10.0])
    current_folder = pwd()
    Lumberjack.debug("current_folder: $current_folder")
    hash_key = string(hash(rand()))
    image_folder = "/tmp/image"*hash_key
    mkdir(image_folder)
    cd(image_folder)
    h5write(image_folder*"/data.h5", "output", output_data)
    run(`h5totxt -o data.txt data.h5`)
    cp(current_folder*"/gnuplot/output_image_gp.gp",image_folder*"/output_image_gp.gp")
    if v_min==0.0 && v_max==0.0
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);xres=$xres;yres=$yres" output_image_gp.gp`)
    else
        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);xres=$xres;yres=$yres;set cbrange [$v_min:$v_max]" output_image_gp.gp`)
    end
    cd(current_folder)
    cp(image_folder*"/data.png",filename,remove_destination=true)
    if save_data
        cp(image_folder*"/data.h5",data_filename,remove_destination=true)
    end
    rm(image_folder,recursive=true)
end

function output_image_gp_traj(t,range,res_x,res_y,filename;v_min=0.0,v_max=0.1,tdiv=0.0)
    if (v_max-v_min <= eps())
        v_max = v_max + 0.1
    end
    tmp = (Trajs[t,:])[1:2,:]
    dots = zeros(Float64,3,size(tmp,2))
    dots[1:2,:] = tmp[:,:]
    #make atom red before vanishing.
    if tdiv>0.0
        tmp_next = (Trajs[t+tdiv,:])[1:2,:]
        dots[3,:] = [(isnan(tmp[1,i])==false && isnan(tmp_next[1,i])==true)?1.0:0.0 for i in 1:size(tmp,2)]
    end
    output_data = Fields.composite_with_position(range,t,[res_x,res_y])
    current_folder = pwd()
    hash_key = string(hash(rand()))
    image_folder = "/tmp/image"*hash_key
    mkdir(image_folder)
    cd(image_folder)
    h5write(image_folder*"/data.h5", "output", output_data)
    h5write(image_folder*"/dots.h5", "output", dots)
    run(`h5totxt -o data.txt data.h5`)
    run(`h5totxt -o dots.txt dots.h5`)
    cp(current_folder*"/gnuplot/output_image_gp_traj.gp",image_folder*"/output_image_gp_traj.gp")
    if v_min==0.0 && v_max==0.0
        # stupid hack to fix v_min and v_max
        v_min = minimum(output_data[3,:,:])
        v_max = maximum(output_data[3,:,:])
#        run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);time=$t" output_image_gp_traj.gp`)
    end
    run(`gnuplot -e "xstart=$(range[1]);xend=$(range[2]);ystart=$(range[3]);yend=$(range[4]);time=$t;set cbrange [$v_min:$v_max]" output_image_gp_traj.gp`)
    cd(current_folder)
    cp(image_folder*"/data.png",filename,remove_destination=true)
    rm(image_folder,recursive=true)
end
