using PyCall
using HDF5
function single_scan_scaling(config::Dict,sfn::ScalarFieldNode,output_file)
    range = config["range"]
    field_name = config["field"]
    scaling = config["scaling"]
    for i = 1:range
        s = replace(scaling,"@i",float(i))
        println("change scaling of field $field_name to ",s)
        s_exp = eval(parse(s))
        Fields.setscaling!(Fields.find_field(x->x.name==ascii(field_name),sfn),s_exp)
        Fields.init_parallel!(sfn)
        println("start calculation...")
        result = calculate()
        println("save results...")
        matwrite(output_file*string(i)*".mat",Dict(
                                                   "result"=>result,
                                                   "tspan" =>TrajSolver.get_tspan(),
                                                   "pos"=>sfn.position,
                                                   "siz"=>sfn.size
                                                   ))
        h5write(output_file*string(i)*".h5", "result", result)
        output_image(sfn,0.0,[69200.0, 69650.0, 100.0, 49900.0],output_file*string(i)*".png")
        output = Fields.composite_slow([69200.0, 69650.0, 100.0, 49900.0],0.0)
        savemat(output_file*string(i)*"_usmall.mat",output,"output")
        calc_score(result)
#        println("outputing movie...")
#        output_movie(sfn,collect(0:0.1:10),[5000, 20000.0, 20000.0, 30000.0],output_file*string(i)*"_u.mp4")
    end
end

function output_movie(sfn,tspan,range,filename)
    try:
        rm("/tmp/movie",recursive=true)
    end
    mkdir("/tmp/movie")
    output_0 = Fields.composite_slow(range,0.0)
    v_min = minimum(output_0)
    v_max = maximum(output_0)
    #TODO: parallel image output
    for t in enumerate(tspan)
        println(t)
        output_image(sfn,t[2],range,"/tmp/movie/img"*@sprintf("%04d",t[1])*".png",v_min=v_min,v_max=v_max)
    end
    hash_key = string(hash(time()))

    current_folder = pwd()
    movie_folder = "/tmp/movie"*hash_key
    cd(movie_folder)
    run(`ffmpeg -framerate 10 -i img%04d.png -s:v 1280x720 -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p out.mp4`)
    cd(current_folder)
    cp("/tmp/movie/out.mp4",filename,remove_destination=true)
    rm(movie_folder,recursive=true)
end

function traj_plot(sfn,result,tspan,range,filename)
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

function output_image_gp(sfn,t,range,filename;v_min=0.0,v_max=0.0)

end

function output_image(sfn,t,range,filename;v_min=0.0,v_max=0.0)#range x1 x2 y1 y2
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
    pp = Polygon([9890.0 10110.0 10110.0 9890.0;25338.0 25338.0 24662 24662])
    score = 0
    for i = 1:size(result,3)
        for j = 1:size(result,2)
            if pointInPolygon(pp,result[1:2,j,i])
                score+=1
            end
        end
    end
    println("score:",score)
end
