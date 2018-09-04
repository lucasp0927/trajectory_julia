using MAT
using HDF5
using ArgParse
function dicttoh5(filename, data)
    @assert filename[end-2:end] == ".h5"
    varnames = collect(keys(data))
    fid = h5open(filename,"w")
    for name in varnames
        write(fid,name,data[name])
    end
    close(fid)
end

function h5todict(filename)
    #read variables under /
    @assert filename[end-2:end] == ".h5"
    data = Dict()
    fid = h5open(filename,"r")
    varnames = names(fid)
    for name in varnames
        data[name] = read(fid,name)
    end
    close(fid)
    return data
end

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "infile"
            help = "input files"
            required = true
        "outfile"
            help = "output files"
            required = true        
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    # println("Parsed args:")
    # for (arg,val) in parsed_args
    #     println("  $arg  =>  $val")
    # end
    infile = parsed_args["infile"]
    outfile = parsed_args["outfile"]    
    data = matread(infile)
    dicttoh5(outfile,data)
    data2 = h5todict(outfile)
    @assert data == data2
end

main()

