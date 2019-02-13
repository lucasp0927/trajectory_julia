using MAT
using HDF5
using ArgParse
function dicttoh5(filename, data)
    @assert filename[end-2:end] == ".h5"
    varnames = collect(keys(data))
    fid = h5open(filename,"w")
    for name in varnames
        @info name
        @info typeof(data[name])
        write_dataset(fid,name,data[name])
    end
    close(fid)
end

function write_dataset(fid,name,data::T) where T<:Number
    @info name
    write(fid,name,data)
end

function write_dataset(fid,name,data::Array{T}) where T<:Real
    @info name
    #write(fid,name,data)
    fid[name,"shuffle",(),"compress",5] = data    
end

function write_dataset(fid,name,data::Array{T}) where T<:Complex
    @info name
    #write(fid,name*"/real",real(data))
    #write(fid,name*"/imag",imag(data))
    g = g_create(fid,name)
    g["real","shuffle",(),"compress",5] = real(data)
    g["imag","shuffle",(),"compress",5] = imag(data)    
end

function h5todict(filename)
    #read variables under /
    @assert filename[end-2:end] == ".h5"
    data = Dict()
    fid = h5open(filename,"r")
    varnames = names(fid)
    for name in varnames
        if isa(fid[name],HDF5Group)
            data[name] = read(fid,name*"/real") .+ (read(fid,name*"/imag").*1im)
        elseif isa(fid[name],HDF5Dataset)
            data[name] = read(fid,name)
        else
            @error "unexpected HDF5 file structure"
        end
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

function check_data(data1,data2)
    @info "Check data consistency"
    if isequal(data1,data2)
        @info "pass"
    else
        @error "data error"
    end
    # keys1 = keys(data1)
    # keys2 = keys(data2)
    # @info keys1
    # @info keys2
    # @assert keys1 == keys2

    # for key in keys1
    #     @info "check "*key
    #     @info size(data1[key])
    #     @info size(data2[key])        
    #     @assert isequal(data1[key],data2[key])
    # end
end

function main()
    parsed_args = parse_commandline()
    # println("Parsed args:")
    # for (arg,val) in parsed_args
    #     println("  $arg  =>  $val")
    # end
    infile = parsed_args["infile"]
    @assert occursin(".mat", infile)
    outfile = parsed_args["outfile"]
    @assert occursin(".h5", outfile)    
    data = matread(infile)
    dicttoh5(outfile,data)
    # check if data consistency
    data2 = h5todict(outfile)
    check_data(data,data2)
end

main()

