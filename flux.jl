#using Iterators
# Given three colinear points p, q, r, the function checks if
# point q lies on line segment 'pr'
@everywhere @inbounds function onSegment(p::Vector{Float64},q::Vector{Float64},r::Vector{Float64})
    return  (q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1]) && q[2] <= max(p[2], r[2]) && q[2] >= min(p[2], r[2]))
end

# To find orientation of ordered triplet (p, q, r).
# The function returns following values
# 0 --> p, q and r are colinear
# 1 --> Clockwise
# 2 --> Counterclockwise
@everywhere @inbounds @fastmath function orientation(p::Vector{Float64},q::Vector{Float64},r::Vector{Float64})
    #See http://www.geeksforgeeks.org/orientation-3-ordered-points/
    val = (q[2] - p[2]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[2] - q[2])
    if val==0
        return 0
    else
        return val>0?1:2
    end
end

# The main function that returns true if line segment 'p1q1'
# and 'p2q2' intersect.
@everywhere function doIntersect(p1::Vector{Float64},q1::Vector{Float64},p2::Vector{Float64},q2::Vector{Float64})
    #return (b,d) tuple, b = true if p1q1 and p2q2 intersect,
    #d = 1 if p2 is on the left side of p1q1 line.
    #d = -1 if p2 is on the right side of p1q1 line.
    #else d = 0
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)
    if (o1 != o2 && o3 != o4)
        if (o1==2 && o2 ==1)
            return (true,1.0)
        else
            return (true,-1.0)
        end
    elseif (o1 == 0 && onSegment(p1, p2, q1))
        return (true,0.0)
    elseif (o2 == 0 && onSegment(p1, q2, q1))
        return (true,0.0)
    elseif (o3 == 0 && onSegment(p2, p1, q2))
        return (true,0.0)
    elseif (o4 == 0 && onSegment(p2, q1, q2))
        return (true,0.0)
    end
    return (false,0.0)
end

@everywhere function cat_ignore_empty(x::Array{Float64},y::Array{Float64})
    if isempty(x) && isempty(y)
        return Float64[]
    elseif isempty(x)
        return cat(2,y)
    elseif isempty(y)
        return cat(2,x)
    else
        return cat(2,x,y)
    end
end

#TODO optimize this
@inbounds function calc_flux2(traj,tspan,config,filename)
    flux = Dict{Any,Any}([ascii(k)=>[] for k in keys(config)])
    tmp = zeros(Float64,4)
    traj_s = copy_to_sharedarray!(traj)
    for (k,v) in config
        v = [promote(v...)...]
        @assert length(v)==4 "wrong array length"
        k = ascii(k)
        flux[k] = @parallel ((x,y)->cat_ignore_empty(x,y)) for i = 1:size(traj_s,3) #loop over trajectories
            tmp = map(1:length(tspan)-1)do j
                (b,d) = doIntersect(v[1:2],v[3:4],traj_s[1:2,j,i],traj_s[1:2,j+1,i])
                b?Float64[tspan[j],traj_s[1:4,j,i]...,d]:Float64[]
            end
            reduce((x,y)->cat_ignore_empty(x,y),tmp)
        end
    end
    matwrite(filename,flux)
    for k in keys(flux)
        println(k," ",size(flux[k]))
    end
    return flux
end

@inbounds function calc_flux(traj,tspan,config,filename)
    flux = Dict{Any,Any}([ascii(k)=>[] for k in keys(config)])
    tmp = zeros(Float64,4)
    traj_s = copy_to_sharedarray!(traj)
    tspan_s = copy_to_sharedarray!(tspan)
    for (k,v) in config
        v = [promote(v...)...]
        @assert length(v)==4 "wrong array length"
        k = ascii(k)
        flux[k] = @parallel ((x,y)->cat_ignore_empty(x,y)) for i = 1:size(traj_s,3) #loop over trajectories
            tmp = zeros(Float64,7)
            for j = 1:length(tspan_s)-1
            (b,d) = doIntersect(v[1:2],v[3:4],traj_s[1:2,j,i],traj_s[1:2,j+1,i])
                if b
                    tmp = cat(2,tmp,Float64[Float64(i),tspan_s[j],traj_s[1:4,j,i]...,d])
                end
            end
            tmp[:,2:end]
        end
    end
    matwrite(filename,flux)
    for k in keys(flux)
        println(k," ",size(flux[k]))
    end
    return flux
end
