# Given three colinear points p, q, r, the function checks if
# point q lies on line segment 'pr'
@inbounds function onSegment(p::Vector{Float64},q::Vector{Float64},r::Vector{Float64})
    return  (q[1] <= max(p[1], r[1]) && q[1] >= min(p[1], r[1]) && q[2] <= max(p[2], r[2]) && q[2] >= min(p[2], r[2]))
end

# To find orientation of ordered triplet (p, q, r).
# The function returns following values
# 0 --> p, q and r are colinear
# 1 --> Clockwise
# 2 --> Counterclockwise
@inbounds function orientation(p::Vector{Float64},q::Vector{Float64},r::Vector{Float64})
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
function doIntersect(p1::Vector{Float64},q1::Vector{Float64},p2::Vector{Float64},q2::Vector{Float64})
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
#TODO optimize this
@inbounds function calc_flux(traj,tspan,config,filename)
    flux = Dict{Any,Any}([ascii(k)=>[] for k in keys(config)])
    tmp = zeros(Float64,4)
    for (k,v) in config
        v = [promote(v...)...]
        k = ascii(k)
        for i = 1:size(traj,3),t_idx = 1:length(tspan)-1#loop over trajectories and time
            (b,d) = doIntersect(v[1:2],v[3:4],traj[1:2,t_idx,i],traj[1:2,t_idx+1,i])
            if b
                tmp[1] = tspan[t_idx]
                tmp[2:3] = traj[1:2,t_idx,i]
                tmp[4] = d
                if isempty(flux[k])
                    flux[k] = tmp
                else
                    flux[k] = cat(2,flux[k],tmp)
                end
            end
        end
    end
    matwrite(filename,flux)
    return flux
end
