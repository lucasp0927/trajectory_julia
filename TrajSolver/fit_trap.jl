using HDF5
try
    using MATLAB
end

function fit_trap(data,axial_t,radial_t)
    if axial_t == radial_t
        Lumberjack.info("axial temperature = radial temperature, no need for matlab.")
#        prob,xstart,xend,ystart,yend = mxcall(:fit_trap,5,data,axial_t,radial_t)
        prob,xstart,xend,ystart,yend = iso_temp_prob_dist(data,axial_t)
    else
        prob,xstart,xend,ystart,yend = mxcall(:fit_trap,5,data,axial_t,radial_t)
    end
    return prob,xstart,xend,ystart,yend
end

function iso_temp_prob_dist(data,temp)
    xx = squeeze(data[1,:,:],1)
    yy = squeeze(data[2,:,:],1)
    uu = squeeze(data[3,:,:],1)
    uu_norm = uu-minimum(uu)
    prob = exp(-uu_norm/temp)
    prob = prob/(maximum(prob)*1.01)
    xstart = xx[1,1]
    xend = xx[end,1]
    ystart = yy[1,1]
    yend = yy[1,end]
    return prob,xstart,xend,ystart,yend
end
