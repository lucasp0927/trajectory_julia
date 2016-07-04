using HDF5
using MATLAB

function fit_trap_matlab(data,axial_t,radial_t)
    prob,xstart,xend,ystart,yend = mxcall(:fit_trap,5,data,axial_t,radial_t)
    println(size(prob))
    println(xstart)
    println(xend)
    println(ystart)
    println(yend)
end

#coeff = [fitr.p00, fitr.p10, fitr.p01, fitr.p20, fitr.p11, fitr.p02, fitr.p30, fitr.p21, fitr.p12, fitr.p03, fitr.p40, fitr.p31, fitr.p22, fitr.p13, fitr.p04];
