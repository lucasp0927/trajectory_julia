@fastmath function wg_transfer_matrix( k::Float64 , l::Float64 )
    M = [exp(1im*k*l) zero(Complex{Float64});zero(Complex{Float64}) exp(-1im*k*l)]
    @assert abs(1-abs(det(M)))<1e-10 "det(M)!=1"
    return M
end

function atom_transfer_matrix(f, f_0, gamma_1d, gamma_prime)
    gamma = gamma_1d::Float64+gamma_prime::Float64
    delta::Float64 = f-f_0
    r::Complex{Float64} = -gamma_1d/(gamma-2.0*im*delta)
    t::Complex{Float64} = one(Complex{Float64})+r
    #t_rec = one(Complex{Float64}) #make it slower

    M11::Complex{Float64} = t-(r^2)/t;
    M12::Complex{Float64} = r/t
    M21::Complex{Float64} = -M12#-r/t
    M22::Complex{Float64} = one(Complex{Float64})/t
    #=
    M22::Complex{Float64} = one(Complex{Float64})/t
    M11::Complex{Float64} = t-(r^2)/t;
    M12::Complex{Float64} = r/t
    M21::Complex{Float64} = -M12#-r/t
    =#
    M = [M11 M12;M21 M22]
    @assert abs2(t)<=1.0 "abs2(t)>1.0! f: $f, f_0: $f_0, gamma_1d: $gamma_1d, gamma_prime: $gamma_prime"
    @assert abs(1-abs(det(M)))<1e-10 "det(M)!=1"
    return M
end

@fastmath function atom_transfer_matrix_deutsch(f, f_0, p, gamma)
    delta_norm::Float64 = (f-f_0)/gamma
    zeta::Complex{Float64} = p*(-2.0*delta_norm+1.0im)/(1.0+4.0*delta_norm^2)
    t::Complex{Float64} = one(Complex{Float64})/(1.0-im*zeta)
    r::Complex{Float64} = 1im*zeta/(1.0-im*zeta)
    #t_rec = one(Complex{Float64}) #make it slower
    M22::Complex{Float64} = one(Complex{Float64})/t
    M11::Complex{Float64} = conj(M22)#1/conj(t)#t-r^2/t;
    M12::Complex{Float64} = r/t
    M21::Complex{Float64} = -M12#-r/t
    M = [M11 M12;M21 M22]
    @assert abs(1-abs(det(M)))<1e-10 "det(M)!=1"
    return M
end
