function benchmark_smp(f)
    for i = 1:1000000
        Fields.sample2!(f,[50000.0+rand()*10.0,25000.0+rand()*10.0],1.0*rand())
    end
end

function benchmark_value()
    for i = 1:1000000
        Fields.value([50000.0+rand()*10.0,25000.0+rand()*10.0],1.0*rand())
    end
end

@everywhere function itp_test()
    N = 1000
    xx = linspace(48000-1000,52000+1000,N)
    yy = linspace(23000,27000,N)
    output = zeros(Float64,(N,N))
    for x in enumerate(xx), y in enumerate(yy)
        output[x[1],y[1]] = Fields.value([x[2],y[2]],0.75)
    end
    return output
end

@everywhere function itp_grad_test()
    N = 1000
    xx = linspace(48000-1000,52000+1000,N)
    yy = linspace(23000,27000,N)
    output = zeros(Float64,(2,N,N))
    for x in enumerate(xx)
        for y in enumerate(yy)
            output[:,x[1],y[1]] = Fields.gradient2([x[2],y[2]],0.5)
        end
    end
    return output
end

@everywhere function itp_grad_test2()
    N = 1000
    xx = linspace(48000-1000,52000+1000,N)
    yy = linspace(23000,27000,N)
    output = zeros(Float64,(4,N,N))
    grad = Array(Float64,4)
    for x in enumerate(xx)
        for y in enumerate(yy)
            Fields.gradient!(0.5,[x[2],y[2],rand(),rand()],grad)
            output[:,x[1],y[1]] = grad
        end
    end
    return output
end
