type Polygon
    n::Int64
    x::Vector{Float64}
    y::Vector{Float64}
    constant::Vector{Float64}
    multiple::Vector{Float64}
    function Polygon(poly::Array{Float64,2})
        @assert size(poly,1) == 2
        n = size(poly,2)
        x = squeeze(poly[1,:],1)
        y = squeeze(poly[2,:],1)
        constant = Array(Float64,n)
        multiple = Array(Float64,n)
        j = n
        for i = 1:n
            if y[j]==y[i]
                constant[i]=x[i]
                multiple[i]=0.0;
            else
                constant[i]=x[i]-(y[i]*x[j])/(y[j]-y[i])+(y[i]*x[i])/(y[j]-y[i]);
                multiple[i]=(x[j]-x[i])/(y[j]-y[i]);
            end
            j=i
        end
        new(n,x,y,constant,multiple)
    end
end

@inbounds function pointInPolygon(p::Polygon,pos::Vector{Float64})
    n = p.n::Int64
    j = n
    oddNodes = false
    for i = 1:n
        if ((p.y[i] < pos[2] <= p.y[j] || p.y[j] < pos[2] <= p.y[i]))
            oddNodes $= (pos[2]*p.multiple[i]+p.constant[i]<pos[1]);
        end
        j = i
    end
    return oddNodes
end
