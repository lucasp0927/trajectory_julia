struct Polygon
    n::Int64
    x::Vector{Float64}
    y::Vector{Float64}
    constant::Vector{Float64}
    multiple::Vector{Float64}
    function Polygon(poly::Array{Float64})
        if size(poly,1)>2
            poly = reshape(poly,(round(Int64,length(poly)/2)),2)'
        end
        @assert size(poly,1) == 2
        n = size(poly,2)
        x = poly[1,:]
        y = poly[2,:]
        constant = Array{Float64}(undef,n)
        multiple = Array{Float64}(undef,n)
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
            oddNodes = xor(oddNodes,(pos[2]*p.multiple[i]+p.constant[i]<pos[1]));
        end
        j = i
    end
    return oddNodes
end

function anyPointInPolygon(p::Polygon,pos::Array{Float64,2})
    @assert size(pos,1) == 2
    result = false
    for i = 1:size(pos,2)
        result = result || pointInPolygon(p,pos[:,i])
    end
    return result
end
