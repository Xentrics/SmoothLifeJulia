#=
#Hack: in the paper, λ is different for n and m! The smoothness allows us to use smaller masks for similar precision
function σ1(x::Float64, a::Float64)
    return 1.0 ./ (1.0 + exp(-4.0.*(x-a)))
end

function σ2(x::Float64, a::Float64, b::Float64)
    return σ1(x,a).*(1.0 - σ1(x,b))
end

function σm(x::Float64, y::Float64, m::Float64)
    return x.*(1.0 - σ1(m, 0.5)) + y.*σ1(m, 0.5)
end

"""
Evaluate inner and outer filling of a cell. The used formulars and constants were proposed for smooth, stable gliders.
These can be found in the original smoothlife paper from 2011.

The function is a mapping of [0,1) x [0,1) ⟶ [0,1).

#Arguments
* `n::Float64`: outer filling ∈ [0,1]
* `m::Float64`: inner filling ∈ [0,1]
"""
function s(n::Float64, m::Float64)
    #Note: birth interval [b1,b2] and death interval [d1, d2]. Values based on original paper
    b1 = 0.278
    b2 = 0.365
    d1 = 0.267
    d2 = 0.445
    return σ2(n, σm(b1,d1,m), σm(b2,d2,m))
end
=#

"""
smoothed integral step into the next generation.
    f(x, t + dt) = f(x) + dt * (2*s(n,m)-1) * f(x)

#Arguments
* `f::Float64`: value of the current cell (x,y) with f = f(x,y), n = mask_outer(x,y) , m = mask_inner(x,y)
* `n::Float64`: outer filling ∈ [0,1)
* `m::Float64`: inner filling ∈ [0,1)
* `dt::Float64`: distance in time
"""
function smoothStep(f::Float64, n::Float64, m::Float64, dt::Float64)
    return f + dt*(2*snm(n,m)-1)
end


function smoothStep(f::Float64, n::Float64, m::Float64)
    return s(n,m)*f
end

# alternative solution
sigmoid_a(x, a, b) =  1 ./(1+exp(-4(x-a)/b))

sigmoid_b(x, b, eb) = 1 - sigmoid_a(x, b, eb)

sigmoid_ab(x, a, b, ea, eb) = sigmoid_a(x, a, ea) .* sigmoid_b(x, b, eb)

function sigmoid_mix(x, y, m, em)
    return x.*(1-sigmoid_a(m, 0.5, em)) + y.*sigmoid_a(m, 0.5, em)
end

function snm(n, m)
    b1 = 0.278
    b2 = 0.365
    d1 = 0.267
    d2 = 0.445
    alphan = 0.028
    alpham = 0.147
    return sigmoid_ab(n, sigmoid_mix(b1, d1, m, alpham), sigmoid_mix(b2, d2, m, alpham), alphan, alphan)
end
