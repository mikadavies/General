# Visualisation of the Hydrogen orbitals

# Includes
using ForwardDiff
using Plots

# Constants
aₒ = 5.2917721090380*10^(-11) # m

# Mathematical functions
D(f) = x -> ForwardDiff.derivative(f, float(x))
D(f, n) = n > 1 ? D(D(f), n-1) : D(f)

function Laguerre(n, a, x)
    if n == 0
        return 1
    elseif n == 1
        return 1+a-x
    else
        k = n-1
        L = ((2*k+1+a-x)*Laguerre(k, a, x)-(k+a)*Laguerre(k-1, a, x))/(k+1)
        return L
    end
end

function Binomial(n, k)
    if n == 1/2
        return Binomial(2*k, k)*(((-1)^(k+1))/((2^(2*k))*(2*k-1)))
    else
        i = 1
        p = 1
        while i <= k
            p *= (n+1-i)/i
            i += 1
        end
        return p
    end
end

function Legendre(l, m, x)
    k = m
    Summation = 0
    while k <= l
        Summation += (factorial(big(k))/factorial(big(k-m)))*(x^(k-m))*Binomial(l, k)*Binomial((l+k-1)/2, l)
        k += 1
    end
    return (2^l)*((1-x^2)^(m/2))*Summation
end

function polar2cartesian(r, θ)
    x = r*cos(θ)
    y = r*sin(θ)
    return x, y
end

function cartesian2polar(x, y)
    r, θ = sqrt(x^2 + y^2), atan(x, y) # mathematically incorrect, but adjusted to switched axes in our functions
    return r, θ
end

# Quantum mechanical functions

function Harmonic(l, m, θ, ϕ)
    return ((-1)^m)*sqrt(((2*l+1)/(4*π))*(factorial(big(l-m))/factorial(big(l+m))))*Legendre(l, m, cos(θ))*exp(im*m*ϕ)
end

function Hydrogen_ψ(n, l, m, r, θ, ϕ)
    h = 1.054571817*10^(-34)
    ve = 299792458/137
    me = 9.109383701528*10^(-31)
    εo = 8.854187812913*10^(-12)
    e = 1.602176634*10^(-19)
    α = (1/(4*π*εo))*((e^2)/(h*299792458))
    a = (h*sqrt(1-(ve/299792458)^2))/(me*299792458*α)

    ρ = (2*r)/(n*a)
    L = Laguerre(n-l-1, 2*l+1, ρ)
    Y = Harmonic(l,m,θ,ϕ)

    ψ = sqrt(Complex(((2/(n*a))^3)*(factorial(big(n-l-1))/(2*n*factorial(big(n+l))))))*exp(-ρ/2)*(ρ^l)*L*Y

    return ψ

end

Hydrogen_P(n, l, m, r, θ, ϕ) = abs2(Hydrogen_ψ(n,l,m,r,θ,ϕ))

# Utilities
function nlm(max)
    all = []
    for i in 1:max
        n = i
        for j in 0:(n-1)
            l = j
            for k in 0:l
                m = k
                append!(all, [(n,l,m)])
            end
        end
    end
    return all
end

# Visualisation
function plotP(n, l, m, rmax, path, steps=100, c_lim=NaN, legend=false) #rmax in terms of aₒ
    function f(x, y)
        r, θ = cartesian2polar(x*aₒ, y*aₒ)
        P = Hydrogen_P(n, l, m, r, θ, 0)
        isnan(P) ? P=0 : P=P
        return P
    end
    rs = range(-rmax, rmax, length=steps)
    if isnan(c_lim)
        legend == false ? heatmap(rs, rs, f, c=:default, axis=nothing, size=(1000,1000), legend=:none) : heatmap(rs, rs, f, c=:default, axis=nothing, size=(1000,1000))
    else
        legend == false ? heatmap(rs, rs, f, c=:default, axis=nothing, clim=(0,c_lim), size=(1000,1000), legend=:none) : heatmap(rs, rs, f, c=:default, axis=nothing, clim=(0,c_lim), size=(1000,1000))
    end
    annotate!(0.9rmax, -0.9rmax, text("($n,$l,$m)", :white, "TimesNewRoman", 15))
    fname = string(path, "\\$n$l$m")
    png(fname)
    println("$n,$l,$m")
end

function showAll(max, path, steps=100)
    indices = nlm(max)
    for i in indices
        if i[1] > 4
            rmax = 10*i[1]+i[1]*i[2]+i[3]
        elseif i[1] < 3
            rmax = 10
        else
            rmax = 10*i[1]+i[2]+0.1*i[3]
        end
        if sum(i) > 8
            c_lim = sum(i)/3*(0.3^sum(i))*0.2e31
        elseif sum(i) <= 5
            c_lim = sum(i)/3*((0.15^(i[1]+i[2]))*0.2e31)
        else
            c_lim = sum(i)/3*((0.2^(i[1]+i[2]))*0.2e31)
        end
        fullP(i[1],i[2],i[3],rmax,steps,c_lim,false,path)
    end
    println("FINISHED!")
end
