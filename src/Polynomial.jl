"""
    Polynomial

one of functions mapping the set of real numbers is the algebraic polynomials, the set of functions of the form:

``p_n(x) = a_n x^n + a_{n-1} x^{n-1} + ... +a_{1}x + a_0``
"""
module Polynomial
using SymPy
using LinearAlgebra
using Latexify

@inline function product(x::Vector)
    res = 1
    for i in 1:size(x)[1]
        res = res*x[i]
    end
    return res
end

"""
    Lagrange(X::Vector, Y::Vector; a::String="x")

use the ``(x_0, y_0), (x_1, y_1), ..., (x_n, y_n)``, to caculate a polynomial.
input X and Y, default use x to express functions.
output: function.
"""
@inline function Lagrange(X::Vector, Y::Vector; a::String="x")
    a = symbols(a)
    m = size(X)[1]
    A = ones(m) * X'
    An = A - A'
    An[An .== 0] .= 1
    up = a .- X
    upr = product(up)
    sum = 0
    for i in 1:m
        down = product(An[:,i])
        u = upr/up[i]
        sum = sum + Y[i] * (u/down)
    end
    return simplify(sum)
end


"""
    Neville(x::Vector, y::Vector, x₀::Real; tab::Bool=false)

input two vector, x and y. satisfy the function `y=f(x)`, then input x₀, which you want
to  approximate f(x₀). The default tab setting is false, just output the estimate. you can set true
to output a result table.
"""
@inline function Neville(x::Vector, y::Vector, x₀::Real; tab::Bool=false)
    n, = size(x)
    res = zeros(n, n)
    for i in 1:n
        res[i, 1] = y[i]
    end
    for i in 2:n
        for j in 2:i
            res[i,j] = ((x₀ - x[i-j+1])*res[i, j-1] - (x₀ - x[i])*res[i-1,j-1])/(x[i] - x[i-j+1])
        end
    end
    if tab
        return res
    else
        return res[n, n]
    end
end

"""
    NDDF(x::Vector, y::Vector; a::String="x", simple::Bool=true, tab::Bool=false, backward::Bool=false)

input ``x_1, x_2, x_3, ...x_n```,values ``f(x_1), f(x_2), ...,f(x_n)``, set Args `simpel = true`
output simplify ans, set tab to output a table. set `backward=true` use the Newton Backward–Difference Formula.
"""
@inline function NDDF(x::Vector, y::Vector; a::String="x", simple::Bool=true, tab=false, backward::Bool=false)
    a = symbols(a)
    n, = size(x)
    res = zeros(n, n)
    for i in 1:n
        res[i,1] = y[i]
    end
    for i in 2:n
        for j in 2:i
            res[i,j] = (res[i, j-1] - res[i-1, j-1])/(x[i] - x[i-j+1])
        end
    end
    sum = res[1,1]
    for i in 2:n
        p = 1
        for j in 1:i-1
            p = p*(a - x[j])
        end
        if backward
            sum = sum + res[n, i]*p
        else
            sum = sum + res[i, i]*p
        end
    end
    if tab
        return res
    else
        if simple
            return simplify(sum)
        else
            return sum
        end
    end
end

"""
    NCSpline( NCSpline(x::Vector, y::Vector; Latex::Bool=false, a::String="t")
input ``x_1, x_2, x_3, ...x_n```,values ``f(x_1), f(x_2), ...,f(x_n)``, set Args `Latex = true`  to
output information about the function.
"""
@inline function NCSpline(x::Vector, y::Vector; Latex::Bool=false, a::String="t")
    a = symbols(a)
    n, = size(x)
    h, α = zeros(n), zeros(n-1)

    for i in 1:n-1
        h[i] = x[i+1] - x[i]
    end
    for i in 2:n-1
        α[i] = 3/h[i] * (y[i+1] - y[i]) - 3/h[i-1] *(y[i] - y[i-1])
    end

    l = ones(n)
    μ = zeros(n)
    z = zeros(n)

    for i in 2:n-1
        l[i] = 2*(x[i+1] - x[i-1]) - h[i-1]*μ[i-1]
        μ[i] = h[i]/l[i]
        z[i] = (α[i] - h[i-1]*z[i-1])/l[i]
    end

    c, b, d = zeros(n), zeros(n), zeros(n)
    for j in reverse(1:n-1)
        c[j] = z[j] - μ[j]*c[j+1]
        b[j] = (y[j+1] - y[j])/h[j] - h[j]*(c[j+1]+2*c[j])/3
        d[j] = (c[j+1] - c[j])/(3*h[j])
    end
    A = [b[1:end-1] c[1:end-1] d[1:end-1]]
    X = []
    for i in 2:4
        for j in 1:3
            push!(X,(a - x[i-1])^(j))
        end
    end
    res = diag(A * reshape(X, (n-1,n-1))) + y[1:end-1]

    if Latex
        for i in 1:n-1
            print(res[i], ",","t ∈","[$(x[i]), $(x[i+1])]")
            print("\n")

        end
    else
        return res
    end
end





end
