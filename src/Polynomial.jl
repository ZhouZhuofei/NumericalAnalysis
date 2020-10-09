"""
    Polynomial

one of functions mapping the set of real numbers is the algebraic polynomials, the set of functions of the form:

``p_n(x) = a_n x^n + a_{n-1} x^{n-1} + ... +a_{1}x + a_0``
"""
module Polynomial
using SymPy

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

end
