#just use for test
my_f(x, y) = 2x + 4y
derivative_my_f(x, y) = ForwardDiff.derivative(x -> my_f(x, y), x)

"""
    absolute_error(x, x̂)

return the Absolute error of x and x̂
"""
absolute_error(x, x̂) = abs(x .- x̂)




"""
    relative_error(x, x̂)


The relative error between X1 and X2
"""
relative_error(x, x̂) = abs((x .- x̂) ./ x)


"""
    ξ(x₀, x)

return the number ξ(x) bwteen x₀ and x.
"""

ξ(x₀, x) = rand(1)[1] * (x - x₀)/x + x₀


"""
    NthDrivative(f, x, n)

return the nth derivative value at x
"""

function NthDerivative(f, x, n)
    if n == 0
        return f(x)
    else
        return ForwardDiff.derivative(x->NthDerivative(f, x, n-1), x)
    end
end




"""
    TaylorPolynomials(f, x, n)

return the value of nth Taylor Ploynomial and the result.
"""


function TaylorPolynomials(f, x, x₀, n)
    result = []
    for i in 0:n
        push!(result, NthDerivative(f, x₀, i)*(x - x₀)^i/factorial(i))
    end
    return sum(result)

end
