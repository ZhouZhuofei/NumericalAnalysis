#Basic function to use

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
@inline ξ(x₀, x) = rand(1)[1] * (x - x₀)/x + x₀


"""
    NthDrivative(f::Function, x::Real, n::Int)

return the nth derivative value at x
"""
@inline function NthDerivative(f::Function, x::Real, n::Int)
    if n == 0
        return f(x)
    else
        return ForwardDiff.derivative(x->NthDerivative(f, x, n-1), x)
    end
end




"""
    TaylorPolynomials(f::Function, x::Real, x₀::Real, n::Int)

return the value of nth Taylor Ploynomial at x.
"""
@inline function TaylorPolynomials(f::Function, x::Real, x₀::Real, n::Int)
    result = []
    for i in 0:n
        push!(result, NthDerivative(f, x₀, i)*(x - x₀)^i/factorial(i))
    end
    return sum(result)

end
