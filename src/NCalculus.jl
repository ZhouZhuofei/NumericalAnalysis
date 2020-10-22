"""
    NCalculus

Numerical Differentiation and Integration Module.
"""
module NCalculus

"""
    ThreePoint(f::Function, x₀::Real, h::Real; method::String="Endpoint")

input a function `f`, x₀ and h.Args is method (default: "Endpoint", if you set method is't "Endpoint", it will
return MidPoint ans.)
"""
@inline function ThreePoint(f::Function, x₀::Real, h::Real; method::String="Endpoint")
    if method == "Endpoint"
        return 1/h*(-3/2*f(x₀) + 2*f(x₀+h)- 1/2*f(x₀+2*h))
    else
        return 1/(2*h) * (f(x₀+h) - f(x₀- h))
    end
end

"""
    FivePoint(f::Function, x₀::Real, h::Real; method::String="Endpoint")

input a function `f`, x₀ and h.Args is method (default: "Endpoint", if you set method is't "Endpoint", it will
return MidPoint ans.)
"""
@inline function FivePoint(f::Function, x₀::Real, h::Real; method::String="Endpoint")
    if method == "Endpoint"
        return 1/(12*h)*(-25*f(x₀)+48f(x₀ +h)-36f(x₀+2h)+16f(x₀+3h)-3f(x₀+4h))
    else
        return 1/(12*h) * (f(x₀ - 2h) - 8f(x₀ -h)+8f(x₀+h)-f(x₀ + 2h))
    end
end

"""
    Trapezoidal(f::function, a::Real, b::Real)

input a function `f`, and ``x in [a, b]`` calculate ∫f(x)dx.
"""
@inline function Trapezoidal(f::Function, a::Real, b::Real)
    return (b-a)/2 * (f(a) + f(b))
end

"""
    Simpson(f::Function, a::Real, b::Real)

input a function `f`, and ``x in [a, b]`` calculate ∫f(x)dx.
"""
@inline function Simpson(f::Function, a::Real, b::Real)
    return (b-a)/6 * (f(a) + 4*f(a+(b-a)/2) + f(b))
end


"""
    Newton_Cotes(f::Function, a::Real, b::Real; n::Int=4, method::String="Midpoint")

input a function `f`, and ``x in [a, b]``` calculate ``int_a^b f(x)dx``. args ``n∈[1,4]``
"""
@inline function Newton_cotes(f::Function, a::Real, b::Real; n::Int=0, method::String="Midpoint")
    if method == "Open"
        o1 = (b-a)/2 * (f(a) + f(b))
        o2 = (b-a)/6 * (f(a) + 4*f(a+(b-a)/2) + f(b))
        o3 = 3*(b-a)/24 * (f(a) + 3*f(a+(b-a)/3) + 3*f(a+2*(b-a)/3) + f(b))
        o4 = 2*(b-a)/180 * (7*f(a) + 32*f(a+(b-a)/4) + 12*f(a+(b-a)/2) + 32*f(a + 3*(b-a)/4) + 7*f(b))
        o = [o1, o2, o3, o4]
        if 0<n≤4
            return o[n]
        else
            return o
        end
    else
        c1 = 2 * (b-a)/2 * f((b+a)/2)
        c2 = 3/2 * (b-a)/3 * (f((b-a)/3) + f(2*(b-a)/3))
        c3 = 4/3 * (b-a)/4 * (2*f((b-a)/4) - f((b-a)/2) + 2*f(3*(b-a)/4))
        c4 = 5/24 * (b-a)/5 * (11*f((b-a)/5)+f(2*(b-a)/5)+f(3*(b-a)/5)+11f(4*(b-a)/5))
        c = [c1, c2, c3, c4]
        if 0<n≤4
            return c[n]
        else
            return c
        end
    end
end

"""
    Romberg(f::Function, a::Real, b::Real, n::Int; seq_tab::Bool=false)
input a function `f`, and ``x in [a, b]``, set `n`.if Args `seq_tab=true`, will return a table, else return a value.
"""
@inline function Romberg(f::Function, a::Real, b::Real, n::Int; seq_tab::Bool=false)
    h = (b-a)
    result = zeros(n,n)
    result[1,1] = (f(a)+f(b))*h/2
    for i = 2:n
        w = f(a)+f(b)
        for j = 1:(2^(i-1)-1)
            w = w + 2*f(a+j*(h/2^(i-1)))
        end
        result[i,1] = w*(h/2^i)
    end
    for i = 2:n
        for j = i:n
            result[j,i] = result[j,i-1]+(1/(4^(i-1)-1))*(result[j,i-1]-result[j-1,i-1])
        end
    end
    if seq_tab
        return result
    else
        return result[n,n]
    end
end

"""
    Gaussian_Quad(f::Function; n::Int=3, a::Real=-1, b::Real=1)

Gaussian Quadrature must set `n<6`,default calculate `int f(x)dx` in `[-1,1]`. you can set interval `[a, b]` to change default interval.
"""
@inline function Gaussian_Quad(f::Function;n::Int=3,a::Real=-1, b::Real=1)
    #Gaussian Args
    Root = [0.5773502692,-0.5773502692,
        0.7745966692,0.0,-0.7745966692,
        0.8611363116,0.3399810436,-0.3399810436,-0.8611363116,
        0.9061798459,0.5384693101,0.0,-0.5384693101,-0.9061798459]
    Coef = [1.0,1.0,
        0.5555555556,0.8888888889,0.5555555556,
        0.3478548451,0.6521451549,0.6521451549,0.3478548451,
        0.2369268850,0.4786286705,0.5688888889,0.4786286705,0.2369268850]

    if a == -1 && b == 1
        return Coef[sum(1:n-1):sum(1:n-1)+n-1]' * f.(Root[sum(1:n-1):sum(1:n-1)+n-1])

    else
        return Coef[sum(1:n-1):sum(1:n-1)+n-1]' * (b-a)/2 * f.((b-a)/2 .* Root[sum(1:n-1):sum(1:n-1)+n-1] .+ (b+a)/2)
    end
end



end
