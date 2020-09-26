"""
    SEq1

Solutions of Equations in One Variable
"""
module SEq1



"""
    Bisection(f::Function, a::Real, b::Real, TOL::Float64=10^(-4), N::Int=10)

Use Bisection To find a solution to.
``f(x) = 0``
You should use the function, and a interval [a, b], tolerance TOL, maximum number of iterations N₀.
Then you will get approximate solution ``p`` or can't find in maximum N and return -1.
"""
@inline function Bisection(f::Function, a::Real, b::Real, TOL::Float64=0.0001, N::Int=20)
    i = 1
    fa = f(a)
    if f(a)*f(b) > 0
        return -1
    else
        while i <= N
            p = a + (b - a)/2
            fp = f(p)
            if fp == 0 || (b - a)/2 < TOL
                return p
                break
            end
            i += 1
            if fa*fp > 0
                a = p
                fa = fp
            else
                b = p
            end
        end
        return -1
    end
end




end
