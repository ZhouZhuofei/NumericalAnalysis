"""
    SEq1

Solutions of Equations in One Variable
here are functions to find root:
1. Bisection
2. Fixed-Point Iteration
3. Newton Method
4. The Secant Method
"""
module SEq1
using ForwardDiff

"""
    Bisection(f::Function, a::Real, b::Real; TOL::Float64=0.0000001, N::Int=40, output_seq::Bool=false)

Use Bisection To find a solution to.

``f(x) = 0``

You should use the function, and a interval ``[a, b]``, tolerance TOL, maximum number of iterations N.

Then you will get approximate solution ``p`` or can't find in maximum N and return -1.

you can set output_seq=true to output sequence.
"""
@inline function Bisection(f::Function, a::Real, b::Real; TOL::Float64=0.0000001, N::Int=40, output_seq::Bool=false)
    i = 1
    fa = f(a)
    seq = []
    if f(a)*f(b) > 0
        return -1
    else
        while i <= N
            p = a + (b - a)/2
            fp = f(p)
            push!(seq, [p, f(p)])
            if fp == 0 || (b - a)/2 < TOL
                if output_seq
                    return seq
                else
                    return p
                end
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

"""
    fixed_point(f::Function, p₀::Real; TOL::Float64=0.0000001, N::Int=30, output_seq::Bool=false)

Use fixed point Iteration to find root.

The number p is a **fixed point** for a given function g if ``g(p) = p``.

input a function, initial p₀, optional Args: tolerancs TOL, maximum number of iterations N, output_seq.

Output a root p or sequence.
"""
@inline function fixed_point(f::Function, p₀::Real; TOL::Float64=0.0000001, N::Int=30, output_seq::Bool=false)
    i = 1
    seq = []
    while i <= N
        p = f(p₀)
        push!(seq, p₀)
        if abs(p - p₀) < TOL
            if output_seq
                return seq
            else
                return p
            end
            break
        end
        i += 1
        p₀ = p
    end
    return -1
end

"""
    Newton(f::Function, p₀::Real; TOL::Float64=0.0000001, N::Int=20, output_seq=false)

Use Newton method to find root.

you should use a function and initial p₀, then you can set Args or use default.

"""
@inline function Newton(f::Function, p₀::Real; TOL::Float64=0.0000001, N::Int=20, output_seq::Bool=false)
    i = 1
    seq = []
    f1(x) = ForwardDiff.derivative(x -> f(x), x)
    while i <= N
        push!(seq, p₀)
        p = p₀ - f(p₀)/f1(p₀)
        if abs(p - p₀) < TOL
            if output_seq
                return seq
            else
                return p
            end
            break
        end
        p₀ = p
        i += 1
    end
    return -1
end

"""
    Secant(f::Function, p₀::Real, )

Use the Secant Method to find root
A function and p₀, p₁(p₀<p₁), then you can set Args or use default.
"""
@inline function Secant(f::Function, p₀::Real, p₁::Real; TOL::Float64=0.0000001, N::Int=30, output_seq::Bool=false)
    i = 2
    seq = []
    a₀ = f(p₀)
    a₁ = f(p₁)
    while i <= N
        p = p₁ - a₁*(p₁ - p₀)/(a₁ - a₀)
        push!(seq, p)
        if abs(p - p₁) < TOL
            if output_seq
                return seq
            else
                return p
            end
            break
        end
        i += 1
        p₀ = p₁
        a₀ = a₁
        p₁ = p
        a₁ = f(p)
    end
end





end
