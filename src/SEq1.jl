"""
    SEq1

Solutions of Equations in One Variable
here are functions to find root:
1. Bisection
2. Fixed-Point Iteration
3. Newton Method
4. The Secant Method
5. The False Position Method
6.  Müller’s Method
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
    push!(seq, p₀)
    push!(seq, p₁)
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

"""
    FalsePos(f::Function, p₀::Real, p₁::Real; TOL::Float64=0.00000001, N::Int=40, output_seq::Bool=false)

Use the false position method to find root.
a function and interval ``[a, b]``, then you can set Args or use default.
"""
@inline function FalsePos(f::Function, p₀::Real, p₁::Real; TOL::Float64=0.00000001, N::Int=40, output_seq::Bool=false)
    i = 2
    seq = []
    q₀ = f(p₀)
    q₁ = f(p₁)
    push!(seq, p₀)
    push!(seq, p₁)
    while i <= N
        p = p₁ - q₁*(p₁ - p₀)/(q₁ - q₀)
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
        q = f(p)
        if q*q₁ < 0
            p₀ = p₁
            q₀ = q₁
        end
        p₁ = p
        q₁ = q
    end
    return -1
end

"""
    ModifiedNewton(f::Function, p₀::Real; TOL::Float64=0.00000001, N::Int=20, output_seq::Bool=false)

Newton and ModifiedNewton, Both methods are rapidly convergent to the actual zero.
"""
@inline function ModifiedNewton(f::Function, p₀::Real; TOL::Float64=0.00000001, N::Int=20, output_seq::Bool=false)
    i = 1
    seq = []
    push!(seq, p₀)
    f1(x) = ForwardDiff.derivative(x -> f(x), x)
    f11(x) = ForwardDiff.derivative(x -> f1(x), x)
    while i <= N
        p = p₀ - f(p₀)*f1(p₀)/((f1(p₀))^2 - f(p₀)*f11(p₀))
        push!(seq, p)
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
    Muller(f::Function, p₀::Real, p₁::Real, p₂::Real; TOL::Float64=0.00000001, N::Int=20, output_seq::Bool=false)
"""
@inline function Muller(f::Function, p₀::Real, p₁::Real, p₂::Real; TOL::Float64=0.00000001, N::Int=20, output_seq::Bool=false)
    h₁ = p₁ - p₀
    h₂ = p₂ - p₁
    δ₁ = (f(p₁) - f(p₀))/h₁
    δ₂ = (f(p₂) - f(p₁))/h₂
    d = (δ₂ - δ₁)/(h₁ + h₂)
    i = 3
    seq = []
    while i <= N
        b = δ₂ + h₂*d
        D = sqrt(complex(b^2 - 4*f(p₂)*d))
        if abs(b - D) < abs(b + D)
            E = b + D
        else
            E = b - D
        end
        h = -2 * f(p₂)/E
        p = p₂ + h
        push!(seq, p)
        if abs(h) < TOL
            if output_seq
                return seq
            else
                return p
            end
            break
        end
        p₀ = p₁
        p₁ = p₂
        p₂ = p
        h₁ = p₁ - p₀
        h₂ = p₂ - p₁
        δ₁ = (f(p₁) - f(p₀))/h₁
        δ₂ = (f(p₂) - f(p₁))/h₂
        d = (δ₂ - δ₁)/(h₁ + h₂)
        i += 1
    end
    return -1
end

end
