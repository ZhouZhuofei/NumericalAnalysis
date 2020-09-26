"""
    SEq1

Solutions of Equations in One Variable
here are functions to find root:
1. Bisection
2. ...
"""
module SEq1



"""
    Bisection(f::Function, a::Real, b::Real; TOL::Float64=10^(-4), N::Int=10, output_seq::Bool=false)

Use Bisection To find a solution to.
``f(x) = 0``
You should use the function, and a interval [a, b], tolerance TOL, maximum number of iterations N₀.
Then you will get approximate solution ``p`` or can't find in maximum N and return -1.
you can set output_seq=true to output seqencen.
"""
@inline function Bisection(f::Function, a::Real, b::Real; TOL::Float64=0.0001, N::Int=20, output_seq::Bool=false)
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





end
