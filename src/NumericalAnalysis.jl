module NumericalAnalysis
using Documenter
using ForwardDiff
using LinearAlgebra
using PyCall
using Conda
using SymPy
include("basic.jl")
export my_f, derivative_my_f, ξ, absolute_error, relative_error, NthDerivative, TaylorPolynomials

# Write your package code here.

include("SEq1.jl")

export Bisection, fixed_point, Newton, Secant, FalsePos, ModifiedNewton

include("Polynomial.jl")

export Lagrange

end
