using NumericalAnalysis
using Test


@testset "NumericalAnalysis.jl" begin
    #2x + 4y
    # Write your tests here.
end

@testset "basic test" begin
    @test 1<=ξ(1, 2)<=2
    @test absolute_error(1.0, 2) == 1.0
    @test relative_error(1, 2) == 1.0
    @test NthDerivative(exp, 4, 4)  == exp(4)
    @test NthDerivative(x->x^3, 4, 4) == 0
    @test 0.9<=TaylorPolynomials(cos, 0.1, 0, 6)<=1.0




end
