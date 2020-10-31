using NumericalAnalysis: Polynomial
using Test

x = [1.0, 1.3, 1.6, 1.9, 2.2]
y = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]

Polynomial.NCSpline([0,1,2,3], [1, exp(1), exp(2), exp(3)])
Polynomial.CCSpline([0,1,2,3], [1, exp(1), exp(2), exp(3)],[1, exp(3)])
@testset "Polynomial test" begin
    @test abs(Polynomial.Lagrange(x,y)(1.5) - 0.5118199)<=0.0000001
    @test abs(Polynomial.NDDF(x,y)(1.5) - 0.5118199)<=0.0000001
    @test abs(abs(Polynomial.Neville(x,y,1.5) - 0.5118199)<=0.0000001)
end
