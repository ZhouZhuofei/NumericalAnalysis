using NumericalAnalysis: SEq1
using Test


@testset "SEq1 test" begin
    @test 3.13<=SEq1.Bisection(sin, π/2, 3π/2)<=3.15
    @test 1.3<=SEq1.fixed_point(x -> x - (x^3+4x^2-10)/(3x^2+8x), 1.5)<=1.37
    @test 0.73<=SEq1.Newton(x->cos(x)-x, 0.74)<=0.74
    @test 0.73<=SEq1.Secant(x -> cos(x) - x, 0.5, π/4)<=0.74
    @test 0.73<=SEq1.FalsePos(x->cos(x)-x, 0.5, π/4)<=0.74
end
