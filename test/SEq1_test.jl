using NumericalAnalysis: SEq1
using Test


@testset "SEq1 test" begin
    @test 3.13<=SEq1.Bisection(sin, π/2, 3π/2)<=3.15
end
