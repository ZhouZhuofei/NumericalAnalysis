using NumericalAnalysis
using Test


@testset "NumericalAnalysis.jl" begin
    #2x + 4y
    @test my_f(1, 2) == 10
    @test my_f(1.0, 3) == 14.0
    # Write your tests here.
end
