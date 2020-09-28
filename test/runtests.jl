using SafeTestsets
@safetestset "Basic function test" begin include("basic_test.jl") end
@safetestset "SEq for one variable" begin include("SEq1_test.jl") end
@safetestset "Polynomial" begin include("Polynomial_test.jl") end
