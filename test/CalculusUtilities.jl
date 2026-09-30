@testset "CalculusUtilities" begin
    x = 0.0:0.5:2.0
    f = x .^ 2

    @testset "Centered and backward derivatives" begin
        @test ∂_cen(f, x)[2:end-1] ≈ [1.0, 2.0, 3.0]
        @test ∂_back(f, x)[2:end] ≈ [0.5, 1.5, 2.5, 3.5]
    end

    @testset "Trapezoid integrals" begin
        @test int_def_trap(f, 0.5) ≈ 2.75
        @test int_trap(f, x) ≈ [0.0, 0.0625, 0.375, 1.1875, 2.75]
    end
end
