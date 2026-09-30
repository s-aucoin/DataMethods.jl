@testset "CoordinateTransformations" begin
    @testset "simple range conversions" begin
        @test r_2_slant_r(10.0, 30.0) ≈ 10.0 / cosd(30.0)
        @test z_2_brange(5.0, 30.0, 0.0, 0.0, 0.0) ≈ 5.0 * sqrt(1.0 + tand(30.0)^2)
        @test brange_2_z(10.0, 30.0, 0.0, 0.0, 0.0) ≈ 10.0 / sqrt(1.0 + tand(30.0)^2)
    end

    @testset "align_z0 re-centers data" begin
        data = [1.0 2.0;
                3.0 4.0]
        zs = [0.0 1.0;
              1.0 2.0]

        result = align_z0(data, zs, 0.5, 1.0)
        @test length(result.ref_z) >= 2
        @test size(result.data) == (length(result.ref_z), 2)
        @test all(isfinite, result.z_true)
    end
end
