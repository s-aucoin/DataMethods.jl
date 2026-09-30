@testset "DataMethods root API" begin
    @testset "IntervalAverage on vectors" begin
        y = [1.0, 3.0, 5.0, 7.0]
        result = IntervalAverage(y, 2.0, 1.0, 1.0)

        @test result.x ≈ [1.0, 3.0]
        @test vec(result.data) ≈ [2.0, 6.0]
    end

    @testset "IntervalAverage on matrices" begin
        A = [1.0 2.0 3.0 4.0;
             2.0 4.0 6.0 8.0]
        result = IntervalAverage(A, 2.0, 1.0, 1.0)

        @test result.x ≈ [1.0, 3.0]
        @test result.data[:, 1] ≈ [1.5, 3.0]
        @test result.data[:, 2] ≈ [3.5, 7.0]
    end

    @testset "IntervalAverage with partial-duty chunks" begin
        y = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
        result = IntervalAverage(y, 2.0, 1.0, 1.0, 0.5, 0.5)

        @test result.x ≈ [1.0, 3.0, 5.0, 7.0]
        @test vec(result.data) ≈ [1.0, 2.5, 4.5, 6.0]

        M = [1.0 2.0 3.0 4.0 5.0 6.0;
             2.0 4.0 6.0 8.0 10.0 12.0]
        matrix_result = IntervalAverage(M, 2.0, 1.0, 1.0, 0.5, 0.5)
        @test size(matrix_result.data) == (2, 4)
        @test matrix_result.x ≈ [1.0, 3.0, 5.0, 7.0]
    end

    @testset "interp_2d fills NaN holes" begin
        A = [1.0 2.0 3.0;
             4.0 NaN 6.0;
             7.0 8.0 9.0]
        B = interp_2d(A)

        @test all(isfinite, B)
        @test B[2, 2] ≈ 5.0 atol = 1e-2
    end
end
