using DimensionalData
using Rasters

@testset "CTDUtilities" begin
    @testset "lag correction" begin
        T = collect(1.0:10.0)
        fs = 2
        result = correct_CT_lag(T, fs; τ=1.0)
        @test result[1:end-1] == [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
        @test isnan(result[end])
    end

    @testset "remove_loops" begin
        test_p = [0.1, 1.0, 2.0, 1.5, 0.5, 3.0]
        @test FindLoops(test_p) == [4, 5]

        z = DimArray(-test_p, (Dim{:z}(1:6),))
        x = DimArray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], (Dim{:z}(1:6),))
        ds = DimStack((z=z, x=x))
        cleaned_ds = remove_loops(ds; tol=0.0)
        @test Array(cleaned_ds.z) == [-0.1, -1.0, -2.0, -3.0]
        @test Array(cleaned_ds.x) == [1.0, 2.0, 3.0, 6.0]

        z = Raster(-test_p, dims=(Dim{:z}(1:6),))
        x = Raster([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dims=(Dim{:z}(1:6),))
        rs = RasterStack((z=z, x=x))
        cleaned_rs = remove_loops(rs; tol=0.0)
        @test Array(cleaned_rs.z) == [-0.1, -1.0, -2.0, -3.0]
        @test Array(cleaned_rs.x) == [1.0, 2.0, 3.0, 6.0]
    end

    @testset "find_ctd_profiles returns cast ranges" begin
        pressure = vcat(0.0:0.1:10.1, 10.0:-0.1:0.1, 0.2, 0.1, 0.2, 0.1, 0.1)
        casts = find_ctd_profiles(pressure; sm=3, Pthresh=1.0, stdthresh=0.5, nomotion_tol=0.05)
        @test length(casts.downcast_inds) >= 1
        @test all(length.(casts.downcast_inds) .> 1)
        @test all(length.(casts.upcast_inds) .> 1)
    end
end
