struct FakeAttrs
    attrib::Dict
end

struct FakeDataGroup
    group::Dict
end

struct FakeDataset
    group::Dict
end

@testset "SignalAnalysis" begin
    @testset "simple SNR formulas" begin
        @test calculate_SNR(20.0, 10.0; mode=:log) == 10.0
        @test calculate_SNR(100.0, 10.0; mode=:lin) ≈ 20.0
    end

    @testset "pulse and range calculations" begin
        @test PulseLag(2.0, 1500.0) ≈ 2 * 2 / 1500
        @test IntrinsicVelocityRange(1500.0, 300000.0, 0.5) ≈ 1500 / (4 * 300000 * 0.5)
        @test ExtendedVelocityRange(1500.0, 300000.0, 0.5; nmax=2) ≈ 2 * 2 * IntrinsicVelocityRange(1500.0, 300000.0, 0.5)
    end

    @testset "correlation and turbulence relationships" begin
        d = 0.5
        Vr = 1.5
        ε = 1e-4
        cor = CorrelationFromTurbulence(ε, d, Vr)
        @test cor > 0.0
        @test cor < 1.0
        @test TurbulenceFromCorrelation(cor, d, Vr) ≈ ε rtol=1e-6
    end

    @testset "variance and correlation overloads" begin
        @test HR_V_variance(0.7, 1.0, 100_000.0, 10_000.0, 0.1; cs=1500.0) > 0.0
        @test HR_c_from_V_variance(1.0, 1.0, 100_000.0, 10_000.0, 0.1; cs=1500.0) > 0.0
        @test HR_c_from_V_variance(1.0, 1.0, 100_000.0, 10_000.0, 0.1; cs=1500.0) <= 1.0

        ds = FakeDataset(Dict(
            "Config" => FakeAttrs(Dict(
                "bursthr_cellSize" => 1.0,
                "beamConfiguration5_frequency" => 100.0,
                "beamConfiguration5_bandwidth" => 10.0,
                "bursthr_lag" => 0.1,
            )),
            "Data" => FakeDataGroup(Dict(
                "IBurstHR" => Dict("SpeedOfSound" => [1500.0]),
            )),
        ))

        var = HR_V_variance([0.7], ds)
        corr = HR_c_from_V_variance([1.0], ds)
        @test all(isfinite, var)
        @test all(isfinite, corr)
    end
end
