@testset "WaveAnalysis" begin
    @testset "projection utilities" begin
        @test DataMethods.unit_spherical_vector_projection_u(0.0, 0.0, 0.0) == 0.0
        @test DataMethods.unit_spherical_vector_projection_w(0.0, 0.0) == 1.0
        @test DataMethods.rfactor2x(0.0, 0.0, 0.0, 0.0) == 0.0
        @test DataMethods.rfactor2z(0.0, 0.0, 0.0) == -1.0
    end

    @testset "wave field helpers" begin
        r = [1.0, 2.0]
        h = 5.0
        θ = 20.0
        θ₀ = 20.0
        β = 0.0
        α = 0.0
        A = [1.0, 2.0]
        ω = [0.5, 0.7]
        ϕ = [0.0, 0.2]
        k = [0.1, 0.12]

        v = DataMethods.wave_from_velocity(A, ω, ϕ, r, h, k, θ, θ₀, β, α)
        a = DataMethods.wave_from_acceleration(A, ω, ϕ, r, h, k, θ, θ₀, β, α)

        @test size(v) == (2, 2)
        @test size(a) == (2, 2)
        @test all(isfinite, v)
        @test all(isfinite, a)
    end

    @testset "wave_U_ofz integrates IMF contributions" begin
        IA = [1.0 0.5 0.2;
              2.0 1.0 0.4;
              3.0 1.5 0.6]
        IF = [0.10 0.15 0.20;
              0.12 0.18 0.22;
              0.14 0.21 0.24]
        IP = [0.0 0.5 0.9;
              0.2 0.9 1.4;
              0.4 1.3 1.8]
        r = [1.0, 2.0, 3.0]
        h = 5.0

        out = wave_U_ofz(IA, IF, IP, r, h, 20.0, 20.0, 0.0, 0.0; method=:velocity, incl_trend=false, IMFs2incl=1:2)
        @test size(out) == (3, 3)
        @test all(isfinite, out)
        @test maximum(abs, out) > 0.0
    end
end
