@testset "VelocityTransformations" begin
    @testset "beam geometry matrices" begin
        T = transformation_matrix(20.0)
        @test size(T) == (4, 4)
        @test isfinite(T[1, 1])

        H = heading_R_matrix(0.0)
        P = pitchroll_R_matrix(0.0, 0.2)
        @test all(isfinite, H)
        @test all(isfinite, P)
    end

    @testset "beam to xyz and ENU conversions" begin
        xyz = beam2xyz(1.0, 2.0, 3.0, 4.0; θ=20.0)
        @test all(isfinite, (xyz.u, xyz.v, xyz.w1, xyz.w2))

        V1 = [1.0 2.0; 3.0 4.0]
        V2 = [2.0 3.0; 4.0 5.0]
        V3 = [3.0 4.0; 5.0 6.0]
        V4 = [4.0 5.0; 6.0 7.0]
        heading = [90.0, 90.0]
        pitch = [0.0, 0.0]
        roll = [0.0, 0.0]
        R = rotation_matrix(heading, pitch, roll)
        enu = beam2ENU(V1, V2, V3, V4, heading, pitch, roll; θ=20.0)
        AHRSraw = zeros(Float64, 3, 3, 3)
        enu_ahrs = beam2ENU(V1, V2, V3, V4, AHRSraw; T=transformation_matrix(20.0))

        @test size(R) == (4, 4, 2)
        @test all(isfinite, enu.uE)
        @test all(isfinite, enu.vN)
        @test all(isfinite, enu.U1)
        @test all(isfinite, enu.U2)
        @test all(isfinite, enu_ahrs.uE)
        @test all(isfinite, enu_ahrs.vN)
        @test all(isfinite, enu_ahrs.U1)
        @test all(isfinite, enu_ahrs.U2)
    end

    @testset "AHRS rotation matrix" begin
        AHRSraw = zeros(Float64, 3, 3, 3)
        R = rotation_matrix(AHRSraw)
        @test size(R) == (4, 4, 3)
        @test all(isfinite, R)
    end

    @testset "RDI transforms" begin
        @test isfinite(RDI_pitch2pitch(10.0, 5.0))

        R = RDI_rotation_matrix(10.0, 5.0, 30.0)
        @test size(R) == (3, 3)
        @test all(isfinite, R)

        uE = [1.0 2.0; 3.0 4.0]
        vN = [2.0 3.0; 4.0 5.0]
        U = [3.0 4.0; 5.0 6.0]
        xyz = RDI_ENU2XYZ(uE, vN, U, 10.0, 5.0, 30.0)
        enu = RDI_XYZ2ENU([1.0 2.0; 3.0 4.0], [2.0 3.0; 4.0 5.0], [3.0 4.0; 5.0 6.0], 10.0, 5.0, 30.0)

        @test all(isfinite, xyz.u)
        @test all(isfinite, xyz.v)
        @test all(isfinite, xyz.w)
        @test all(isfinite, enu.uE)
        @test all(isfinite, enu.vN)
        @test all(isfinite, enu.U)
    end
end
