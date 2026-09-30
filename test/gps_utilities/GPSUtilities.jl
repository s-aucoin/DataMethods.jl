using CSV
using DataFrames
using Dates
using TimeZones

@testset "GPSUtilities" begin
    @testset "time alignment helpers" begin
        t0 = ZonedDateTime(DateTime(2020, 1, 1, 0, 0, 0), tz"UTC")
        gps_time = [t0, t0, t0 + Second(1), t0 + Second(1)]
        with_ms = add_ms(gps_time; fs_gps=2)

        @test with_ms[1] == t0
        @test with_ms[2] == t0 + Millisecond(500)
        @test with_ms[3] == t0 + Second(1)
        @test with_ms[4] == t0 + Second(1) + Millisecond(500)

        times = collect(t0:Second(1):(t0 + 4*Second(1)))
        idx = index_by_time(times, [t0 + Second(1), t0 + Second(3)], [t0 + Second(2), t0 + Second(4)])
        @test idx == [2:3, 4:5]
    end

    @testset "CSV cleanup and averaging" begin
        mktempdir() do dir
            gps1 = DataFrame(
                "DATE" => ["200101", "200101"],
                "TIME" => ["000000", "000001"],
                "LATITUDE N/S" => ["47.0N", "47.2N"],
                "LONGITUDE E/W" => ["10.0W", "10.2W"],
            )
            gps2 = DataFrame(
                "DATE" => ["200101", "200101"],
                "TIME" => ["000000", "000001"],
                "LATITUDE N/S" => ["47.1N", "47.3N"],
                "LONGITUDE E/W" => ["10.1W", "10.3W"],
            )

            p1 = joinpath(dir, "gps1.csv")
            p2 = joinpath(dir, "gps2.csv")
            CSV.write(p1, gps1)
            CSV.write(p2, gps2)

            cleaned1_pth = splitext(p1)[1] * "_cleaned.csv"
            clean_gps_data(p1; fs_gps=2, datatypes=[String, String, String, String])
            cleaned1 = CSV.read(cleaned1_pth, DataFrame)

            cleaned2_pth = splitext(p2)[1] * "_cleaned.csv"
            clean_gps_data(p2; fs_gps=2, datatypes=[String, String, String, String])
            cleaned2 = CSV.read(cleaned2_pth, DataFrame)
            @test nrow(cleaned1) >= 2
            @test "time" in names(cleaned1)

            mean_gps_positions([cleaned1_pth, cleaned2_pth])
            mean_file = joinpath(dir, "GPS_mean.csv")
            @test isfile(mean_file)

            mean_df = CSV.read(mean_file, DataFrame)
            @test "lat" in names(mean_df)
            @test "lon" in names(mean_df)
            @test nrow(mean_df) >= 1

            clean_and_mean([p1, p2]; fs_gps=2, datatypes=[String, String, String, String])
            @test isfile(joinpath(dir, "GPS_mean.csv"))
        end
    end
end
