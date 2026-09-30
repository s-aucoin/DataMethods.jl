@testset "SignalProcessing" begin
    @testset "Welch and FT_params" begin
        segs = [collect(0.0:0.1:0.9), collect(0.1:0.1:1.0)]
        I, F = Welch(segs, 10)

        @test size(I) == (10,)
        @test size(F) == (10, 2)

        params = FT_params(10.0, 0.5)
        @test params.N == 5
        @test length(params.freqs) == 5
        @test params.freqs[3] ≈ 0.0
    end

    @testset "Spike detection" begin
        signal = [0.0, 0.0, 0.0, 100.0, 0.0, 0.0]
        spikes = findspikes(signal; sm=3, std_thresh=1.5)
        ref = DataMethods.moving_average(signal, 3)
        despiked = despike(copy(signal); sm=3, std_thresh=1.5)

        @test spikes[4] == true
        @test despiked[spikes] == ref[spikes]
        @test despiked[.!spikes] == signal[.!spikes]
    end
end
