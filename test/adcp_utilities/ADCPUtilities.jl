@testset "ADCPUtilities" begin
    mktempdir() do dir
        files = [
            joinpath(dir, "sample_001.nc"),
            joinpath(dir, "sample_002.nc"),
            joinpath(dir, "sample.ad2cp.nc")
        ]
        for f in files
            write(f, "placeholder")
        end

        names = get_sig_fnames(dir, "sample")
        @test names[end] == files[3]
        @test names[1:end-1] == sort(names[1:end-1])
    end
end
