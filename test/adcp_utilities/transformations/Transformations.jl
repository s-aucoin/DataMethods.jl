@testset "Transformations" begin
    data = reshape([1.0, 2.0, 3.0, 4.0, 5.0], :, 1)
    result = rm_sidelobe_contam(data, collect(0.0:1.0:4.0), 2.5)

    @test size(result) == size(data)
    @test any(isnan, result)
end
