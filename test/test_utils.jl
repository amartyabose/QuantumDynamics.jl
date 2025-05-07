@testitem "Matrix2Vector and Back" begin
    for _ = 1:10
        random_matrix = rand(ComplexF64, 5, 5)
        @test random_matrix == Utilities.density_matrix_vector_to_matrix(Utilities.density_matrix_to_vector(random_matrix))
    end
end

@testitem "Trapezoid" begin
    x = 0:1e-6:π
    sin_x = sin.(x)
    cos_x = cos.(x)
    sin2_x = sin.(x) .^ 2
    cos2_x = cos.(x) .^ 2
    @testset "Single-threaded" begin
        @test abs(Utilities.trapezoid(x, sin_x) - 2.0) < 1e-5
        @test abs(Utilities.trapezoid(x, cos_x)) < 1e-5
        @test abs(Utilities.trapezoid(x, sin2_x) - π / 2) < 1e-5
        @test abs(Utilities.trapezoid(x, cos2_x) - π / 2) < 1e-5
    end
    @testset "Multithread" begin
        @test abs(Utilities.trapezoid(x, sin_x; exec=FLoops.ThreadedEx()) - 2.0) < 1e-5
        @test abs(Utilities.trapezoid(x, cos_x; exec=FLoops.ThreadedEx())) < 1e-5
        @test abs(Utilities.trapezoid(x, sin2_x; exec=FLoops.ThreadedEx()) - π / 2) < 1e-5
        @test abs(Utilities.trapezoid(x, cos2_x; exec=FLoops.ThreadedEx()) - π / 2) < 1e-5
    end
end

@testitem "Fourier transform" begin
    ts = 0:0.01:100
    f = sin.(2ts) .* exp.(-ts .^ 2 ./ 50)
    ωs, F = Utilities.fourier_transform(ts, f; neg=false)
    t, ftilde = Utilities.inverse_fourier_transform(ωs, F; x0=ts[1], neg=false)
    @test sum(abs2.(ts .- t)) ≤ 1e-15
    @test sum(abs2.(f .- ftilde)) ≤ 1e-15
end