using Test

include(joinpath(@__DIR__, "..", "src", "RootDepth", "module_rootdepth.jl"))

@testset "RootDepth Fortran→Julia translations" begin
    @test length(THETA_SAT) == nstyp
    @test length(KSAT) == nstyp

    @testset "init_soil" begin
        nzg = 5
        slz = zeros(Float64, nzg + 1)
        dz = zeros(Float64, nzg)
        init_soil_depth(nzg, slz, dz)
        @test slz[end] == 0.0
        @test all(dz .> 0)
        @test all(slz[1:end-1] .< slz[2:end])

        slz2 = zeros(Float64, nzg + 1)
        dz2 = zeros(Float64, nzg)
        init_soil_depth_clm(nzg, slz2, dz2)
        @test slz2[end] == 0.0
        @test all(dz2 .> 0)
    end

    @testset "interception" begin
        intercepstore, ppdrip, et_i = interception(1.0, 2.0, 3.0, 0.2, 0.0, 0.3, 0.0)
        @test isapprox(intercepstore, 0.6; atol=1e-12)
        @test isapprox(ppdrip, 1.6; atol=1e-12)
        @test et_i == 0.0
    end

    @testset "tridag" begin
        n = 3
        a = [0.0, -1.0, -1.0]
        b = [2.0, 2.0, 2.0]
        c = [-1.0, -1.0, 0.0]
        r = [0.0, 0.0, 4.0]
        u = zeros(3)
        tridag(a, b, c, r, u, n)
        @test u ≈ [1.0, 2.0, 3.0]
    end

    @testset "potevap" begin
        pet = potevap_priestly_taylor(1, 1, 300.0, 10.0, 1013.0, 0.0)
        @test pet > 0.0

        vals = potevap_shutteworth_wallace(1, 1, 3600.0, 300.0, 150.0, 100.0, 101325.0, 0.01, 2.0, 2.0, 5.0, 1.0,
                                           0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                           0.0, 0.0, 0.0, 0.0, 0)
        @test length(vals) == 12
        @test all(isfinite, vals)
    end

    @testset "update functions smoke" begin
        nzg = 3
        slz = [-1.0, -0.6, -0.2, 0.0]
        dz = [0.4, 0.4, 0.2]
        smoieq = [0.12, 0.12, 0.12]
        smoi = [0.2, 0.2, 0.2]
        soiltxt = [1, 2]

        wtd, rech = update_shallow_wtd(1, 1, nzg, 0, slz, dz, soiltxt, smoieq, 0.3, smoi, -0.7, 0.0, 2.0)
        @test isfinite(wtd)
        @test isfinite(rech)

        qlatflux = zeros(Float64, nzg + 2)
        wtd2, qspring = update_wtb_qlat(nzg, slz, dz, -0.7, 0.0, 0.01, smoi, smoieq, soiltxt, 0.3, qlatflux, 2.0)
        @test isfinite(wtd2)
        @test isfinite(qspring)
    end
end
