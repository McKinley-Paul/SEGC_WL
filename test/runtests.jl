using Test
using StaticArrays
using Random

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using segc_wl   # or the module name inside segc_wl.jl
using segc_wl.utils_module
using segc_wl.initialization_module
using segc_wl.lj_module


@testset "Utils" begin
    r1 = 100*rand(3)
    r2 = 100*rand(3)
    @test euclidean_distance(r1,r2) ≈ sqrt(sum((r1.-r2).^2))
    r1_box = [0,0,0.4] 
    r2_box = [0,0,-0.4]
    @test euclidean_distance_squared_pbc(r1_box,r2_box) ≈ 0.04

    r = [0.0,0.0,0.0]
    r1 = translate_by_random_vector(r,0.1)
    @test sqrt(sum((r1).^2)) ≤ sqrt(3)*0.1 # test that it moves it by no more than \sqrt(3)*δr_max_box from sqrt(0.1^2 + 0.1^2 + 0.1^2)

    @test metropolis(100.0,1.0) == false
    @test metropolis(-0.1,1.0) == true
    @test metropolis(0.0,1.0) == true

end

@testset "Metropolis acceptance rates" begin
    rng = MersenneTwister(1234)
    
    # For ΔE/T = 1, acceptance should be ≈ exp(-1) ≈ 36.8%
    ΔE = 1.0
    T_σ = 1.0
    n_trials = 10000
    n_accepted = sum(metropolis(ΔE, T_σ, rng) for _ in 1:n_trials)
    acceptance_rate = n_accepted / n_trials
    expected_rate = exp(-ΔE / T_σ)
    
    # Test within reasonable statistical bounds (±3σ for binomial)
    #σ = sqrt(n*p*(1-p)) / n
    σ = sqrt(n_trials * expected_rate * (1 - expected_rate)) / n_trials
    # @test acceptance_rate ≈ expected_rate atol=3*σ
    @test acceptance_rate ≈ expected_rate atol = σ
    println("Regular metropolis testing: ΔE/T=1: Expected $(expected_rate*100)% acceptance rate, Got $(acceptance_rate*100)%")
end