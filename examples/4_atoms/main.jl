using Test
using StaticArrays
using Random

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using segc_wl   # or the module name inside segc_wl.jl

input_path = joinpath(@__DIR__, "4_atom_cnf.inp")
μstate = init_microstate(filename=input_path)
T_σ = 1.0 
Λ_σ = argon_deBroglie(T_σ)
println(Λ_σ)
sim = SimulationParams(
    N_max=4,
    N_min=0,
    T_σ=T_σ,
    Λ_σ = Λ_σ,
    λ_max = 99,
    r_cut_σ = 3.,
    input_filename=input_path,
    save_directory_path= @__DIR__ , 
    maxiter=100_000_000)

wl = init_WangLandauVars(sim.λ_max,sim.N_max,sim.L_σ)

cache = init_cache(sim,μstate)

initialization_check(sim,μstate,wl)

@time run_simulation!(sim,μstate,wl,cache)

post_run(sim,μstate,wl)

logQ = correct_Q(wl)

println(logQ)