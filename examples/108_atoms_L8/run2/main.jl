using Test
using StaticArrays
using Random

import Pkg
Pkg.activate("/Users/mckinleypaul/Documents/montecarlo/segc_wl/src/")
Pkg.activate("/Users/mckinleypaul/Documents/montecarlo/segc_wl")
using segc_wl   # or the module name inside segc_wl.jl

input_path =  "/Users/mckinleypaul/Documents/montecarlo/segc_wl/initial_configs/N108_L8.inp"
T_σ = 1.0 
Λ_σ = argon_deBroglie(T_σ)
println(Λ_σ)
sim = SimulationParams(
    N_max=108,
    N_min=0,
    T_σ=T_σ,
    Λ_σ = Λ_σ,
    λ_max = 99,
    r_cut_σ = 2.5, # 3 sigma for replication of results, partition function is similar for 2.5 and 3 sigma
    input_filename=input_path,
    save_directory_path= @__DIR__ , 
    maxiter=100_000_000_000_000) # 100 trillion iters 

μstate = init_microstate(sim=sim,filename=input_path)

wl = init_WangLandauVars(sim)

cache = init_cache(sim,μstate)

initialization_check(sim,μstate,wl)

cl = CellList(sim.N_max + 1, sim.r_cut_box)
make_list!(cl, μstate.N, μstate.r_box)

run_simulation!(sim,cl,μstate,wl,cache)

post_run(sim,μstate,wl)

logQ = correct_Q(wl)

println(logQ)