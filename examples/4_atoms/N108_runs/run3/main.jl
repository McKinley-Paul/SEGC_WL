using Test
using StaticArrays
using Random
using Dates

import Pkg
Pkg.activate("/Users/mckinleypaul/Documents/montecarlo/segc_wl/src/")
Pkg.activate("/Users/mckinleypaul/Documents/montecarlo/segc_wl")

using segc_wl   # or the module name inside segc_wl.jl

input_path = "/Users/mckinleypaul/Documents/montecarlo/segc_wl/initial_configs/N108_DD_Density.inp"
T_σ = 1.0 
Λ_σ = argon_deBroglie(T_σ)

sim = SimulationParams(
        N_max=108,
        N_min=0,
        T_σ=T_σ,
        Λ_σ = Λ_σ,
        λ_max = 99,
        r_cut_σ = 3.,
        input_filename=input_path,
        save_directory_path= @__DIR__ , 
        maxiter=100_000_000_000_000,# 100 trillion iters
        dynamic_δr_max_box  = true)

μstate = init_microstate(sim=sim,filename=input_path)

wl = init_WangLandauVars(sim) 

cache = init_cache(sim,μstate)

initialization_check(sim,μstate,wl)

# starting at 
println("Starting run now, time is ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) 
flush(stdout)

run_simulation!(sim,μstate,wl,cache)

# stopping at 
println("Finished run now, time is ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS")) 
flush(stdout)

post_run(sim,μstate,wl)

logQ = correct_Q(wl)

println(logQ)
