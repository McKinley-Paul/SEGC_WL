using Test
using StaticArrays
using Random

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using segc_wl   # or the module name inside segc_wl.jl

input_path = joinpath(@__DIR__, "../../initial_configs/4_atom_cnf.inp")
T_σ = 1.0 
Λ_σ = argon_deBroglie(T_σ)
println(Λ_σ)


#initialization_check(sim,μstate,wl)
logQ_avg = [0.,0.,0.,0.,0.]
seconds_avg = 0
for _ in 1:5
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

    μstate = init_microstate(sim=sim,filename=input_path)

    wl = init_WangLandauVars(sim.λ_max,sim.N_max,sim.L_σ)

    cache = init_cache(sim,μstate)

    seconds = @elapsed run_simulation!(sim,μstate,wl,cache)


    global seconds_avg += seconds
    post_run(sim,μstate,wl)
    logQ = correct_Q(wl)
    for ii in 1:length(logQ)
        logQ_avg[ii] += logQ[ii]
    end
end 
logQ_avg = logQ_avg.*(1/5) # average value of partition functions over five monte carlo runs
seconds_avg =  seconds_avg*(1/5)
println(seconds_avg)
println(logQ_avg)
# 3.8866949751999997 seconds was the average
# [0.0, 13.889340269565583, 25.2001935005188, 29.506129717826845, 42.1611690312624]
# average total iterations: 51,937,000
# average translation acceptance ratio: 0.958823438
# average lambda acceptance ratio (calculated by hand in excel): 0.954698777