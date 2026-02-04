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
    λ_max = 9,
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
# avg seconds 0.34566888300000004 
# avg logQ [0.0, 14.029298728704454, 26.978791612386704, 39.389910909533505, 52.09403860867024]
# top value is pretty significantly different 52 here for M =10 vs 42 for M = 100 ...
# avg iters: 2,119,600
# avg translatio 0.962837162
# avg lambda: 0.905370343 