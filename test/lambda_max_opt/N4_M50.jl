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
    λ_max = 49,
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
# average time: 1.4816433414
# log Q average: [0.0, 13.931493091583253, 26.825919499993326, 31.961642715334893, 44.716592928767206]
# qualitatively similar to M100 values but a little different
# avg iters: 18,075,600
# avg translation acceptance rate: 0.959512476
# avg lambda acceptance rate: 0.943231223