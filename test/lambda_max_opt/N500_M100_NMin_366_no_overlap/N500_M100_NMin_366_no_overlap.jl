using Test
using StaticArrays
using Random

import Pkg
Pkg.activate(joinpath(@__DIR__, "../../../"))
using segc_wl   # or the module name inside segc_wl.jl

input_path = joinpath(@__DIR__, "../../../initial_configs/N500_L8.inp")
T_σ = 1.0 
Λ_σ = argon_deBroglie(T_σ)
println(Λ_σ)

# logQ_avg = [0.,0.,0.,0.,0.]
# seconds_avg = 0
# for _ in 1:5
sim = SimulationParams(
N_max=500,
N_min=366,
T_σ=T_σ,
Λ_σ = Λ_σ,
λ_max = 99,
r_cut_σ = 3.,
input_filename=input_path,
save_directory_path= @__DIR__ , 
maxiter=100_000_000_000_000 ) # 100 trilli iters

μstate = init_microstate(sim=sim,filename=input_path)

wl = init_WangLandauVars(sim.λ_max,sim.N_max,sim.L_σ)

cache = init_cache(sim,μstate)

initialization_check(sim,μstate,wl)

seconds = @elapsed run_simulation!(sim,μstate,wl,cache) # started at 7:07pm 

post_run(sim,μstate,wl)

logQ = correct_Q(wl)

println(seconds) 
println(logQ)
