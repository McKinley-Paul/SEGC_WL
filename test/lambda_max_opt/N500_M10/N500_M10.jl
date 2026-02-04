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
N_min=0,
T_σ=T_σ,
Λ_σ = Λ_σ,
λ_max = 9,
r_cut_σ = 3.,
input_filename=input_path,
save_directory_path= @__DIR__ , 
maxiter=100_000_000_000_000 ) # 100 trilli iters

μstate = init_microstate(sim=sim,filename=input_path)

wl = init_WangLandauVars(sim.λ_max,sim.N_max,sim.L_σ)

cache = init_cache(sim,μstate)

initialization_check(sim,μstate,wl)

seconds = @elapsed run_simulation!(sim,μstate,wl,cache) # start 11:15pm


#global seconds_avg += seconds
post_run(sim,μstate,wl)

logQ = correct_Q(wl)
#     for ii in 1:length(logQ)
#         logQ_avg[ii] += logQ[ii]
#     end
# end 
# logQ_avg = logQ_avg.*(1/5) # average value of partition functions over five monte carlo runs
# seconds_avg =  seconds_avg*(1/5)
println(seconds) # seconds was 
println(logQ)
# THIS RAN FOR OVER 16 HOURS AND DID NOT FINISH ONE EPOCH OF WANG LANDAU
