using Test
using StaticArrays
using Random
import Pkg
Pkg.activate("/Users/mckinleypaul/Documents/montecarlo/segc_wl/src/")
Pkg.activate("/Users/mckinleypaul/Documents/montecarlo/segc_wl")
using segc_wl 

input_path = "/Users/mckinleypaul/Documents/montecarlo/segc_wl/initial_configs/N350_L8_random_deletion.inp"
T_σ = 1.2 
Λ_σ = argon_deBroglie(T_σ)
println(Λ_σ)
sim = SimulationParams(
    N_max=350,
    N_min=0,
    T_σ=T_σ,
    Λ_σ = Λ_σ,
    λ_max = 99,
    r_cut_σ = 3.,
    input_filename=input_path,
    save_directory_path= @__DIR__ , 
    maxiter=100_000_000_000_000) # 100 trillion iters

μstate = load_microstate_jld2("/Users/mckinleypaul/Documents/montecarlo/segc_wl/examples/N350_L8_T1.2/run1/microstate_checkpoint.jld2")

wl = load_wanglandau_jld2("/Users/mckinleypaul/Documents/montecarlo/segc_wl/examples/N350_L8_T1.2/run1/wl_checkpoint_before_rezeroing.jld2")

cache = init_cache(sim,μstate)

# initialization_check(sim,μstate,wl)
flush(stdout)
run_simulation!(sim,μstate,wl,cache)

post_run(sim,μstate,wl)

logQ = correct_logQ(wl)

println(logQ)