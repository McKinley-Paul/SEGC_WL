include("../src/initialization_module.jl")
using initialization_module

input_path = "./4_atoms_cnf.inp"
μstate = init_microstate(input_path)
T_σ = 1 
Λ_σ = argon_deBroglie(T_σ)
sim = SimulationParams(
    N_max=4,
    N_min=0,
    T_σ=T_σ,
    Λ_σ = Λ_σ,
    λ_max = 99,
    r_cut_σ = 3,
    input_filename=input_path)

wl = init_WangLandauVars(sim.λ_max,sim.N_max,sim.L_σ)

initialization_check(sim,μstate,wl)

run_simulation!(sim,μstate,wl)

post_run(sim,μstate,wl)