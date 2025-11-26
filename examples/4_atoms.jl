include("../src/initialization_module.jl")
using initialization_module

input_path = "./4_atoms_cnf.inp"
μstate = init_microstate(input_path)
Λ_σ = 3.0
sim = SimulationParams(
    N_max=4,
    N_min=0,
    T_σ=1,
    Λ_σ = Λ_σ,
    λ_max = 99,
    r_cut_σ = 3,
    input_filename=input_path)

if Λ_σ == 3.0
    throw("you need to use debroglie wavelength fr not just the placeholder you put in originally")
end