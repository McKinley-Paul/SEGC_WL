using Test
using StaticArrays
using Random

import Pkg
Pkg.activate(joinpath(@__DIR__, "../.."))
using segc_wl   

using Revise
using BenchmarkTools
using Profile
using PProf

# see this video for how to use REPL: https://www.youtube.com/watch?v=mZOYM48YEFg ctrl+L to clear REPL

function init_4atom_state()
    input_path = joinpath(@__DIR__, "4_atom_cnf.inp")
    μstate = init_microstate(filename=input_path)
    T_σ = 1.0 
    Λ_σ = argon_deBroglie(T_σ)
    sim = SimulationParams(
        N_max=4,
        N_min=0,
        T_σ=T_σ,
        Λ_σ = Λ_σ,
        λ_max = 99,
        r_cut_σ = 3.,
        input_filename=input_path,
        save_directory_path= @__DIR__ , 
        maxiter=10_000)

    wl = init_WangLandauVars(sim.λ_max,sim.N_max,sim.L_σ)
    c = init_cache(sim,μstate)
    return(μstate,sim,wl,c)
end

function init_cube_state()
   # initializing some vars for a test
    # r_σ = [-1  1  1  -1  -1 -1  1  1;  
    #         -1 -1  1   1   1 -1 -1  1;
    #         1  1  1   1  -1 -1 -1  -1  ] # 8 particles on in cube of sidelength 2 centered at (0,0,0), lets just say "box" is sidelength 5
    # r_box = r_σ ./ L_σ
    # r_box .= r_box .- round.(r_box)
    T_σ=1.
    Λ_σ=argon_deBroglie(T_σ)
    λ_max = 99
    N_max = 8
    r_cut_σ=3
    input_path = joinpath(@__DIR__, "cube_vertices_home_made.inp") 
    sim = SimulationParams(N_max=N_max,N_min=0,T_σ=T_σ,Λ_σ=Λ_σ,
                            λ_max=λ_max,r_cut_σ=r_cut_σ,
                            input_filename=input_path,# hommade input is the same as examples above 8 atoms on cube of sidelenghth 2 in simulation box of length 5 in σ units
                            save_directory_path= @__DIR__ , rng=MersenneTwister(1),maxiter=100_000)
    wl = WangLandauVars(1,zeros(λ_max+1,N_max+1),zeros(λ_max+1,N_max+1),0,0,0,0,0,0.15)
    μ = init_microstate(filename=input_path)
    c = init_cache(sim,μ)
    return(μ,sim,wl,c)
end
μ4,sim4,wl4,c4 = init_4atom_state()
μcube,simcube,wlcube,ccube = init_cube_state()
#with @bechmark $var 'inlines' variable, not sure what that means but is supposed to be necessary for accurate benchmarking
function run_simulation_allocs_benchmark(μ::microstate,sim::SimulationParams,wl::WangLandauVars,c::SimCache)
    Profile.Allocs.clear()

    Profile.Allocs.@profile sample_rate=1 run_simulation!(sim,μ,wl,c)
    # get link to flamegraph with PProf.Allocs.pprof(from_c=false) in the REPL after this runs

end

#run_simulation_allocs_benchmark(μcube,simcube,wlcube,ccube)

function run_simulation_time_benchmark()
    μ4,sim4,wl4,c = init_4atom_state()
    @benchmark run_simulation!(sim4,μ4,wl4,c)
    # @time output:  16.137074 seconds (533.01 M allocations: 14.336 GiB, 5.20% gc time)
end




############### warntype #########################
function warntype_translation_move(μ::microstate,sim::SimulationParams,wl::WangLandauVars,c::SimCache)
    @code_warntype translation_move!(sim,μ,wl) # this looks pretty good and only shows red for AbstractRNG, which i fixed
end

function warntype_λ_move(μ::microstate,sim::SimulationParams,wl::WangLandauVars,c::SimCache)
    @code_warntype λ_move!(sim,μ,wl,c) # nothing in red but some unions of concrete types highlighted in yellow
end

function warntype_run_simulation(μ::microstate,sim::SimulationParams,wl::WangLandauVars)
    @code_warntype run_simulation!(sim,μ,wl) # nothing in red or yellow
end

simcube.maxiter


function E_12_LJ_old(rij_squared_σ::Float64)::Float64 #  ✅ 
    #= Computes the interaction energy between two normal lennard jones particles in LJ units 
    rij_squared_\sigma = squared distance between two particles in lennard jones \sigma =1 units
    =# 
    E_int = (1/rij_squared_σ)^6 - (1/rij_squared_σ)^3
    E_int = 4*E_int 
    return(E_int)

end #E12_LJ

 function E_12_LJ_new(rij_squared_σ::Float64)::Float64
    inv_r6 = (inv(rij_squared_σ))^3 # inv(x) = 1/x for x::Float64
    return 4.0 * (inv_r6^2 - inv_r6)
end
 
