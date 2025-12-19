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
        maxiter=100)

    wl = init_WangLandauVars(sim.λ_max,sim.N_max,sim.L_σ)
    c = init_cache()
    return(μstate,sim,wl,c)
end

function init_cube_state()
   # initializing some vars for a test
    r_σ = [-1  1  1  -1  -1 -1  1  1;  
            -1 -1  1   1   1 -1 -1  1;
            1  1  1   1  -1 -1 -1  -1  ] # 8 particles on in cube of sidelength 2 centered at (0,0,0), lets just say "box" is sidelength 5
    T_σ=1.
    Λ_σ=argon_deBroglie(T_σ)
    L_σ=5.
    r_box = r_σ ./ L_σ
    r_box .= r_box .- round.(r_box)
    λ_max = 99
    N_max = 8
    r_frac_box = [0.,0.,0.,]
    r_cut_σ=3
    input_path = joinpath(@__DIR__, "cube_vertices_home_made.inp") 
    sim = SimulationParams(N_max=N_max,N_min=0,T_σ=T_σ,Λ_σ=Λ_σ,
                            λ_max=λ_max,r_cut_σ=r_cut_σ,
                            input_filename=input_path,# hommade input is the same as examples above 8 atoms on cube of sidelenghth 2 in simulation box of length 5 in σ units
                            save_directory_path= @__DIR__ , rng=MersenneTwister(1))
    wl = WangLandauVars(1,zeros(λ_max+1,N_max+1),zeros(λ_max+1,N_max+1),0,0,0,0,0,0.15)
    μ = microstate(size(r_box,2),34,r_box,r_frac_box)
    c = init_cache()

    return(μ,sim,wl,c)
end

μ4,sim4,wl4,c = init_4atom_state()

μcube,simcube,wlcube,c = init_cube_state()
#with @bechmark $var 'inlines' variable, not sure what that means but is supposed to be necessary for accurate benchmarking
function run_simulation_allocs_benchmark(μ::microstate,sim::SimulationParams,wl::WangLandauVars,c::SimCache)
    Profile.Allocs.clear()

    Profile.Allocs.@profile sample_rate=1 run_simulation!(sim,μ,wl,c)
    # get link to flamegraph with PProf.Allocs.pprof(from_c=false) in the REPL after this runs

end

function run_simulation_time_benchmark()
    μ4,sim4,wl4,c = init_4atom_state()
    @benchmark run_simulation!(sim4,μ4,wl4,c)
    # @time output:  16.137074 seconds (533.01 M allocations: 14.336 GiB, 5.20% gc time)
end


println("bechmarking with removed particle")
copy_microstate!(c.μ_prop,μcube)
c.μ_prop.r_box = μcube.r_box[:,2:end]
c.μ_prop.N = μcube.N-1
#@benchmark λ_metropolis_pm1($μcube,$c.μ_prop,$1,$wlcube,$simcube)

# println("benchmarking with added particle")
# copy_microstate!(c.μ_prop,μcube)
# c.μ_prop.r_box = hcat(μcube.r_box,[0. ,0., 0.])
# c.μ_prop.N = μcube.N+1
# @benchmark λ_metropolis_pm1($μcube,$c.μ_prop,$1,$wlcube,$simcube)

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