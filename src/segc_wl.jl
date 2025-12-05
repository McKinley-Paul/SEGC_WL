module segc_wl
using StaticArrays
include("initialization_module.jl")
using .initialization_module

include("utils_module.jl") # even though only initialization_module is used here, i had to do these to export them in this file so they are public facing and can be tested in /test/
using .utils_module
include("lj_module.jl")
using .lj_module
# ✅  = unit tests for function exist in /tests/
# this module contains the meat of the package including run_simulation and the various
# High Level Wang Landau Monte Carlo functions. 
export run_simulation!, translation_move!,λ_move!,update_wl!, post_run

export 
    utils_module,
    initialization_module,
    lj_module



function run_simulation!(sim::SimulationParams, μ::microstate,wl::WangLandauVars)
    f_convergence_threshold = 10^(-8)
    logf_convergence_threshold = log(f_convergence_threshold)
    while (wl.logf ≥ logf_convergence_threshold ) && ( (wl.λ_moves_proposed + wl.translation_moves_proposed) < sim.maxiter)
        ζ = rand(sim.rng)
        if ζ < 0.75                         # propose translational moves 75% of the time 
            translation_move!(sim,μ,wl)
        else                                # this means  ζ ≥ 0.75 and so we propose λ moves 25% of the time 
            λ_move!(sim,μ,wl)
        end

        update_wl!(wl,μ)

        # check flatness and save state:
        if (wl.λ_moves_proposed + wl.translation_moves_proposed) % 1000 ==0 # check flatness only every 1000 cycles
            min = minimum(wl.H_λN)
            if min ≥  1000 # this is our flatness criteria
                #save_wanglandau_jld2(wl,sim, "wl_checkpoint_before_rezeroing.jld2") # jld2 is quick save binary, to inspect checkpoint, open up julia ipynb and use wl_loaded = load_wanglandau_jld2("checkpoint.jld2")
                wl.H_λN = zeros(Int64,sim.λ_max+1,sim.N_max+1)
                wl.logf = 0.5*wl.logf
            end
            if (wl.λ_moves_proposed + wl.translation_moves_proposed) % 50000 ==0 # save checkpoint every 50,000 moves 
                save_wanglandau_jld2(wl,sim, "wl_checkpoint.jld2") # jld2 is quick save binary, to inspect checkpoint, open up julia ipynb and use wl_loaded = load_wanglandau_jld2("checkpoint.jld2")
                save_microstate_jld2(μ,sim, "microstate_checkpoing.jld2")
            end
        end #flatness/printing

    end # while logf ≥ logf_convergence_threshold

    #return(wl) is this necessary?
end #run_simulation

function translation_move!(sim::SimulationParams,μ::microstate,wl::WangLandauVars) #✅
    # wrote the body of all these functions before refactoring into structs, could rewrite at some point because ugly
    wl.translation_moves_proposed += 1

    i = rand(sim.rng,1:(μ.N+1)) # randomly pick atom to move, includes fractional particle
    if i < μ.N # move normal particle
        ri_box = @view μ.r_box[:,i]
    
        ri_proposed_box = translate_by_random_vector(ri_box, wl.δr_max_box, sim.rng) # Trial move to new position (in box=1 units)
        ri_proposed_box .= ri_proposed_box .- round.(ri_proposed_box)   # PBC
        E_proposed =potential_1_normal(μ.r_box,ri_proposed_box,i,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box) 

        if E_proposed == typemax(Float64) # check overlap
            nothing # reject the move and recount this state for histogram and partition function 
        else   
            E_old =potential_1_normal(μ.r_box,ri_box,i,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box) 
            ΔE = E_proposed - E_old 
            accept = metropolis(ΔE,sim.T_σ,sim.rng) 
            if accept
                μ.r_box[:,i] = ri_proposed_box
                wl.translation_moves_accepted += 1
            end
        end

    else # move fractional particle
        r_frac_proposed_box = translate_by_random_vector(μ.r_frac_box,wl.δr_max_box,sim.rng)
        r_frac_proposed_box .= r_frac_proposed_box .- round.(r_frac_proposed_box)
        if μ.λ == 0 # auto accept because if λ =0 , translational move of the fractional particle doesn't change energy
            μ.r_frac_box = r_frac_proposed_box
            wl.translation_moves_accepted += 1
        else 
            E_proposed = potential_1_frac(μ.r_box, r_frac_proposed_box  ,μ.λ ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box)

            E_old =potential_1_frac(μ.r_box,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box)
            ΔE = E_proposed - E_old 
                accept = metropolis(ΔE,sim.T_σ,sim.rng) 
                if accept
                    μ.r_frac_box = r_frac_proposed_box
                    wl.translation_moves_accepted += 1
                end
        end

    end # if i < N deciding to move normal or translational particle

    if (wl.translation_moves_accepted/wl.translation_moves_proposed > 0.55) && wl.δr_max_box < 1.0 # tune δr_max_box to get ~50% acceptance, pg 159 Allen Tildesly
        # added the wl.δr_max_box < 1.0 because for dilute systems or ideal gas conditions you accept every move and the δr_max_box grows riducously and unphysically for a periodic system using  box=1 units
        wl.δr_max_box = wl.δr_max_box * 1.05

    elseif wl.translation_moves_accepted/wl.translation_moves_proposed  < 0.45
        wl.δr_max_box = wl.δr_max_box*0.95
    end 
end# translation move

function λ_move!(sim::SimulationParams,μ::microstate,wl::WangLandauVars)
    # again wrote the body of this and all subfunctions before doing structure refactor, could rewrite because as of now it's ugly and verbose but if it works...
    # currently only implementing Δλ = ±1, can do CFMC/translational move style drawing from scaled uniform distribution if this doesnt work well
    wl.λ_moves_proposed += 1

    λ_proposed = μ.λ + 2*rand(sim.rng,Bool) - 1 # change λ by ±1; can go out of our range, and can be -1 or 100 which is our signal to change N 

    if -1 < λ_proposed ≤ sim.λ_max # λ_proposed < 99 for λ_max = 99 no change to N
        N_proposed = μ.N 
        r_proposed_box = @view μ.r_box[:,:]
        r_frac_proposed_box  = @view μ.r_frac_box[:,:]
        idx_deleted = 0
    elseif λ_proposed == -1 # decrement N 
        # keep the fractional particle in r_frac_box but delete a random particle from the list, need to pick a random particle I think
        λ_proposed = sim.λ_max     
        N_proposed = μ.N-1
        if N_proposed < sim.N_min # return early/reject move if negative particles appear i.e. if λ_proposed = -1 initially when N=0
            # have to do this here otherwise line below idx_deleted will throw error when μ.N == 0
            return()
        end
        
        idx_deleted = rand(sim.rng,1:μ.N)
        cols = [1:idx_deleted-1; idx_deleted+1:size(μ.r_box,2)] # columns (particles) to keep 
        r_proposed_box = @view μ.r_box[:,cols]
        r_frac_proposed_box  = @view μ.r_frac_box[:,:]
        

    elseif λ_proposed ≥ sim.λ_max # increment N 
        λ_proposed = 0 
        N_proposed = μ.N+1
        if N_proposed > sim.N_max # return early/reject move if we go out of bounds with too many particles (i.e. λ_prop = 100 when N=N_max)
            return()
        end
        r_proposed_box = hcat(μ.r_box, μ.r_frac_box) # add the fractional particle to the end of the list
        r_frac_proposed_box = @MArray rand(sim.rng, 3) # 3 componenets betwewen [0,1)
        r_frac_proposed_box .-= 0.5 # shift so three components in correct domain [-0.5,0.5]
        idx_deleted = 0
    else
        throw("Error in λ move ΔN control flow")
    end # ΔN control flow

    accept = λ_metropolis_pm1(μ.λ,μ.N,μ.r_box,μ.r_frac_box,
                            λ_proposed, N_proposed, r_proposed_box, r_frac_proposed_box,idx_deleted,
                            wl.logQ_λN, sim.Λ_σ,sim.V_σ,sim.T_σ,
                            sim.λ_max,sim.L_squared_σ,sim.r_cut_squared_box,sim.rng)
    if accept == true
        wl.λ_moves_accepted += 1
        μ.λ = λ_proposed
        μ.N = N_proposed
        μ.r_box = r_proposed_box
        μ.r_frac_box = r_frac_proposed_box
    end
end

function update_wl!(wl::WangLandauVars,μ::microstate)
    wl.logQ_λN[μ.λ+1,μ.N+1] += wl.logf
    wl.H_λN[μ.λ+1,μ.N+1] += 1
end #update_wl!

function post_run(sim::SimulationParams,μ::microstate,wl::WangLandauVars)
    println("Wang Landau converged or reached max iterations, logf has reached ", wl.logf, " and convergence is achieved when logf reaches ", log( 10^(-8) ) )
    println("Iterations: ", (wl.translation_moves_proposed+wl.λ_moves_proposed), " with maxiters: ", sim.maxiter )
    println("Total translation moves proposed: ", wl.translation_moves_proposed, ", translation moves accepted: ", wl.translation_moves_accepted, ", Acceptance ratio: ", wl.translation_moves_accepted/wl.translation_moves_proposed)
    println("Total λ moves proposed: ", wl.λ_moves_proposed, ", λ moves accepted: ", wl.λ_moves_accepted, ", Acceptance ratio: ", wl.λ_moves_accepted/wl.λ_moves_proposed)
    save_wanglandau_jld2(wl,sim,"final_wl.jld2")
    save_microstate_jld2(μ,sim,"final_microstate.jld2")
    println("Final Microstate and Wang Landau variables saved")
end

end