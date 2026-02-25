module segc_wl
using StaticArrays
using Dates # mostly for debugging/monitoring long calculations
include("initialization.jl")
include("cell_list.jl")
include("utils.jl")
include("lj.jl")
include("thermo.jl")

# ✅  = unit tests for function exist in /tests/
# this module contains the meat of the package including run_simulation and the various High Level Wang Landau Monte Carlo functions. 


# segc_wl.jl functions:
export run_simulation!, translation_move!,λ_move!,update_wl!, post_run
# initialization functions exports:
export microstate,SimulationParams,init_microstate,check_inputs, print_simulation_params, print_microstate,print_wl, WangLandauVars,init_WangLandauVars, initialization_check, save_wanglandau_jld2, save_microstate_jld2,load_microstate_jld2, load_wanglandau_jld2, load_configuration, SimCache,init_cache,copy_microstate!
# utils exports:
export euclidean_distance, min_config_distance, euclidean_distance_squared_pbc, translate_by_random_vector, metropolis
# lj exports:
export argon_deBroglie, E_12_LJ, E_12_frac_LJ, potential_1_normal, potential_1_frac, λ_metropolis_pm1
# thermo stuff:
export correct_Q
#CellList:
export CellList, make_list!, neighbours!

function run_simulation!(sim::SimulationParams, cl::CellList,             # ← CL
                         μ::microstate, wl::WangLandauVars, c::SimCache)

    log_path     = joinpath(sim.save_directory_path, "wl_progress_log.txt")
    progress_log = open(log_path, "a")

    f_convergence_threshold    = exp(10^(-8))
    logf_convergence_threshold = log(f_convergence_threshold)

    while (wl.logf ≥ logf_convergence_threshold) && (wl.iters < sim.maxiter)

        if rand(sim.rng) < 0.75
            translation_move!(sim, cl, μ, wl, c)                          # ← CL
        else
            λ_move!(sim, cl, μ, wl, c)                                    # ← CL
        end

        update_wl!(wl, μ)

        if wl.iters % 1000 == 0
            min_H = minimum(wl.H_λN[:, (sim.N_min+1):(sim.N_max+1)])
            if min_H ≥ 1000
                wl.H_λN  = zeros(Int64, sim.λ_max+1, sim.N_max+1)
                wl.logf  = 0.5 * wl.logf
                println(progress_log, "New WL epoch! logf = ", wl.logf,
                        " at ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
                flush(progress_log)
            end

            if wl.iters % 100_000 == 0
                save_wanglandau_jld2(wl, sim, "wl_checkpoint.jld2")
                save_microstate_jld2(μ,  sim, "microstate_checkpoint.jld2")
            end
        end

    end # while

    close(progress_log)
end # run_simulation!

function translation_move!(sim::SimulationParams, cl::CellList,          # ← CL
                           μ::microstate, wl::WangLandauVars, c::SimCache)

    wl.translation_moves_proposed += 1

    c.ζ_idx = rand(sim.rng, 1:(μ.N + 1))   # includes fractional particle slot

    if c.ζ_idx <= μ.N   # ── move a normal particle ──────────────────────────
        # NB: original had `c.ζ_idx < μ.N` which accidentally excluded the
        # last normal particle; corrected to `<= μ.N` here.

        c.ri_proposed_box .= μ.r_box[c.ζ_idx]
        translate_by_random_vector!(c.ri_proposed_box, wl.δr_max_box, sim.rng, c)
        pbc_wrap!(c.ri_proposed_box)

        E_proposed = potential_1_normal(cl, μ.r_box, c.ri_proposed_box, c.ζ_idx,  # ← CL
                         μ.r_frac_box, μ.λ, sim.λ_max, μ.N,
                         sim.L_squared_σ, sim.r_cut_squared_box, μ.ϵ_ξ, μ.σ_ξ_squared)

        if E_proposed == typemax(Float64)
            nothing   # reject: overlap
        else
            E_old = potential_1_normal(cl, μ.r_box, μ.r_box[c.ζ_idx], c.ζ_idx,   # ← CL
                         μ.r_frac_box, μ.λ, sim.λ_max, μ.N,
                         sim.L_squared_σ, sim.r_cut_squared_box, μ.ϵ_ξ, μ.σ_ξ_squared)

            ΔE = E_proposed - E_old
            if metropolis(ΔE, sim.T_σ, sim.rng)
                # 1. compute new cell BEFORE writing position                 # ← CL
                ci_new = c_index(cl, c.ri_proposed_box)                      # ← CL

                # 2. update positions
                μ.r_box[c.ζ_idx]        .= c.ri_proposed_box
                c.μ_prop.r_box[c.ζ_idx] .= c.ri_proposed_box

                # 3. update cell list — move atom to its new cell             # ← CL
                move_in_list!(cl, c.ζ_idx, ci_new)                           # ← CL

                wl.translation_moves_accepted += 1
            end
        end

    else   # ── move the fractional particle ───────────────────────────────
        # Fractional particle is NEVER in the cell list, so no list update needed.

        c.ri_proposed_box .= μ.r_frac_box
        translate_by_random_vector!(c.ri_proposed_box, wl.δr_max_box, sim.rng, c)
        pbc_wrap!(c.ri_proposed_box)

        if μ.λ == 0
            # λ=0 → fractional particle decoupled → auto-accept
            μ.r_frac_box            .= c.ri_proposed_box
            c.μ_prop.r_frac_box     .= c.ri_proposed_box
            wl.translation_moves_accepted += 1
        else
            E_proposed = potential_1_frac(cl, μ.r_box, c.ri_proposed_box,    # ← CL
                             μ.λ, sim.λ_max, μ.N,
                             sim.L_squared_σ, sim.r_cut_squared_box, μ.ϵ_ξ, μ.σ_ξ_squared)

            E_old = potential_1_frac(cl, μ.r_box, μ.r_frac_box,              # ← CL
                        μ.λ, sim.λ_max, μ.N,
                        sim.L_squared_σ, sim.r_cut_squared_box, μ.ϵ_ξ, μ.σ_ξ_squared)

            ΔE = E_proposed - E_old
            if metropolis(ΔE, sim.T_σ, sim.rng)
                μ.r_frac_box        .= c.ri_proposed_box
                c.μ_prop.r_frac_box .= c.ri_proposed_box
                wl.translation_moves_accepted += 1
            end
        end
    end

    # ── dynamic step-size tuning (unchanged) ────────────────────────────────
    if sim.dynamic_δr_max_box && wl.logf == 1
        ratio = wl.translation_moves_accepted / wl.translation_moves_proposed
        if ratio > 0.55 && wl.δr_max_box < 1.0
            wl.δr_max_box *= 1.05
        elseif ratio < 0.45
            wl.δr_max_box *= 0.95
        end
    end
end # translation_move!

function λ_move!(sim::SimulationParams, cl::CellList,                    # ← CL
                 μ::microstate, wl::WangLandauVars, c::SimCache)

    wl.λ_moves_proposed += 1
    λ_proposed = μ.λ + 2*rand(sim.rng, Bool) - 1   # ±1

    # ── Case 1: pure λ change, N unchanged ──────────────────────────────────
    if -1 < λ_proposed ≤ sim.λ_max

        c.μ_prop.λ            = λ_proposed
        c.μ_prop.ϵ_ξ          = (λ_proposed / (sim.λ_max + 1))^(1/3)
        c.μ_prop.σ_ξ_squared  = (λ_proposed / (sim.λ_max + 1))^(1/2)

        # No particle positions change → no cell-list update needed.         # ← CL note
        if λ_metropolis_pm1(cl, μ, c.μ_prop, 0, wl, sim)                   # ← CL
            wl.λ_moves_accepted += 1
            μ.λ           = c.μ_prop.λ
            μ.ϵ_ξ         = c.μ_prop.ϵ_ξ
            μ.σ_ξ_squared = c.μ_prop.σ_ξ_squared
        else
            c.μ_prop.λ           = μ.λ
            c.μ_prop.ϵ_ξ         = μ.ϵ_ξ
            c.μ_prop.σ_ξ_squared = μ.σ_ξ_squared
        end

    # ── Case 2: λ_proposed == -1  →  destroy a normal particle, N decreases ─
    elseif λ_proposed == -1

        c.μ_prop.λ = sim.λ_max
        c.μ_prop.N = μ.N - 1

        if c.μ_prop.N < sim.N_min
            c.μ_prop.λ = μ.λ
            c.μ_prop.N += 1
            return nothing
        end

        c.μ_prop.ϵ_ξ         = (sim.λ_max / (sim.λ_max + 1))^(1/3)
        c.μ_prop.σ_ξ_squared = (sim.λ_max / (sim.λ_max + 1))^(1/2)

        idx_deleted = rand(sim.rng, 1:μ.N)

        # Proposed state: swap last particle into deleted slot (in cache only).
        if idx_deleted != μ.N
            c.μ_prop.r_box[idx_deleted] .= μ.r_box[μ.N]
        end

        if λ_metropolis_pm1(cl, μ, c.μ_prop, idx_deleted, wl, sim)         # ← CL
            wl.λ_moves_accepted += 1

            # ── Apply deletion to the cell list ─────────────────────────────
            # Step A: remove idx_deleted from the list.                       # ← CL
            ci_del = SVector{3,Int}(cl.cell[:, idx_deleted])                 # ← CL
            destroy_in_list!(cl, idx_deleted, ci_del)                        # ← CL

            if idx_deleted != μ.N
                # Step B: particle μ.N is being moved into slot idx_deleted.
                # Remove μ.N from the list under its old index …             # ← CL
                ci_last = SVector{3,Int}(cl.cell[:, μ.N])                   # ← CL
                destroy_in_list!(cl, μ.N, ci_last)                           # ← CL
                # … then re-insert it under the new index idx_deleted.       # ← CL
                create_in_list!(cl, idx_deleted, ci_last)                    # ← CL
                # (ci_last is still correct because the position didn't move) # ← CL

                # Update actual microstate position array.
                μ.r_box[idx_deleted] .= μ.r_box[μ.N]
            end

            μ.N           -= 1
            μ.λ            = c.μ_prop.λ
            μ.ϵ_ξ          = c.μ_prop.ϵ_ξ
            μ.σ_ξ_squared  = c.μ_prop.σ_ξ_squared

            # Cell-list entry for index μ.N+1 is now stale but unreachable   # ← CL
            # because make_list! / all loops only go up to the current μ.N.  # ← CL

        else   # rejected
            c.μ_prop.λ           = μ.λ
            c.μ_prop.N          += 1
            c.μ_prop.r_box[idx_deleted] .= μ.r_box[idx_deleted]
            c.μ_prop.ϵ_ξ         = μ.ϵ_ξ
            c.μ_prop.σ_ξ_squared = μ.σ_ξ_squared
        end

    # ── Case 3: λ_proposed > λ_max  →  create a normal particle, N increases ─
    elseif λ_proposed > sim.λ_max

        c.μ_prop.λ = 0
        c.μ_prop.N = μ.N + 1

        if c.μ_prop.N > sim.N_max
            c.μ_prop.λ = μ.λ
            c.μ_prop.N -= 1
            return nothing
        end

        # The old fractional particle becomes the new full particle at index N+1.
        c.μ_prop.r_box[c.μ_prop.N] .= μ.r_frac_box

        # Draw a new random position for the fractional particle.
        @inbounds for i in 1:3
            c.μ_prop.r_frac_box[i] = rand(sim.rng) - 0.5
        end

        c.μ_prop.ϵ_ξ         = (0 / (sim.λ_max + 1))^(1/3)   # = 0.0
        c.μ_prop.σ_ξ_squared = (0 / (sim.λ_max + 1))^(1/2)   # = 0.0

        if λ_metropolis_pm1(cl, μ, c.μ_prop, 0, wl, sim)                   # ← CL
            wl.λ_moves_accepted += 1

            μ.N          += 1
            μ.λ           = 0
            μ.r_box[μ.N] .= c.μ_prop.r_box[μ.N]   # materialise new particle
            μ.r_frac_box .= c.μ_prop.r_frac_box

            # Register the new particle in the cell list.                    # ← CL
            ci_new = c_index(cl, μ.r_box[μ.N])                              # ← CL
            create_in_list!(cl, μ.N, ci_new)                                 # ← CL

            μ.ϵ_ξ         = c.μ_prop.ϵ_ξ
            μ.σ_ξ_squared = c.μ_prop.σ_ξ_squared

        else   # rejected
            c.μ_prop.λ          = μ.λ
            c.μ_prop.N         -= 1
            c.μ_prop.r_frac_box .= μ.r_frac_box
            c.μ_prop.ϵ_ξ        = μ.ϵ_ξ
            c.μ_prop.σ_ξ_squared = μ.σ_ξ_squared
        end

    else
        error("λ_move!: unexpected λ_proposed = $λ_proposed")
    end
end # λ_move!

function update_wl!(wl::WangLandauVars,μ::microstate)
    wl.logQ_λN[μ.λ+1,μ.N+1] += wl.logf
    wl.H_λN[μ.λ+1,μ.N+1] += 1
    wl.iters +=1
end #update_wl!

function post_run(sim::SimulationParams,μ::microstate,wl::WangLandauVars)
    println("Wang Landau converged or reached max iterations, logf has reached ", wl.logf, " and convergence is achieved when logf reaches ", log( exp(10^(-8))) )
    println("Iterations: ", (wl.translation_moves_proposed+wl.λ_moves_proposed), " with maxiters: ", sim.maxiter )
    println("Total translation moves proposed: ", wl.translation_moves_proposed, ", translation moves accepted: ", wl.translation_moves_accepted, ", Acceptance ratio: ", wl.translation_moves_accepted/wl.translation_moves_proposed)
    println("Total λ moves proposed: ", wl.λ_moves_proposed, ", λ moves accepted: ", wl.λ_moves_accepted, ", Acceptance ratio: ", wl.λ_moves_accepted/wl.λ_moves_proposed)
    save_wanglandau_jld2(wl,sim,"final_wl.jld2")
    save_microstate_jld2(μ,sim,"final_microstate.jld2")
    println("Final Microstate and Wang Landau variables saved")
end

end