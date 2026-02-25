# this module contains system specific stuff to the argon lennard jones system
using Random
#  ✅  == checked in \test 

function argon_deBroglie(T_σ::Float64)::Float64 #  ✅ 
    # computes de broglie wavelength for the LJ model of Argon at a given input reduced temperature T_σ = kB T / ϵ

    # -- physical constants (SI)
    h_Js  = 6.62607015e-34      # Planck constant, J s
    kB_J_K = 1.380649e-23        # Boltzmann constant, J/K
    amu_kg = 1.66053906660e-27  # atomic mass unit, kg

    # -- Argon LJ parameters (typical)
    ϵ_over_kB_K = 117.05              # epsilon / k_B (K)
    σ_Å = 3.4                       # sigma in Angstroms
    mass_amu = 39.948               # atomic mass of Ar in amu

    # -- convert to SI
    T_K = ϵ_over_kB_K * T_σ             # temperature in K
    σ_m = σ_Å * 1e-10                 # sigma in meters
    m_kg= mass_amu * amu_kg             # mass in kg

    # -- de Broglie wavelength (SI)
    Λ_m = h_Js / sqrt(2*pi * m_kg * kB_J_K * T_K)

    # -- dimensionless quantities in LJ reduced units
    Λ_σ = Λ_m / σ_m


    return (Λ_σ)
end # argon debroglie

function E_12_LJ(rij_squared_σ::Float64)::Float64 #  ✅ 
    #= Computes the interaction energy between two normal lennard jones particles in LJ units 
    rij_squared_\sigma = squared distance between two particles in lennard jones \sigma =1 units
    =# 
    E_int = (1/rij_squared_σ)^6 - (1/rij_squared_σ)^3
    E_int = 4*E_int 
    return(E_int)
end 

function E_12_frac_LJ(rij_squared_σ::Float64,ϵ_ξ::Float64,σ_ξ_squared::Float64)::Float64 #  ✅ 
    #= computes interaction between fractional particle and normal LJ particle according to equation 16 of Desgranges 2016. Note that `M` in Desgranges = λ_max + 1 in our notation
    returns energy in lennard jones units
    # takes in ϵ_ξ and σ_ξ_squared instead of precomputing because they only need to be computed 1 time every time λ changes
    # assumes that rij_squared_σ was computed using the correct minimum image PBC 
    =# 

    E_int = (σ_ξ_squared/rij_squared_σ)^6 - (σ_ξ_squared/rij_squared_σ)^3
    E_int = 4*ϵ_ξ*E_int
    return(E_int)
end #E12_frac_lj 

function potential_1_normal(
    cl           ::CellList,
    r_box        ::Vector{<:AbstractVector},
    ri_box       ::AbstractVector,
    i            ::Int,
    r_frac_box   ::AbstractVector,
    λ            ::Int,
    λ_max        ::Int,
    N            ::Int,
    L_squared_σ  ::Float64,
    r_cut_sq_box ::Float64,
    ϵ_ξ          ::Float64,
    σ_ξ_sq       ::Float64,
)::Float64

    E  = 0.0
    ci = c_index(cl, ri_box)

    # ── normal–normal via cell list ──────────────────────────────────────────
    neighbours!(cl, i, ci; half=false)          # fills cl.nbr_buf[1:cl.nbr_n]
    @inbounds for k in 1:cl.nbr_n
        j          = cl.nbr_buf[k]
        rij_sq_box = euclidean_distance_squared_pbc(ri_box, r_box[j])
        rij_sq_box >= r_cut_sq_box && continue
        rij_sq_box == 0.0          && return typemax(Float64)
        E += E_12_LJ(rij_sq_box * L_squared_σ)
    end

    # ── normal–fractional (not in list, always O(1)) ─────────────────────────
    rij_sq_box = euclidean_distance_squared_pbc(ri_box, r_frac_box)
    if rij_sq_box < r_cut_sq_box
        rij_sq_box == 0.0 && return typemax(Float64)
        E += E_12_frac_LJ(rij_sq_box * L_squared_σ, ϵ_ξ, σ_ξ_sq)
    end

    return E
end

function potential_1_frac(
    cl           ::CellList,
    r_box        ::Vector{<:AbstractVector},
    r_frac_box   ::AbstractVector,
    λ            ::Int,
    λ_max        ::Int,
    N            ::Int,
    L_squared_σ  ::Float64,
    r_cut_sq_box ::Float64,
    ϵ_ξ          ::Float64,
    σ_ξ_sq       ::Float64,
)::Float64

    N == 0 && return 0.0

    E  = 0.0
    ci = c_index(cl, r_frac_box)

    neighbours!(cl, 0, ci; half=false)          # i=0 never matches a real atom
    @inbounds for k in 1:cl.nbr_n
        j          = cl.nbr_buf[k]
        rij_sq_box = euclidean_distance_squared_pbc(r_frac_box, r_box[j])
        rij_sq_box >= r_cut_sq_box && continue
        rij_sq_box == 0.0          && return typemax(Float64)
        E += E_12_frac_LJ(rij_sq_box * L_squared_σ, ϵ_ξ, σ_ξ_sq)
    end

    return E
end

function λ_metropolis_pm1(cl::CellList,                                  # ← CL
                           μ::microstate, μ_prop::microstate,
                           idx_deleted::Int64,
                           wl::WangLandauVars, sim::SimulationParams)::Bool

    # Partition-function ratio term (unchanged logic).
    logQ_diff       = wl.logQ_λN[μ.λ+1, μ.N+1] - wl.logQ_λN[μ_prop.λ+1, μ_prop.N+1]
    partition_ratio = exp(logQ_diff)

    if (μ.λ > 0 && μ_prop.λ > 0) || (μ.λ == 0 && μ_prop.λ == 0)
        V_Λ_prefactor = sim.V_σ^(μ_prop.N - μ.N) * sim.Λ_σ^(3*μ.N - 3*μ_prop.N)
    elseif μ.λ == 0 && μ_prop.λ > 0
        V_Λ_prefactor = sim.V_σ^(μ_prop.N+1 - μ.N) * sim.Λ_σ^(3*μ.N - 3*(μ_prop.N+1))
    elseif μ.λ > 0 && μ_prop.λ == 0
        V_Λ_prefactor = sim.V_σ^(μ_prop.N - μ.N - 1) * sim.Λ_σ^(3*(μ.N+1) - 3*μ_prop.N)
    end

    # Energy terms — now forwarding cl to every potential call.             # ← CL
    local E_old::Float64, E_proposed::Float64, factorial_prefactor::Float64

    if μ.N == μ_prop.N   # ── pure λ change ───────────────────────────────
        factorial_prefactor = 1.0
        E_old      = potential_1_frac(cl, μ.r_box, μ.r_frac_box,           # ← CL
                         μ.λ, sim.λ_max, μ.N,
                         sim.L_squared_σ, sim.r_cut_squared_box, μ.ϵ_ξ, μ.σ_ξ_squared)
        E_proposed = potential_1_frac(cl, μ_prop.r_box, μ_prop.r_frac_box, # ← CL
                         μ_prop.λ, sim.λ_max, μ_prop.N,
                         sim.L_squared_σ, sim.r_cut_squared_box, μ_prop.ϵ_ξ, μ_prop.σ_ξ_squared)

    elseif μ.N < μ_prop.N   # ── particle created ────────────────────────
        # Old fractional → becomes new full particle at index μ_prop.N.
        factorial_prefactor = 1.0 / μ_prop.N
        E_old      = potential_1_frac(cl, μ.r_box, μ.r_frac_box,           # ← CL
                         μ.λ, sim.λ_max, μ.N,
                         sim.L_squared_σ, sim.r_cut_squared_box, μ.ϵ_ξ, μ.σ_ξ_squared)
        # The new particle at μ_prop.N is not yet in the cell list, so we   # ← CL note
        # pass half=false (the default in potential_1_normal) which uses    # ← CL note
        # c_index on the probe position — works correctly without pre-      # ← CL note
        # inserting the particle.                                            # ← CL note
        E_proposed = potential_1_normal(cl, μ_prop.r_box,                   # ← CL
                         μ_prop.r_box[μ_prop.N], μ_prop.N,
                         μ_prop.r_frac_box, μ_prop.λ, sim.λ_max, μ_prop.N,
                         sim.L_squared_σ, sim.r_cut_squared_box,
                         μ_prop.ϵ_ξ, μ_prop.σ_ξ_squared)

    else   # μ.N > μ_prop.N  ── particle destroyed ────────────────────────
        factorial_prefactor = Float64(μ.N)
        # Compute energy of the particle-to-be-deleted using the CURRENT list
        # (it is still registered under idx_deleted, position unchanged).
        E_old      = potential_1_normal(cl, μ.r_box,                        # ← CL
                         μ.r_box[idx_deleted], idx_deleted,
                         μ.r_frac_box, μ.λ, sim.λ_max, μ.N,
                         sim.L_squared_σ, sim.r_cut_squared_box, μ.ϵ_ξ, μ.σ_ξ_squared)
        E_proposed = potential_1_frac(cl, μ_prop.r_box, μ_prop.r_frac_box, # ← CL
                         μ_prop.λ, sim.λ_max, μ_prop.N,
                         sim.L_squared_σ, sim.r_cut_squared_box,
                         μ_prop.ϵ_ξ, μ_prop.σ_ξ_squared)
    end

    ΔE         = E_proposed - E_old
    prob_ratio = partition_ratio * V_Λ_prefactor * factorial_prefactor * exp(-ΔE / sim.T_σ)

    return prob_ratio > 1.0 || prob_ratio > rand(sim.rng)
end # λ_metropolis_pm1
