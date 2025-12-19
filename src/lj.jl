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

end #E12_LJ

function E_12_frac_LJ(rij_squared_σ::Float64,λ::Int64,λ_max::Int64)::Float64 #  ✅ 
    #= computes interaction between fractional particle and normal LJ particle according to equation 16 of Desgranges 2016. Note that `M` in Desgranges = λ_max + 1 in our notation
    returns energy in lennard jones units
    # assumes that rij_squared_σ was computed using the correct minimum image PBC 
    =# 
    rij_σ = sqrt(rij_squared_σ) # sqrt is unavoidable here I think 
    M = λ_max + 1
    ϵ_ξ = (λ/M)^(1/3) # fractional energy  coupling parameter in lennard jones units (so really this might be named ϵ_ξ_σ but that just looked like too much)
    σ_ξ = (λ/M)^(1/4) # \sigma_\xi is the fractional distance coupling parameter in LJ σ units  
    E_int = (σ_ξ/rij_σ)^12 - (σ_ξ/rij_σ)^6
    E_int = 4*ϵ_ξ*E_int
    return(E_int)

end #E12_frac_lj 


function potential_1_normal(r_box::Vector{SVector{3,Float64}},ri_box::SVector{3,Float64},i::Int64,r_frac_box::SVector{3,Float64},λ::Int64,λ_max::Int64,N::Int64,L_squared_σ::Float64,r_cut_squared_box::Float64)::Float64 #  ✅ 
    #= Calculates the sum of the interaction potential between 1 special particle labelled i (e.g. the particle involved in a proposal move) in the r list and all others including the fractional particle
    Inputs:
        -r_box (Array of Float64s): 3xN position matrix of N normal particles in box units
        -ri_box (static array of Float64s): 3x1 vector of the position of ri you want to use for distances between all other particles in box units
        -i (int): index of particle in r list you want to compute interactions between, used for avoiding self interactions.
        - r_frac_box (static array of Float64s): location of fractional particle in box units
        - λ (int): coupling constant 
        - λ_max (int): M-1 in desgranges paper, number of unique values λ can take 
        - N (int): current number of particles, could be read from _,N = size(r_box) but want to minimize all activity in this loop and we already have to keep track of N outside this function call
        - L_squared_σ (Float64): squared length of box in σ units (LJ this time)
        - r_cut_squared_box  (Float64): cutoff length squared  for the potential in box=1 units
    Returns:
        - E_int_σ (Float64): total interaction energy between particle i at location ri and all other particles including the fractional one in natural units for the potential (i.e. lennard jones units)
    =# 

    E_int_σ = 0

    # first computing interaction with normal particles
    for j in 1:N
        if j != i # avoid double counting
            rj_box = r_box[j]
            rij_squared_box = euclidean_distance_squared_pbc(ri_box,rj_box)
            if rij_squared_box < r_cut_squared_box # only evaluate potential if pair is inside cutoff distance
                rij_squared_σ = rij_squared_box  * L_squared_σ   # convert to potential natural units (i.e. lennard jones) to compute the potentials see allen and Tildesly one Note note for explanation why
                if rij_squared_σ >  0.5 # close to overlap threshold given by allen and tildesly of potential E > 100, ours is E > 224
                    E_int_σ = E_int_σ + E_12_LJ(rij_squared_σ)
                else # if particles overlap quit early and return a huge energy
                    return(typemax(Float64))
                end 
            end
        end
    end 

    # adding fractional interaction
    rij_squared_box = euclidean_distance_squared_pbc(ri_box,r_frac_box)
    if (rij_squared_box < r_cut_squared_box) && (rij_squared_box != 0.0)# if the distance is zero, E_12_frac_LJ returns NaN
        rij_squared_σ = rij_squared_box  * L_squared_σ  
        E_int_σ = E_int_σ + E_12_frac_LJ(rij_squared_σ,λ,λ_max)
    elseif rij_squared_box == 0.0 # if this is the case need to just reject this config because LJ involves 1/rij which if rij=0 is undefined
            return(typemax(Float64))
    end
    
    return(E_int_σ)

end # potential_1_normal

function potential_1_frac(r_box::Vector{SVector{3,Float64}},r_frac_box::SVector{3,Float64},λ::Int64,λ_max::Int64,N::Int64,L_squared_σ::Float64,r_cut_squared_box::Float64)::Float64
    #= this function calculates the sum of all interactions between the fractional particle and all the normal particles
        -r_box (Array of Float64s): 3xN position matrix of N normal particles in box units
        - r_frac_box (static array of Float64s): location of fractional particle in box units
        - λ (int): coupling constant 
        - λ_max (int): M-1 in desgranges paper, number of unique values λ can take 
        - N (int): current number of particles, could be read from _,N = size(r_box) but want to minimize all activity in this loop and we already have to keep track of N outside this function call
        - L_squared_σ (Float64): squared length of box in σ units (LJ this time)
        - r_cut_squared_box  (Float64): cutoff length squared  for the potential in box=1 units
    Returns:
        - E_int_σ (Float64): total interaction energy between particle i at location ri and all other particles including the fractional one in natural units for the potential (i.e. lennard jones units)
    =# 
    E_int_σ = 0

    for j in 1:N
        rij_squared_box = euclidean_distance_squared_pbc(r_frac_box,r_box[j])
        if (rij_squared_box < r_cut_squared_box) && (rij_squared_box > 0) # if rij_squared_box ==0 E_12_frac_LJ returns a NaN
            rij_squared_σ = rij_squared_box  * L_squared_σ 
            # here we would check for overlap normally and reject but we will not do this for the fractional particle because if λ=1 or something, could still have low energy if overlapped
            E_int_σ = E_int_σ + E_12_frac_LJ(rij_squared_σ,λ,λ_max)
        elseif rij_squared_box == 0 # if this is the case need to just reject (or get out of) this config because LJ involves 1/rij which if rij=0 is undefined/+infinity
            return(typemax(Float64))
        end
    end
    return(E_int_σ)

end # potential_1_frac


@inline function  λ_metropolis_pm1(μ::microstate,
                           μ_prop::microstate, idx_deleted::Int64, # μ_prop = μ_proposed the next proposed microstate
                           wl::WangLandauVars,sim::SimulationParams)::Bool #  ✅ 

        # MAKES SOME SIMPLIFYING ASSUMPTIONS THAT ONLY WORK WHEN LAMBDA CAN ONLY CHANGE BY ±1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        # the purpose of the function is to handle the complex control flow based on the current value of λ and the proposed one.
        # many things change including the form of the metropolis criterion, what energies you have to compute to get ΔE, and so on.
        # covered by in general by equations 10-12 in Desgranges 2012
        # assumes that you have already checked that N_proposed is in bounds  (N_min ≤ N_proposed ≤ N_max)

        # first we compute the multiplicative prefactor term involving Q,V,Λ in eqns 10-12-- because the N! terms can cause overflow, we only compute them once we know something about N old vs new
        logQ_diff = wl.logQ_λN[μ.λ+1 , μ.N+1] - wl.logQ_λN[μ_prop.λ+1 , μ_prop.N+1]
        partition_ratio = exp(logQ_diff)
        if (μ.λ > 0 && μ_prop.λ > 0) || (μ.λ == 0 && μ_prop.λ == 0)
            V_Λ_prefactor = sim.V_σ^(μ_prop.N -μ.N) * sim.Λ_σ^(3*μ.N - 3*μ_prop.N)
        elseif μ.λ ==0 && μ_prop.λ > 0
            V_Λ_prefactor = sim.V_σ^(μ_prop.N+1 - μ.N) * sim.Λ_σ^(3*μ.N - 3*(μ_prop.N+1))
        elseif μ.λ > 0 && μ_prop.λ == 0 
            V_Λ_prefactor = sim.V_σ^(μ_prop.N - μ.N - 1) * sim.Λ_σ^(3*(μ.N+1)-3*μ_prop.N)
        end

        # now we compute the exponential part of the criterion having to do with the configurational potential energy
        # and the N-dependent factorial prefactor factorial_prefactor = Nold!/Nnew!

        if μ.N == μ_prop.N # λ changed so change in configurational energy only has to do with fractional particle old vs new
            # following equation 10 desgranges
           E_old = potential_1_frac(μ.r_box,μ.r_frac_box,   μ.λ   ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box) 
           E_proposed = potential_1_frac(μ_prop.r_box,μ_prop.r_frac_box,      μ_prop.λ    ,sim.λ_max,μ_prop.N,sim.L_squared_σ,sim.r_cut_squared_box)
           factorial_prefactor = 1
        
        elseif μ.N < μ_prop.N # particle created so λ = λ_max and λ_proposed = 0
            #  change in potential energy comes from fractional particle becoming full particle and the new fractional particle with λ=0 makes no contribution to energy
            # so E_Old = E_old_frac interaction with all others and E_new = E_new_full_particle interaction with all others
            factorial_prefactor = 1/μ_prop.N # only works for ±1
            E_old = potential_1_frac(μ.r_box,μ.r_frac_box,  μ.λ   ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box)
            i = length(μ_prop.r_box)
            E_proposed = potential_1_normal(μ_prop.r_box , μ_prop.r_box[i],i,μ_prop.r_frac_box,μ_prop.λ,sim.λ_max,μ_prop.N,sim.L_squared_σ,sim.r_cut_squared_box)
                    
        elseif μ.N > μ_prop.N # particle destroyed so λ = 0 and λ_proposed = 99
            # old energy is energy of destroyed particle with rest of full particles 
            # new configurational energy is interaction of fractional particle with others
            factorial_prefactor = μ.N # only works for ±1
            # destroyed_particle = @view μ.r_box[:,idx_deleted]
            E_old = potential_1_normal(μ.r_box,μ.r_box[idx_deleted],idx_deleted,μ.r_frac_box,μ.λ,sim.λ_max,μ.N,sim.L_squared_σ,sim.r_cut_squared_box)
            E_proposed = potential_1_frac(μ_prop.r_box,μ_prop.r_frac_box,  μ_prop.λ   ,sim.λ_max,μ_prop.N ,sim.L_squared_σ,sim.r_cut_squared_box)
        end # ΔN logic 

        ΔE = E_proposed - E_old
        exponent = -1*ΔE/sim.T_σ 
        prob_ratio = partition_ratio*V_Λ_prefactor*factorial_prefactor*exp(exponent)
        if prob_ratio > 1
            return(true)
        else
            ζ = rand(sim.rng)
            accept = (prob_ratio > ζ)   #boolean
            return(accept)
        end
end #λ_metropolis_pm1
