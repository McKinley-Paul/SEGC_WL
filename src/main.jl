using StaticArrays

include("utils_module.jl")
include("lj_module.jl")
using .utils_module
using .lj_module

########################################### INITIALIZATION ######################################################################################

#Grand canonical/ SEGC variables
# const N_max= 4 # maximum number of atoms in the simulation
# const N_min = 0 # minimum number of atoms in simulation
# const T_σ = 1 # temperature for Q(N,V,T,λ) in lennard jones units
# # box length and volume are read in from initial configuration file

#λ = 0 # λ is what I am using for the fractional particle number, desgranges et al uses "l" for this
# λ_max = 99 # maximum value lambda can take, equivalent of M-1 in desgranges et al because λ=100 is the same as λ=0
#r_frac_box = @MArray [0.0,0.0,0.0] # location of fractional particle in box units, good use case for static array because it is a short vector, r_frac_box = [x,y,z]. but mutable so instead of @SArray we use @MArray


#= we are using allen and tildesly's initialize.py to generate our starting config
it prints the number of atoms on the first line, and the box length on the second line. 
The atomic coordinates (x,y,z) are row vectors starting on the third line.
=#

# filename = "cnf.inp"
# N,L_σ,r_σ = load_configuration(filename) # loads in N (the current number of atoms), L (the box length), and r (3xN array of particle positions) from allen and tildesly
# the σ in L_σ and r_σ denotes that these quantities are in the units from the input file, in this case the natural Lennard jones units but perhaps in the future HS units or angstroms even
# r_box = r_σ ./ L_σ # converting to "box units" by dividing the positions by the length of the box. Now r_i ∈ [-0.5,0.5]^3 because inputs are centered around 0,0
# r_box .= r_box .- round.(r_box) # periodic conditions, forces any atoms outside the box to wrap around and go back in
# min_distance,particle_1,particle_2 = min_config_distance(r_box)
# println("Minimum distance in box units (L_box=1) between particles in the initial configuration is: ",min_distance, " between particles ",particle_1, " and ", particle_2)
# println("This is compared to length of σ in box units which is: ", 1/L_σ)

#L_squared_σ = L_σ*L_σ # squared length of box in potential natural / lennard jones units, used inside potential energy evaluation
# r_cut_σ  = 3 # cutoff distance of 3σ
# r_cut_box = r_cut_σ/L_σ
# r_cut_squared_box = r_cut_box*r_cut_box

#V_σ = L_σ^3
# if N != N_max
#     throw(ArgumentError("Input mismatch: input config has $N atoms but N_max is $N_max, you want to start the simulation in the densest configuration"))
# end

# println()
# println("           Starting SEGC Wang Landau Simulation        ")
# println("T* = ", T_σ, "")
# println("V = ", V_σ, "σ^3") 
# println("Box Length L = ", L_σ, "σ")
# println("Nmax = ", N_max)
# println("Nmin = ", N_min)
# println("λmax = ", λ_max)

Λ_σ = argon_deBroglie(T_σ) # debroglie wavelength in σ LJ units but this IS SPECIFIC TO the lennard jones model of ARGON
#println("Due to the occurence of Λ in the metropolis criterion, this simulation cannot be fully carried out for lennard jonesium and is here specific to Argon")

# Wang Landau variables
logf = 1 # natural log of f, the modification factor in the Wang Landau scheme which is initiall set to f = the euler constant
H_λN = zeros(Int,λ_max,N_max) # histogram of number of times each (λ,N) pair has been visited.  λ_max rows, N_max columns.
logQ_λN = zeros(Float64,λ_max,N_max) # "density of states" in SEGC is the extended partition function Q(NVT,l) because V and T are fixed for given sim., they are left out of name, λ on rows and N on columns

# Monte Carlo variables
translation_moves_proposed = 0 
translation_moves_accepted = 0  

λ_moves_proposed = 0 
λ_moves_accepted = 0

########################################### MONTE CARLO ######################################################################################

f_convergence_threshold = 10^(-8)
logf_convergence_threshold = log(f_convergence_threshold)
δr_max_box = 0.15/L_σ # inital max distance to move LJ particle in translation moves. Will be dynamically adjusted in the loop. Using the same initial value as Allen and Tildesly's default. Note that this is δr_max per xyz component of the move

while logf ≥ logf_convergence_threshold

    ######################### Proposing and Accepting a Move #########################

    ζ = rand() # \zeta is a random Float64 uniformly sampled from [0,1)

    if ζ < 0.75 # propose translational moves 75% of the time 
        translation_moves_proposed += 1

        i = rand(1:(N+1)) # randomly pick atom to move, includes fractional particle
        if i < N # move normal particle
            ri_box = @view r_box[:,i]
        
            ri_proposed_box = translate_by_random_vector(ri_box,δr_max_box) # Trial move to new position (in box=1 units)
            ri_proposed_box .= ri_proposed_box .- round.(ri_proposed_box)   # PBC
            E_proposed =potential_1_normal(r_box,ri_proposed_box,i,r_frac_box,λ,λ_max,N,L_squared_σ,r_cut_squared_box) 

            if E_proposed == typemax(Float64) # check overlap
                nothing # reject the move and recount this state for histogram and partition function 

            else   

                 E_old =potential_1_normal(r_box,ri_box,i,r_frac_box,λ,λ_max,N,L_squared_σ,r_cut_squared_box) 
                 ΔE = E_proposed - E_old 
                 accept = metropolis(ΔE,T_σ) 
                 if accept
                    r_box[:,i] = ri_proposed_box
                    translation_moves_accepted += 1
                 end
            end

        else # move fractional particle
            r_frac_proposed_box = translate_by_random_vector(r_frac_box,δr_max_box)
            r_frac_proposed_box .= r_frac_proposed_box .- round.(r_frac_proposed_box)
            if λ == 0 # auto accept because if λ =0 , translational move of the fractional particle doesn't change energy
                r_frac_box = r_frac_proposed_box
                translation_moves_accepted += 1
            else 
                E_proposed = potential_1_frac(r_box, r_frac_proposed_box ,λ,λ_max,N,L_squared_σ,r_cut_squared_box)

                E_old =potential_1_frac(r_box,r_frac_box,λ,λ_max,N,L_squared_σ,r_cut_squared_box)
                ΔE = E_proposed - E_old 
                    accept = metropolis(ΔE,T_σ) 
                    if accept
                        r_frac_box = r_frac_proposed_box
                        translation_moves_accepted += 1
                    end
            end

        end # if i < N deciding to move normal or translational particle

        if translation_moves_accepted/translation_moves_proposed > 0.55 # tune δr_max_box to get ~50% acceptance, pg 159 Allen Tildesly
            δr_max_box = δr_max_box * 1.05
        elseif translation_moves_accepted/translation_moves_proposed  < 0.45
            δr_max_box = δr_max_box*0.95
        end 

    else # this means  ζ ≥ 0.75 and so we propose λ moves 25% of the time 
        # currently only implementing Δλ = ±1, can do CFMC/translational move style drawing from scaled uniform distribution if this doesnt work well
        λ_moves_proposed += 1

        λ_proposed = λ + 2*rand(Bool) - 1 # change λ by ±1; can go out of our range, and can be -1 or 100 which is our signal to change N 

        if -1 < λ_proposed ≤ λ_max # < 99 for λ_max = 99 no change to N
            N_proposed = N 
            r_proposed_box = @view r_box[:,:]
            r_frac_proposed_box  = @view r_frac_box[:,:]
            idx_deleted = 0
        elseif λ_proposed == -1 # decrement N 
            # keep the fractional particle in r_frac_box but delete a random particle from the list, need to pick a random particle I think
            λ_proposed = λ_max     
            N_proposed = N-1
            
            idx_deleted = rand(1:N)
            cols = [1:idx_deleted-1; idx_deleted+1:size(r_box,2)] # columns (particles) to keep 
            r_proposed_box = @view r_box[:,cols]
            r_frac_proposed_box  = @view r_frac_box[:,:]
            

        elseif λ_proposed ≥ λ_max # increment N 
            λ_proposed = 0 
            N_proposed = N+1

            r_proposed_box = hcat(r_box, r_frac_box) # add the fractional particle to the end of the list
            r_frac_proposed_box = @MArray (rand(3) .- 0.5)
            idx_deleted = 0
        else
            throw("Error in λ move ΔN control flow")
        end # ΔN control flow

        if N_min ≤ N_proposed ≤ N_max # reject if N goes out of bounds for the sim
            accept = false
        else 
            accept = λ_metropolis_pm1(λ,N,r_box,r_frac_box,
                              λ_proposed, N_proposed, r_proposed_box, r_frac_proposed_box,idx_deleted,
                              logQ_λN, Λ_σ,V_σ,T_σ,
                              λ_max,L_squared_σ,r_cut_squared_box)
            if accept == true
                λ_moves_accepted += 1
                λ = λ_proposed
                N = N_proposed
                r_box = r_proposed_box
                r_frac_box = r_frac_proposed_box
            end
        end
    end # proposal if statement based on ζ

    #### update Wang Landau Variables
    logQ_λN[λ,N] = logQ_λN[λ,N] + logf
    H_λN[λ,N] += 1
    
    if (translation_moves_proposed + λ_moves_proposed) % 1000 ==0 # check flatness only every 1,000 cycles
        min = minimum(H_λN)
        if min ≥  1000 # this is our flatness criteria
            H_λN = zeros(Int,λ_max,N_max)
            logf = 0.5*logf
        end
    end
    break
end

println("Wang Landau converged, logf has reached ", logf, " and convergence is achieved when logf reaches ", logf_convergence_threshold)
println("Total translation moves proposed: ", translation_moves_proposed, ", translation moves accepted: ", translation_moves_accepted, ", Acceptance ratio: ", translation_moves_accepted/translation_moves_proposed)
println("Total λ moves proposed: ", λ_moves_proposed, ", λ moves accepted: ", λ_moves_accepted, ", Acceptance ratio: ", λ_moves_accepted/λ_moves_proposed)
println("Partition function logQ_λN:")
println(logQ_λN)
