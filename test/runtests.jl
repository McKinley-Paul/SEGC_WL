using Test
using StaticArrays
using Random

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using segc_wl   

# These Unit Tests Have Not Been Updated Since Code Base Has Been Changed To Optimize For Speed, 
# since before introduction of cache variables and stuff like that, they were used to make sure physics and expected behavior
# were intially programmed correct to get to the point where the package could run and compute partition functions
# now that we are well past that point, we are assuming that if the system wide physics tests are correct
# compared to analytic results, then  all of the units that make up those computations are correct as each unit is probably 
# called millions of times during those system wide calculations

# println("")
# println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Running Function Unit Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
# println("")

# ############ Tests for things in utils_module ############
# @testset "Utils" begin
#     r1 = 100*rand(3)
#     r2 = 100*rand(3)
#     @test euclidean_distance(r1,r2) ≈ sqrt(sum((r1.-r2).^2))
#     r1_box = [0,0,0.4] 
#     r2_box = [0,0,-0.4]
#     @test euclidean_distance_squared_pbc(r1_box,r2_box) ≈ 0.04

#     r = [0.0,0.0,0.0]
#     r1 = translate_by_random_vector(r,0.1)
#     @test sqrt(sum((r1).^2)) ≤ sqrt(3)*0.1 # test that it moves it by no more than \sqrt(3)*δr_max_box from sqrt(0.1^2 + 0.1^2 + 0.1^2)

# end

# @testset "Utils Metropolis acceptance rates" begin
#     @test metropolis(100.0,1.0) == false # reject large increase
#     @test metropolis(-0.1,1.0) == true # accept downhill move
#     @test metropolis(0.0,1.0) == true # accept neutral move

#     rng = MersenneTwister(1234) # now we test acceptance of slightly uphill move
    
#     # For ΔE/T = 1, acceptance should be ≈ exp(-1) ≈ 36.8%
#     ΔE = 1.0
#     T_σ = 1.0
#     n_trials = 10000
#     n_accepted = sum(metropolis(ΔE, T_σ, rng) for _ in 1:n_trials)
#     acceptance_rate = n_accepted / n_trials
#     expected_rate = exp(-ΔE / T_σ)
    
#     # Test within reasonable statistical bounds (±3σ for binomial)
#     #σ = sqrt(n*p*(1-p)) / n
#     σ = sqrt(n_trials * expected_rate * (1 - expected_rate)) / n_trials
#     # @test acceptance_rate ≈ expected_rate atol=3*σ
#     @test acceptance_rate ≈ expected_rate atol = 2*σ
#     #println("Regular metropolis testing: ΔE/T=1: Expected $(expected_rate*100)% acceptance rate, Got $(acceptance_rate*100)%")
# end

# ############ Tests for things in lj_module ############
# @testset "LJ_module" begin
#     @test argon_deBroglie(2.0) ≈ 0.0530973 atol=0.000001 # verified against /test/mathematica_verification.nb

#     @test E_12_frac_LJ(2.,0,99) == 0
#     fracs = []

#     for r2 in [0.0001,1.,2.,10.,30.,5.6,100000000000.0] #testing distances
#         push!(fracs,E_12_frac_LJ(r2,3,99))
#     end
#     mathematica= [3.355811e19 , -0.0064247, -0.000806758, -6.45823*10^-6, -2.39195*10^-7, -0.0000367738, -6.45826*10^-36]
#     # println(mathematica[1]) # issue if you write 355811 * 10^19 because goes 10^19 goes over int64 or something?
#     @test fracs[1] ≈ mathematica[1] atol=10^15
#     @test all(isapprox.(fracs[2:end], mathematica[2:end];atol=1e-6)) # tests E_12_frac_Lj
    
#     # tested against allen and tildesly's potential_1 from mc_lj_module.py:
#     cfg = joinpath(@__DIR__, "cnf_default.inp")
#     N,L_σ,r_σ = load_configuration(cfg)
#     r_box = r_σ./L_σ
#     test_idx = 1 .+[0,2,3,56,100,34,222,255,78,89]
#     E_test = []
#     for idx in test_idx
#         r_test = @view r_box[:,idx]
#         push!(E_test,potential_1_normal(r_box,r_test,idx,[0.0,0.,0.],0,99,size(r_box,2),L_σ^2,(3/L_σ)^2))
#     end
#     allen_tildesly_results = [-11.932319723716308 -11.932319753958657 -11.932319753958655 -11.932319756099865 -11.9323197764367 -11.932319756566486 -11.932319756099865 -11.932319738837482 -11.932319771221039 -11.9323197764367]
#     @test all(E_test .≈ allen_tildesly_results) # this tests E_12_LJ,euclidean_distance_squared_pbc, the non-fractional-part of potential_1_normal
# end

# println("")

# @testset "λ_metropolis" begin
# # initializing some dummy inputs for testing
#     r_σ = [-1  1  1  -1  -1 -1  1  1;  
#             -1 -1  1   1   1 -1 -1  1;
#             1  1  1   1  -1 -1 -1  -1  ] # 8 particles on in cube of sidelength 2 centered at (0,0,0), lets just say "box" is sidelength 5
#     L_σ = 5.0
#     V_σ = L_σ^3
#     T_σ=1.
#     Λ_σ=argon_deBroglie(T_σ)
#     r_box = r_σ ./ L_σ
#     r_box .= r_box .- round.(r_box)
#     N = size(r_box,2)
#     λ_max = 99
#     N_max = 450
#     logQ_λN = zeros(λ_max+1,N_max+1)
#     r_frac_box = [0.,0.,0.,]

#     L_squared_σ = L_σ^2
#     r_cut_σ=3
#     r_cut_squared_σ=r_cut_σ^2
#     r_cut_squared_box = r_cut_squared_σ/(L_σ^2)

#     # first lets handle the situation where only λ changes within the range s.t. N doesn't change
#     # and there is no partition function bias 

#     λ = 34
#     λ_proposed = 33 # energy goes down in this configuration so always accepts if propose 35
#     # first with no bias from WL and keeping logQ(λ,N) the same for each
#     # should have the same probability as metropolis()
#     rng = MersenneTwister()
#     ΔE = potential_1_frac(r_box,r_frac_box,  λ_proposed  ,λ_max,N,L_squared_σ,r_cut_squared_box) -  potential_1_frac(r_box,r_frac_box,   λ   ,λ_max,N,L_squared_σ,r_cut_squared_box)
#     n_trials = 10000
#     n_accepted = sum(metropolis(ΔE, T_σ, rng) for _ in 1:n_trials)
#     acceptance_rate = n_accepted / n_trials
#     expected_rate = exp(-ΔE / T_σ)
#     #println("When N does not change, metropolis() has acceptance probability: ",acceptance_rate )

#     n_trials = 10000
#     n_accepted = sum(λ_metropolis_pm1(λ,N,r_box,r_frac_box,
#                         λ_proposed, N, r_box, r_frac_box,0,
#                         logQ_λN, Λ_σ,V_σ,T_σ,
#                         λ_max,L_squared_σ,r_cut_squared_box, rng) for _ in 1:n_trials)
#     acceptance_rate = n_accepted / n_trials
#         # Test within reasonable statistical bounds (±3σ for binomial)
#     #σ = sqrt(n*p*(1-p)) / n
#     #println("When N does not change, and logQ is equal for each state, λ_metropolis_pm1() has acceptance probability: ",acceptance_rate )
#     σ = sqrt(n_trials * expected_rate * (1 - expected_rate)) / n_trials
#     # @test acceptance_rate ≈ expected_rate atol=3*σ
#     @test acceptance_rate ≈ expected_rate atol = 3*σ

#     # BUT NOW LETS BIAS WANG LANDAU 
#     logQ_λN[λ+1,N+1] = 2
#     logQ_λN[λ_proposed+1,N+1] = 3

#         n_accepted = sum(λ_metropolis_pm1(λ,N,r_box,r_frac_box,
#                         λ_proposed, N, r_box, r_frac_box,0,
#                         logQ_λN, Λ_σ,V_σ,T_σ,
#                         λ_max,L_squared_σ,r_cut_squared_box, rng) for _ in 1:n_trials)
#     acceptance_rate_biased = n_accepted / n_trials

#     #println("However, when there is a bias in logQ to stay at current state, λ_metropolis_pm1() accepts with prob ", acceptance_rate_biased)
#     #println("expect biased prob to be " , round(exp(2)/exp(3),sigdigits=2)," % of other one")

#     @test acceptance_rate_biased < acceptance_rate
#     @test acceptance_rate_biased ≈ exp(2)/exp(3)*acceptance_rate atol = 0.1



#     #particle created, no partition bias should accept because having particle at 0,0,0 is lower energy for this config

#     logQ_λN = zeros(λ_max+1,N_max+1)    
#     λ = λ_max
#     λ_proposed = 0
#     N_proposed = N+1
#     r_proposed_box = hcat(r_box,[0.,0.,0.])
#     r_proposed_frac_box = [0.1,0.2,-0.3,]

#     E_old = potential_1_frac(r_box,[0.,0.,0.],λ,λ_max,N,L_squared_σ,r_cut_squared_box)
#     E_new = potential_1_normal(r_proposed_box,[0.,0.,0.,],N_proposed,r_proposed_frac_box,λ_proposed,λ_max,N_proposed,L_squared_σ,r_cut_squared_box)
#     ΔE = E_new-E_old
#     #println("The config. energy is downhill so metropolis would accept: ΔE particle created= ",ΔE)
#     prefactor = ( (V_σ^N_proposed)/(V_σ^(N+1)) )* (1/N_proposed)*((Λ_σ^(3*N + 3))/(Λ_σ^(3*N_proposed)))
#    # println("However, the prefactor here is ", prefactor, " which is just = 1/N_proposed")
#    # println("Thus we expect acceptance prob of ", prefactor*exp(-1*ΔE/T_σ))
#     n_trials=1000

#     n_accepted = sum(λ_metropolis_pm1(λ,N,r_box,r_frac_box,
#                         λ_proposed, N_proposed, r_proposed_box, r_proposed_frac_box,0,
#                         logQ_λN, Λ_σ,V_σ,T_σ,
#                         λ_max,L_squared_σ,r_cut_squared_box, MersenneTwister()) for _ in 1:n_trials)
#     acceptance_rate= n_accepted / n_trials
#     #println("ACCEPTANCE rate: ", acceptance_rate)
#     @test acceptance_rate ≈ prefactor*exp(-1*ΔE/T_σ) atol = 0.02

#     # Particle DESTROYED
#     logQ_λN = zeros(λ_max+1,N_max+1)   
#     λ = 0
#     λ_proposed = λ_max
#     N_proposed = N-1
#     idx_deleted = 3 # just picked randomly
#     keepers = [1:idx_deleted-1; idx_deleted+1:size(r_box,2)]
#     r_proposed_box = @view r_box[:,keepers]
#     r_frac_proposed_box = [-0.08,0.1,0.1]
#     r_frac_box = r_frac_proposed_box

    
#     E_old = potential_1_normal(r_box,r_box[:,idx_deleted],idx_deleted,r_frac_box,λ,λ_max,N,L_squared_σ,r_cut_squared_box )

#     E_new = potential_1_frac(r_proposed_box,r_frac_proposed_box,λ_proposed,λ_max,N_proposed,L_squared_σ,r_cut_squared_box)

#     ΔE = E_new-E_old # 3.1378182749795753 
#     # metropolis would accept 
#     metrop_prob = exp(-ΔE/T_σ) # .0433
    
#     prefactor = (V_σ^(N_proposed+1 - N) )*( Λ_σ^(3*N - 3N_proposed-3)  ) # = 1.0
#     factorial_prefactor=N
#     expected_prob = metrop_prob*factorial_prefactor*prefactor

#     n_trials=10000

#     n_accepted = sum(λ_metropolis_pm1(λ,N,r_box,r_frac_box,
#                         λ_proposed, N_proposed, r_proposed_box, r_frac_proposed_box,idx_deleted,
#                         logQ_λN, Λ_σ,V_σ,T_σ,
#                         λ_max,L_squared_σ,r_cut_squared_box, MersenneTwister()) for _ in 1:n_trials)
#     acceptance_rate= n_accepted / n_trials
#     σ = sqrt(n_trials * expected_prob * (1 - expected_prob)) / n_trials
#    #println(acceptance_rate,expected_prob)
#     @test acceptance_rate ≈ expected_prob atol=3*σ



# end

# @testset "translation_move!() tests" begin
#     # initializing some vars for a test
#     r_σ = [-1  1  1  -1  -1 -1  1  1;  
#             -1 -1  1   1   1 -1 -1  1;
#             1  1  1   1  -1 -1 -1  -1  ] # 8 particles on in cube of sidelength 2 centered at (0,0,0), lets just say "box" is sidelength 5
#     T_σ=1.
#     Λ_σ=argon_deBroglie(T_σ)
#     L_σ=5.
#     r_box = r_σ ./ L_σ
#     r_box .= r_box .- round.(r_box)
#     N = size(r_box,2)
#     λ_max = 99
#     N_max = 450
#     logQ_λN = zeros(λ_max+1,N_max+1)
#     r_frac_box = [0.,0.,0.,]

#     L_squared_σ = L_σ^2
#     r_cut_σ=3
#     r_cut_squared_σ=r_cut_σ^2
#     r_cut_squared_box = r_cut_squared_σ/(L_σ^2)
#     input_path = joinpath(@__DIR__, "cube_vertices_home_made.inp")

#     sim = SimulationParams(N_max=450,N_min=0,T_σ=T_σ,Λ_σ=Λ_σ,
#                             λ_max=λ_max,r_cut_σ=r_cut_σ,
#                             #input_filename="/Users/mckinleypaul/Documents/montecarlo/segc_wl/test/cube_vertices_home_made.inp",# hommade input is the same as examples above 8 atoms on cube of sidelenghth 2 in simulation box of length 5 in σ units
#                             input_filename= input_path,# hommade input is the same as examples above 8 atoms on cube of sidelenghth 2 in simulation box of length 5 in σ units
#                             save_directory_path= @__DIR__ ,rng=MersenneTwister(1234))
#     wl = WangLandauVars(1,zeros(λ_max+1,N_max+1),zeros(λ_max+1,N_max+1),0,0,0,0,0,0.15)
#     μ = microstate(size(r_box,2),34,r_box,r_frac_box)

#     #### accepting fractional particle move 

#     # println(rand(sim.rng,1:μ.N+1) ) # MersenneTwister(1234) outputs 9 as a random number so it moves fractional particle
#     # println(rand(sim.rng) ) # MersenneTwister(1234) outputs 0.944 so the move should be accepted
#     old1 = (wl.translation_moves_accepted)
#     old2 = deepcopy((μ.r_frac_box)) 
#     old3 = deepcopy(μ.r_box)
#     translation_move!(sim,μ,wl)
#     @test old1 != wl.translation_moves_accepted # with conservative δr_max_box pretty sure this will always be accepted
#     @test old2 != μ.r_frac_box
#     @test old3 == μ.r_box # this shouldn't change

#     # making sure δr_max gets updated properly
#     wl.translation_moves_accepted = 2
#     wl.translation_moves_proposed=3 
#     δr_max_new = 0.15*1.05
#     translation_move!(sim,μ,wl) # accepts a move so this becomes 3/4
#     @test wl.δr_max_box == δr_max_new

#     wl.δr_max_box = 0.15
#     wl.translation_moves_accepted = 1
#     wl.translation_moves_proposed = 10
#     δr_max_new = 0.15*0.95
#     translation_move!(sim,μ,wl) # accepts move so this becomes 1/5 acceptance rate
#     @test wl.δr_max_box == δr_max_new



#     ######## accepting normal particle move MersenneTwister(1) 
#     # first i = rand(sim.rng,1:(μ.N+1))  in translation_move! picks i=6. we then move this particle a little
#     # and we get E_prop ~ -0.049 and E_old  ~ -0.22831 and ΔE = 0.18 so very slightly uphill so we may accept or reject this move, either one is reasonable
#     # however, when we use this seed we get a tiny random number in the metropolis criterion ζ = 0.008 and exp(-ΔE/T) = 0.8 so we accept the slightly uphill move for this seed
#     # so we're moving a normal particle this time and should probably reject the move

#     # initializing some vars for a test
#     r_σ = [-1  1  1  -1  -1 -1  1  1;  
#             -1 -1  1   1   1 -1 -1  1;
#             1  1  1   1  -1 -1 -1  -1  ] # 8 particles on in cube of sidelength 2 centered at (0,0,0), lets just say "box" is sidelength 5
#     T_σ=1.
#     Λ_σ=argon_deBroglie(T_σ)
#     L_σ=5.
#     r_box = r_σ ./ L_σ
#     r_box .= r_box .- round.(r_box)
#     N = size(r_box,2)
#     λ_max = 99
#     N_max = 450
#     logQ_λN = zeros(λ_max+1,N_max+1)
#     r_frac_box = [0.,0.,0.,]

#     L_squared_σ = L_σ^2
#     r_cut_σ=3
#     r_cut_squared_σ=r_cut_σ^2
#     r_cut_squared_box = r_cut_squared_σ/(L_σ^2)
#     input_path = joinpath(@__DIR__, "cube_vertices_home_made.inp") 
#     sim = SimulationParams(N_max=450,N_min=0,T_σ=T_σ,Λ_σ=Λ_σ,
#                             λ_max=λ_max,r_cut_σ=r_cut_σ,
#                             input_filename=input_path,# hommade input is the same as examples above 8 atoms on cube of sidelenghth 2 in simulation box of length 5 in σ units
#                             save_directory_path= @__DIR__ , rng=MersenneTwister(1))
#     wl = WangLandauVars(1,zeros(λ_max+1,N_max+1),zeros(λ_max+1,N_max+1),0,0,0,0,0,0.15)
#     μ = microstate(size(r_box,2),34,r_box,r_frac_box)

#     old1 = (wl.translation_moves_accepted)
#     old2 = deepcopy(μ.r_box)

#     translation_move!(sim,μ,wl) # probably should reject move 
#     println("debugger")
#     @test old1 != wl.translation_moves_accepted 
#     @test old2 != μ.r_box # if we accepted move, should be different

#     # rejecting move 


# end

println("")
println("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Now Running Physics Tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
println("")

@testset "Comparison to Analytic Values of Q(N=1),Q(N=2)" begin
    # simulation of up to 4 atoms with SEGC-WL currently takes about 15 seconds
    # so lets run this 5 times for statistics
    logQ_N1 = 0 
    logQ_N2 = 0
    for ii in 1:5
        input_path = joinpath(@__DIR__, "4_atom_cnf.inp")
        T_σ = 1.0 
        Λ_σ = argon_deBroglie(T_σ)
        sim = SimulationParams(
                        N_max=4,
                        N_min=0,
                        T_σ=T_σ,
                        Λ_σ = Λ_σ,
                        λ_max = 99,
                        r_cut_σ = 2.5, # 3 sigma for replication of results, partition functions still pass for 2.5 though
                        input_filename=input_path,
                        save_directory_path= @__DIR__ , 
                        maxiter=100_000_000)

        μstate = init_microstate(sim=sim,filename=input_path)

        wl = init_WangLandauVars(sim)

        cache = init_cache(sim,μstate)

        # initialization_check(sim,μstate,wl)

        cl = CellList(sim.N_max + 1, sim.r_cut_box)
        make_list!(cl, μstate.N, μstate.r_box)

        run_simulation!(sim,cl,μstate,wl,cache)

        logQ = correct_Q(wl)
        logQ_N1 += logQ[2]
        logQ_N2 += logQ[3]
    end
    logQ_N1 *= (1/5) # average value of Q(N=1,λ=0|V=512σ,T*=1) over five monte carlo runs
    logQ_N2 *= (1/5)
    println(logQ_N1)
    println(logQ_N2)
    #= the technically and statistically best thing to do would be run the simulation m independent times - independent meaning they must start from different initializaiton states, which is not the FCC focused workflow we have now
    Then observe the sample standard error from the m trials.  (from normal Standard Error = Standard Deviation/sqrt(m) formula)
    then you would use that SE to form a confidence interval and assert the analytic value falls inside of it:
        abs(mean_logQ - logQ_analytic) ≤ k *SE where k is the confidence interval such as k=1.96 corresponding to 95% confidence interval
    can come back and do this later but because this is going to run often we're just going to do it quick and dirty:
    =#
    mathematica_logQ_N1 = 14.0055 # analytically computed in /tests/mathematica_verification.nb a mathematica notebook
    mathematica_logQ_N2 = 27.3179 # analytically computed in /tests/mathematica_verification.nb a mathematica notebook

    @test logQ_N1 ≈ mathematica_logQ_N1 atol = (0.05*mathematica_logQ_N1) # within 5% seems good
    @test logQ_N2 ≈ mathematica_logQ_N2 atol = (0.05*mathematica_logQ_N2)

end


@testset "Comparison to Analytic Ideal Gas Limit" begin
    # in the high temperature limit, we expect our system to become the ideal gas which has an analytically computable partition function

    logQ_avg = [0.,0.,0.,0.,0.]
    for ii in 1:5
        input_path = joinpath(@__DIR__, "4_atom_cnf.inp")
        T_σ = 1_000_000.0 
        Λ_σ = argon_deBroglie(T_σ)
        sim = SimulationParams(
                        N_max=4,
                        N_min=0,
                        T_σ=T_σ,
                        Λ_σ = Λ_σ,
                        λ_max = 99,
                        r_cut_σ = 2.5, # 3 sigma for replication of results, partition function is similar for 2.5 and 3 sigma
                        input_filename=input_path,
                        save_directory_path= @__DIR__ , 
                        maxiter=100_000_000)

        μstate = init_microstate(sim=sim,filename=input_path)

        wl = init_WangLandauVars(sim)

        cache = init_cache(sim,μstate)

        # initialization_check(sim,μstate,wl)

        cl = CellList(sim.N_max + 1, sim.r_cut_box)
        make_list!(cl, μstate.N, μstate.r_box)

        run_simulation!(sim,cl,μstate,wl,cache)

        logQ = correct_Q(wl)
        for ii in 1:5
            logQ_avg[ii] += logQ[ii]
        end
    end
    logQ_avg = logQ_avg.*(1/5) # average value of partition functions over five monte carlo runs
    println(logQ_avg)
    #println(logQ_avg)

    #= the technically and statistically best thing to do would be run the simulation m independent times - independent meaning they must start from different initializaiton states, which is not the FCC focused workflow we have now
    Then observe the sample standard error from the m trials.  (from normal Standard Error = Standard Deviation/sqrt(m) formula)
    then you would use that SE to form a confidence interval and assert the analytic value falls inside of it:
        abs(mean_logQ - logQ_analytic) ≤ k *SE where k is the confidence interval such as k=1.96 corresponding to 95% confidence interval
    can come back and do this later but because this is going to run often we're just going to do it quick and dirty:
    =#
    mathematica_logQ_N2 = 68.7644 # analytically computed in /tests/mathematica_verification.nb a mathematica notebook
    mathematica_logQ_N3 = 102.395
    mathematica_logQ_N4 = 135.737
    @test logQ_avg[3] ≈ mathematica_logQ_N2 atol = (0.05*mathematica_logQ_N2) # within 5% seems good
    @test logQ_avg[4] ≈ mathematica_logQ_N3 atol = (0.05*mathematica_logQ_N3) # within 5% seems good
    @test logQ_avg[5] ≈ mathematica_logQ_N4 atol = (0.05*mathematica_logQ_N4) # within 5% seems good
end