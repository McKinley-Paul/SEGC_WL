# This file contains functions general to doing thermodynamics or stat mech including some thermodynamic post processing of the SEGC-WL found configurational integrals 
# also contains some functions to do analysis
#  ✅  == checked in /test

function correct_logQ(wl::WangLandauVars)::Array{Float64}   
    #❌ This function does not take into account wl.N_min yet, assumes N_min==0 ❌
    #= after doing WL, we solve only for the simplified expanded canonical partition functions  Q(λ,N|V,T) up to a multiplicative constant
     this function solves for this constant and returns the true Q's by recognizing that we know 
        Q(N=0,l=0,|T,V) == 1 
     as such we solve for the arbitrary constant, C, by doing:
         Q_true(N=0,l=0) = 1 = C * Q_wl(N=0,l=0)
         1/Q_wl(N=0,l=0) = C 
     then we can solve for the true partition functions by doing:
         Qtrue(N,l=0) = C * Q_wl(N,l=0)
     but of course all this is done in log space to get most precision possible.

    returns the corrected,true partition functions as a Vector. logQ[1] = logQ(N=0|V,T) and so on
    =#
    logC = -1*wl.logQ_λN[1,1] # natural log of multiplicative constant 
    logQtrue = wl.logQ_λN[1,:] .+ logC
    return(logQtrue)
end #correct_Q

