module segc_wl
export run_simulation



function run_simulation(sim::SimulationParams, μstate::microstate)

    check_inputs(sim,μstate) # tell user the min distance in config and similar things
    print_simulation_params(sim)
    print_microstate(μstate,true)


end

end