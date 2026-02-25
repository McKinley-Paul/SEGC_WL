function highest_N_threshold_entry(A; threshold = 1000)
    # Returns the last N,λ for which there is a non-1 entry. measure of high far right the mc sampler sampled
    #sweeps through the array backward (right to left column wise and down to up row wise) and picks out the first 
    # entry greater than 1 and returns its value and its index. basically just go backwards if you flatten the array 
    # to the vector
    v = vec(A)

    for k in length(v):-1:1
        if v[k] > threshold
            I = CartesianIndices(A)[k]
            return (linear_index = k,
                    row = I[1],
                    col = I[2],
                    value = v[k])
        end
    end

    return nothing
end

function analyze_experiment(wl,bashout_path)
    println("wl.logf: ", wl.logf)
    wl_data = wl.H_λN
    open(bashout_path, "r") do f
        last_line = last(eachline(f))
        println(last_line)
    end
    println("λ move acceptance rate: ", round( 100*(wl.λ_moves_accepted/ (wl.λ_moves_proposed+wl.λ_moves_accepted) ) ), " %")
    println("Translation move acceptance rate: ", round( 100*(wl.translation_moves_accepted/ (wl.translation_moves_proposed+wl.translation_moves_accepted) ) ), " %")

    #wl_data = full histogram wl.H_λN
    println("Minimum of all rows and columns from 401 to 500 (expect min=1):",minimum(wl_data[:,401:500])) # hasnt even visited all states
    # probably it rejects the state with 500 particles and lambda = 1 (in wl_control.H_λN[2,501] ) so don't expect those states to be visited
    # other states from 500 lambda = 0 down probably should be visited though hence the above
    println("Minimum of all rows and columns from 401 to 501 (expect min=0 but higher better):",minimum(wl_data[:,401:501]))
    println("Maximum: ", maximum(wl_data)) 

    μ = mean(wl_data[:,401:501])
    σ = std(wl_data[:,401:501])
    println("Mean of histogram is: ", round(μ))
    println("Standard Deviation of histogram is: ", round(σ))
    println("Coefficient of variation, σ/|μ| = ", σ/abs(μ))
    println("σ/|μ| = 0 is perfectly flat, σ/|μ|=0.05 has some structure and σ/|μ|=0.1 is definitely not flat")
    _,λ_high,Nhigh,first_non_1_value = highest_N_threshold_entry(wl_data[:,401:501])
    println("The last state sampled at least a thousand times was (N,λ) = (",Nhigh+400,",",λ_high-1,") with histogram value: ", first_non_1_value)

    col_avg = mean(wl_data, dims=1)[:]     
    col_idx = 1:size(wl_data, 2)

    p1 = plot(
        col_idx,
        col_avg,
        xlabel = "Column number",
        ylabel = "Average value",
        xlims = (399,502),
        title = "Column-wise average (linear scale)",
        legend = false, 
        xticks = 399:10:502,
        size = (1200, 300),
    )
    vline!(p1, [Nhigh+400], color=:red, linestyle=:dash, label="Threshold",linewidth=3)
    display(p1)

    p2 = plot(
        col_idx,
        col_avg,
        xlabel = "Column number",
        ylabel = "Average value",
        xlims = (399,502),
        ylims = (0,1000),
        title = "Column-wise average (linear scale)",
        legend = false, 
        xticks = 399:10:502,
        size = (1200, 300),    )
    vline!(p2, [Nhigh+400], color=:red, linestyle=:dash, label="Threshold",linewidth=3)

    display(plot(p2))

    display(heatmap(wl_data,
        aspect_ratio = 1,
        colorbar_title = "Value",
        title = "2D field",
        xlims = (400,501)
        ))

    display(heatmap(wl_data .- μ,
        aspect_ratio = 1,
        colorbar_title = "A - mean",
        title = "Deviation from mean",
        color = :balance,
        xlims = (400,501)
        ))
end
