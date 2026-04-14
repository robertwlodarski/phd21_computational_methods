# MIT Shock Results: Plotting
# Contents: 
# 1. Plotting panels 
# 2. Compute key series
# 3. Main plotting call

# 1. Plotting panels 
function fnPanelPlot(series_vec, ss_vec, labels, suptitle;
                     Tᵖˡᵒᵗ = length(series_vec[1]),
                     colors = [:blue, :red, :green, :orange, :purple, :teal, :magenta, :brown, :maroon, :navy])

    t_grid = 1:Tᵖˡᵒᵗ
    N      = length(series_vec)
    ncols  = 4
    nrows  = ceil(Int, N / ncols)

    plots  = []
    for i in 1:N
        p = plot(t_grid, series_vec[i][1:Tᵖˡᵒᵗ],
                 lw     = 2,
                 color  = colors[mod1(i, length(colors))],
                 label  = "Transition",
                 legend = :bottomright,
                 grid   = true,
                 xlabel = "",
                 ylabel = "",
                 title  = labels[i])
        hline!(p, [ss_vec[i]],
               lw        = 1.5,
               ls        = :dash,
               color     = :black,
               label     = "Steady state")
        push!(plots, p)
    end

    # Pad with empty plots if odd number
    if N % ncols != 0
        push!(plots, plot(framestyle=:none, legend=false))
    end

    plt = plot(plots...,
               layout = (nrows, ncols),
               size   = (1700, 350 * nrows),
               margin = 5Plots.mm,
               plot_title = suptitle)
    display(plt)
    return plt
end

# 2. Compute key series
function fnComputeMITPlotSeries(params, mit_endo, ss_endo,  A⃗)

    @unpack Nᶻ, Nᵃ, Nˡ, Nᵘ, Tᴹᴵᵀ, a⃗, z⃗, α, θ = params

    # A. Total credit (external finance) per period
    TotalCredit     = zeros(Tᴹᴵᵀ)
    ss_TotalCredit  = 0.0
    for it in 1:Tᴹᴵᵀ
        ωᵃ              = dropdims(sum(@views(mit_endo.g[:,:,:,:,it]), dims=(3,4)), dims=(3,4))
        TotalCredit[it] = sum(ωᵃ .* max.(@views(mit_endo.𝐤[:,:,it]) .- a⃗', 0.0) .* @views(mit_endo.𝐨[:,:,it]))
    end
    ωᵃ_ss          = dropdims(sum(ss_endo.g, dims=(3,4)), dims=(3,4))
    ss_TotalCredit = sum(ωᵃ_ss .* max.(ss_endo.𝐤 .- a⃗', 0.0) .* ss_endo.𝐨)

    # B. Average entrepreneurial productivity per period
    AvgProd     = zeros(Tᴹᴵᵀ)
    ss_AvgProd  = 0.0
    for it in 1:Tᴹᴵᵀ
        ωᵃ          = dropdims(sum(@views(mit_endo.g[:,:,:,:,it]), dims=(3,4)), dims=(3,4))
        ent_mass    = sum(ωᵃ .* @views(mit_endo.𝐨[:,:,it]))
        AvgProd[it] = ent_mass > 0 ? sum(ωᵃ .* z⃗ .* @views(mit_endo.𝐨[:,:,it])) / ent_mass : 0.0
    end
    ss_ent_mass = sum(ωᵃ_ss .* ss_endo.𝐨)
    ss_AvgProd  = ss_ent_mass > 0 ? sum(ωᵃ_ss .* z⃗ .* ss_endo.𝐨) / ss_ent_mass : 0.0

    # C. Output per period
    Output      = zeros(Tᴹᴵᵀ)
    ss_Output   = 0.0
    for it in 1:Tᴹᴵᵀ
        ωᵃ          = dropdims(sum(@views(mit_endo.g[:,:,:,:,it]), dims=(3,4)), dims=(3,4))
        Output[it]  = sum(ωᵃ .* A⃗[it] .* z⃗ .* @views(mit_endo.𝐤[:,:,it]).^α .* @views(mit_endo.𝐥[:,:,it]).^θ .* @views(mit_endo.𝐨[:,:,it]))
    end
    ss_Output = sum(ωᵃ_ss .* params.A .* z⃗ .* ss_endo.𝐤.^α .* ss_endo.𝐥.^θ .* ss_endo.𝐨)

    # D. Investment: I_t = K_{t+1} - (1-δ)K_t
    Investment      = zeros(Tᴹᴵᵀ)
    for it in 1:Tᴹᴵᵀ-1
        Investment[it] = mit_endo.Kˢ[it] - (1 - params.δ) * mit_endo.Kᵈ[it]
    end
    Investment[Tᴹᴵᵀ] = Investment[Tᴹᴵᵀ-1]
    ss_Investment     = ss_endo.Kˢ - (1 - params.δ) * ss_endo.Kᵈ

    return TotalCredit, ss_TotalCredit, AvgProd, ss_AvgProd, Output, ss_Output, Investment, ss_Investment
end

# 3. Main plotting call
function fnPlotMITResults(params, mit_endo, ss_endo,  A⃗; Tᵖˡᵒᵗ = params.Tᴹᴵᵀ)

    # A. Compute derived series
    TotalCredit, ss_TotalCredit, AvgProd, ss_AvgProd, Output, ss_Output, Investment, ss_Investment = fnComputeMITPlotSeries(params, mit_endo, ss_endo, A⃗)

    # B. Collect all series and steady-state values
    series_vec = [
        Output,
        Investment,
        mit_endo.U,
        mit_endo.Kᵈ,
        TotalCredit,
        mit_endo.τₜ,
        mit_endo.E,
        mit_endo.rₜ,
        mit_endo.wₜ,
        AvgProd
    ]

    ss_vec = [
        ss_Output,
        ss_Investment,
        ss_endo.U,
        ss_endo.Kᵈ,
        ss_TotalCredit,
        ss_endo.τₜ,
        ss_endo.E,
        ss_endo.rₜ,
        ss_endo.wₜ,
        ss_AvgProd
    ]

    labels = [
        "Yₜ",
        "Investment",
        "Uₜ",
        "Kₜ",
        "Credit",
        "τ",
        "𝔼[o(z,a)=1]",
        "rₜ",
        "wₜ",
        "𝔼[z]"
    ]

    # C. Plot
    fnPanelPlot(series_vec, ss_vec, labels, ""; Tᵖˡᵒᵗ = Tᵖˡᵒᵗ)
end