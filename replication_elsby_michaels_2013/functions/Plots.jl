# 1. Beveridge curve

function fnPlotBeveridgeCurve(simu, params)
    # A. Unpacking business 
    @unpack L = params

    # B. Preparing data
    U⃗ = (L .- simu.N⃗) ./ L
    V⃗ = (simu.M⃗ ./ simu.q⃗) ./ L
    Tq = length(U⃗) ÷ 13
    U⃗q = [mean(U⃗[(t-1)*13+1 : t*13]) for t in 1:Tq]
    V⃗q = [mean(V⃗[(t-1)*13+1 : t*13]) for t in 1:Tq]
    Ũ = U⃗q .- mean(U⃗q)
    Ṽ = V⃗q .- mean(V⃗q)

    # C. Plot 
    plt = scatter(Ũ, Ṽ, xlabel  ="Unemployment rate (U/L)", 
                        ylabel  ="Vacancy rate (V/L)",
                        color   = :maroon,
                        legend  =false, 
                        ms      =3, 
                        mc      =:gray60,
                        markerstrokewidth=0.5, 
                        grid=true, size=(600,500))
    savefig(plt, "plots/beveridge_curve.pdf")
    return plt
end

# 2. Impulse responses 
function fnPlotStateDepIRF(simu, lee, params; T_irf=174, T_shock=12)
    @unpack L, Nₓ, x̲, ξ, p̄ₓ, α, λ = params
    T = length(lee.N⃗)

    # A. Select initial conditions: boom (90th pctile) and recession (10th pctile)
    N_p90 = quantile(lee.N⃗, 0.90)
    N_p10 = quantile(lee.N⃗, 0.10)
    t_boom = argmin(abs.(lee.N⃗[1:T-T_irf] .- N_p90))
    t_rec  = argmin(abs.(lee.N⃗[1:T-T_irf] .- N_p10))

    # B. Forward simulate with and without shock
    function forward_irf(t₀)
        # Baseline: just read off converged paths
        U_base = (L .- lee.N⃗[t₀:t₀+T_irf-1]) ./ L
        f_base = simu.f⃗[t₀:t₀+T_irf-1]
        S_base = simu.S⃗[t₀:t₀+T_irf-1]
        V_base = (simu.M⃗[t₀:t₀+T_irf-1] ./ simu.q⃗[t₀:t₀+T_irf-1])

        # Shocked: perturb p by -1% for T_shock weeks, then forward simulate
        N_sh = zeros(T_irf)
        f_sh = zeros(T_irf)
        S_sh = zeros(T_irf)
        M_sh = zeros(T_irf)
        q_sh = zeros(T_irf)
        N_sh[1] = lee.N⃗[t₀]

        for j in 1:T_irf
            t = t₀ + j - 1
            t = mod1(t, T)  # wrap around
            pₜ = simu.p⃗̂[t] * (j <= T_shock ? 0.99 : 1.0)

            # Use converged policy functions at period t
            Nⱼ, _, Sⱼ, Mⱼ, _ = fAggregationAggregateUncertainty(
                params, lee.R⃗ᵥ[t], lee.∂R⃗ᵥ[t], lee.R⃗[t], lee.∂R⃗[t], lee.n⃗[t], pₜ)
            S_sh[j] = Sⱼ
            M_sh[j] = Mⱼ
            f_sh[j] = Mⱼ / (L - N_sh[j] + 1e-8)
            q_sh[j] = fUpdatedJobFindingRateInverse(f_sh[j], params)
            if j < T_irf
                N_sh[j+1] = N_sh[j] + Mⱼ - Sⱼ
            end
        end

        U_sh = (L .- N_sh) ./ L
        V_sh = M_sh ./ q_sh
        θ_base = V_base ./ ((L .- lee.N⃗[t₀:t₀+T_irf-1]))
        θ_sh   = V_sh ./ (L .- N_sh)
        s_base = S_base ./ L
        s_sh   = S_sh ./ L

        return U_base, U_sh, θ_base, θ_sh, f_base, f_sh, s_base, s_sh
    end

    # C. Compute IRFs
    Ub_b, Us_b, θb_b, θs_b, fb_b, fs_b, sb_b, ss_b = forward_irf(t_boom)
    Ub_r, Us_r, θb_r, θs_r, fb_r, fs_r, sb_r, ss_r = forward_irf(t_rec)

    # D. Weeks to months
    T_mo = T_irf ÷ 4  
    avg_mo(x) = [mean(x[(i-1)*4+1:min(i*4, length(x))]) for i in 1:T_mo]

    # E. Plot
    months = 0:T_mo-1
    lw_b, lw_r = 2, 2
    lbl_b, lbl_r = "Boom (P90)", "Recession (P10)"

    p1 = plot(months, avg_mo(Us_b), lw=lw_b, label=lbl_b, color=:navy,
              title="Panel A. Unemployment rate", ylabel="U/L")
    plot!(p1, months, avg_mo(Us_r), lw=lw_r, label=lbl_r, color=:maroon, ls=:dash)

    p2 = plot(months, avg_mo(θs_b), lw=lw_b, label=false, color=:navy,
              title="Panel B. Vacancy-unemployment ratio", ylabel="θ")
    plot!(p2, months, avg_mo(θs_r), lw=lw_r, label=false, color=:maroon, ls=:dash)

    p3 = plot(months, avg_mo(fs_b), lw=lw_b, label=false, color=:navy,
              title="Panel C. Job-finding rate", ylabel="f")
    plot!(p3, months, avg_mo(fs_r), lw=lw_r, label=false, color=:maroon, ls=:dash)

    p4 = plot(months, avg_mo(ss_b), lw=lw_b, label=false, color=:navy,
              title="Panel D. Inflow rate", ylabel="S/L", xlabel="Months since shock")
    plot!(p4, months, avg_mo(ss_r), lw=lw_r, label=false, color=:maroon, ls=:dash)

    plt = plot(p1, p2, p3, p4, layout=(2,2), size=(900,700), margin=5Plots.mm)
    savefig(plt, "plots/state_dependent_irf.pdf")
    return plt
end

# 3. Simulated paths 
function fnPlotSimulatedPaths(simu, lee, params; path = "plots/simulated_paths.pdf", λ_hp = 100_000)
    T = length(simu.p⃗̂)
    t⃗ = 1:T
    trend(x) = hp_filter(x, λ_hp)[2]

    # Identify recession bands: runs of p < 1 lasting > 24 weeks
    low = simu.p⃗̂ .< 1.0
    bands = Tuple{Int,Int}[]
    i = 1
    while i <= T
        if low[i]
            j = i
            while j <= T && low[j]; j += 1; end
            if j - i > 24
                push!(bands, (i, j-1))
            end
            i = j
        else
            i += 1
        end
    end

    shade!(p) = for (a,b) in bands; vspan!(p, [a,b], color=:gray90, alpha=0.5, label=false, lw=0); end

    p1 = plot(t⃗, simu.p⃗̂, title="Aggregate productivity", ylabel="pₜ",
              lw=0.5, color=:gray70, legend=false, grid=true, yformatter = y -> @sprintf("%.3g", y))
    plot!(p1, t⃗, trend(simu.p⃗̂), lw=2, color=:black)
    shade!(p1)

    p2 = plot(t⃗, lee.N⃗, title="Aggregate employment", ylabel="Nₜ",
              lw=0.5, color=:lightblue, legend=false, grid=true, yformatter = y -> @sprintf("%.3g", y))
    plot!(p2, t⃗, trend(lee.N⃗), lw=2, color=:navy)
    shade!(p2)

    p3 = plot(t⃗, lee.q⃗, title="Job-filling rate", ylabel="qₜ", xlabel="t (weeks)",
              lw=0.5, color=:lightsalmon, legend=false, grid=true, yformatter = y -> @sprintf("%.3g", y))
    plot!(p3, t⃗, trend(lee.q⃗), lw=2, color=:maroon)
    shade!(p3)

    plt = plot(p1, p2, p3, layout=(3,1), size=(900,700), margin=5Plots.mm)
    savefig(plt, path)
    return plt
end