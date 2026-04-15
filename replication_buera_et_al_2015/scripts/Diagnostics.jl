# Diagnosing what's not okay with the model 
# Content:
# 1. Number of workers close to the consumption floor
# 2. Checking the productivity process

# 1. Number of workers close to the consumption floor
function fnConsumptionFloor(endo,params)
    @unpack a⃗, z⃗, c̲ = params
    ωᵃ              = dropdims(sum(endo.g, dims=(3,4)), dims=(3,4))
    threshold       = 1   # or pick whatever "close" means
    mass_near_floor = sum((ωᵃ[ij] for ij in eachindex(ωᵃ) if endo.𝐜[ij] < threshold), init=0.0)
    println("Mass near c̲: $(round(mass_near_floor, digits=4)) (threshold = $threshold)")
end 
fnConsumptionFloor(Endo,UsedParameters)

# 2. Checking the productivity process 
function fnPlotProductivityDistribution(endo, params)
    @unpack z⃗, μ⃗, η = params

    # A. Model distribution (marginal over z)
    ωᵃ = dropdims(sum(endo.g, dims=(3,4)), dims=(3,4))
    model_pdf = vec(sum(ωᵃ, dims=2))  # sum over assets

    # B. Theoretical Pareto PDF: f(z) = η * z^(-(η+1))
    pareto_pdf = η .* z⃗ .^ (-(η + 1))
    pareto_pdf ./= sum(pareto_pdf)  # normalize to compare

    p = plot(z⃗, μ⃗, lw=2, color=:maroon, label="μ⃗ (grid masses)")
    plot!(z⃗, pareto_pdf, lw=2, ls=:dash, color=:navy, label="Pareto PDF (normalized)")
    display(p)
    return p
end
fnPlotProductivityDistribution(Endo, UsedParameters)

# 3. Diagnostic statistics 
function fnDiagnostics(params, endo)
    @unpack a⃗, z⃗, δ = params

    # Marginal (z,a) distribution
    ωᵃ = dropdims(sum(endo.g, dims=(3,4)), dims=(3,4))

    # Aggregate quantities
    t_e         = sum(ωᵃ .* endo.𝐨)
    t_w         = endo.W
    avg_firmsize = t_e > 0 ? endo.Lᵈ / t_e : 0.0

    # Worker-to-entrepreneur transition: share of workers (iu=2 or employed, 𝐨=false) who switch
    g_workers   = dropdims(sum(endo.g[:,:,:,1], dims=3), dims=3) .* (.!endo.𝐨) .+
                  dropdims(sum(endo.g[:,:,:,2], dims=3), dims=3)
    new_ent     = sum(ωᵃ .* endo.𝐨 .* (.!endo.𝐨))   # placeholder — see note below
    w2e_rate    = t_e > 0 ? endo.S / max(sum(g_workers), 1e-8) : 0.0

    # Investment rate: I/Y = δ*K/Y
    Y           = sum(ωᵃ .* endo.𝐨 .* (endo.Π .+ (endo.rₜ + δ) .* endo.𝐤 .+ endo.wₜ .* endo.𝐥))
    inv_rate    = Y > 0 ? δ * endo.Kᵈ / Y : 0.0

    # Credit-to-output
    ext_fin     = sum(ωᵃ .* max.(endo.𝐤 .- a⃗', 0.0) .* endo.𝐨)
    credit_Y    = Y > 0 ? ext_fin / Y : 0.0

    # Top assets mass 
    top_asset_mass = sum(@view(ωᵃ[:, end]))

    println("\n", repeat("─", 60))
    println("Diagnostics")
    println(repeat("─", 60))
    @printf("%-40s %10.4f\n", "Average firm size (workers)",     avg_firmsize)
    @printf("%-40s %10.4f\n", "Worker→Entrepreneur rate",        w2e_rate)
    @printf("%-40s %10.4f\n", "Unemployment rate",               endo.U)
    @printf("%-40s %10.4f\n", "Investment rate (δK/Y)",          inv_rate)
    @printf("%-40s %10.4f\n", "Credit-to-output ratio",          credit_Y)
    @printf("%-40s %10.4f\n", "Entrepreneur share",              t_e)
    @printf("%-40s %10.4f\n", "Employed worker share",           t_w)
    @printf("%-40s %10.4f\n", "Wage",                            endo.wₜ)
    @printf("%-40s %10.4f\n", "Interest rate",                   endo.rₜ)
    @printf("%-40s %10.4f\n", "Tax",                             endo.τₜ)
    @printf("%-40s %10.4f\n", "K demand",                        endo.Kᵈ)
    @printf("%-40s %10.4f\n", "K supply",                        endo.Kˢ)
    @printf("%-40s %10.4f\n", "L demand",                        endo.Lᵈ)
    @printf("%-40s %10.4f\n", "L supply",                        endo.Lˢ)
    @printf("%-40s %10.4f\n", "Top assets mass",                 top_asset_mass)
    println(repeat("─", 60))
end
 fnDiagnostics(UsedParameters,Endo)

 # 4. Distribution Diagnostics
 function fnCheckDistribution(params, endo)
    @unpack z⃗, a⃗, l⃗ = params

    total_mass   = sum(endo.g)
    ωᵃ           = dropdims(sum(endo.g, dims=(3,4)), dims=(3,4))
    marg_z       = vec(sum(ωᵃ, dims=2))
    marg_a       = vec(sum(ωᵃ, dims=1))
    mass_neg     = sum(x -> x < 0 ? abs(x) : 0.0, endo.g)
    mass_emp     = sum(endo.g[:,:,:,1])
    mass_unemp   = sum(endo.g[:,:,:,2])

    println("\n", repeat("─", 50))
    println("Distribution checks")
    println(repeat("─", 50))
    @printf("%-35s %10.6f\n", "Total mass (should be 1.0)",  total_mass)
    @printf("%-35s %10.6f\n", "Employed mass (iu=1)",        mass_emp)
    @printf("%-35s %10.6f\n", "Unemployed mass (iu=2)",      mass_unemp)
    @printf("%-35s %10.6f\n", "Negative mass (should be 0)", mass_neg)
    @printf("%-35s %10.6f\n", "Marginal z (sum)",            sum(marg_z))
    @printf("%-35s %10.6f\n", "Marginal a (sum)",            sum(marg_a))
    @printf("%-35s %10.6f\n", "Min g value",                 minimum(endo.g))
    @printf("%-35s %10.6f\n", "Max g value",                 maximum(endo.g))
    println(repeat("─", 50))
end
fnCheckDistribution(UsedParameters,Endo)
ωa = dropdims(sum(Endo.g, dims=(3,4)), dims=(3,4))
frac_constrained = sum(ωa .* Endo.𝕀ᶜ .* Endo.𝐨) / max(sum(ωa .* Endo.𝐨), 1e-8)
@printf("Fraction of entrepreneurs constrained: %.4f\n", frac_constrained)