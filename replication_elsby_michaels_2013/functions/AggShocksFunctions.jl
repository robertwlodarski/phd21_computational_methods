# Important note on notation with aggregate uncertainty 
# 𝒩⃗: caligraphic letters:           State space
# N⃗: standard letters:              Simulation  
# Π̃: greek letters with a tilde:    Value functions

# 1. Aggregation with aggregate uncertainty 
function fAggregationAggregateUncertainty(params,R⃗ᵥ,∂R⃗ᵥ,R⃗,∂R⃗,n⃗,p)

    # A. Unpacking business 
    @unpack x̲, x⃗, ξ, p̄ₓ, α, λ = params 

    # B. Compute CDFs, PDFs, and expectation 
    𝐆R⃗ᵥ     = (1 .- (x̲ ./ R⃗ᵥ).^ξ) ./ p̄ₓ
    𝐠R⃗ᵥ     = ((1 / p̄ₓ) * ξ * x̲^ξ) ./ ((R⃗ᵥ).^(ξ+1)) 
    𝐆R⃗      = (1 .- (x̲ ./ R⃗).^ξ) ./ p̄ₓ
    𝐠R⃗      = ((1 / p̄ₓ) * ξ * x̲^ξ) ./ ((R⃗).^(ξ+1)) 
    𝐇n⃗      = 𝐆R⃗ ./ (1 .- 𝐆R⃗ᵥ .+ 𝐆R⃗)
    𝐡n⃗      = ((1 .- 𝐆R⃗ᵥ) .* 𝐠R⃗ .* ∂R⃗ + 𝐆R⃗ .* 𝐠R⃗ᵥ .* ∂R⃗ᵥ) ./ ((1 .- 𝐆R⃗ᵥ .+ 𝐆R⃗).^2)
    𝔼x      = (1 / p̄ₓ) * x̲^ξ * (ξ /(ξ - 1)) * (R⃗.^(-ξ+1)-R⃗ᵥ.^(-ξ+1)) ./ max.(𝐆R⃗ᵥ .- 𝐆R⃗,1e-8)

    # C. Compute aggregate values 
    N       = fSimpsonRule(n⃗ .*𝐡n⃗, n⃗)                   # Employed
    Y       = fSimpsonRule(𝔼x .* n⃗.^α .*𝐡n⃗,n⃗)           # Production
    S       = λ * fSimpsonRule((1 .- 𝐇n⃗).*𝐆R⃗,n⃗)         # Separations 
    M       = λ * fSimpsonRule(𝐇n⃗.*(1 .- 𝐆R⃗ᵥ),n⃗)        # Matches 
    A       = fSimpsonRule(p .* 𝔼x .* n⃗.^(α-1).*𝐡n⃗,n⃗)   # Total marginal product of labour  
    return N, Y, S, M, A 
end 

# 2. Robust spline for aggregate uncertainty 
function fRobustSplineAggregateUncertainty(x_in, y_in, eval_grid)
    # A. Sort by the x-coordinate 
    p           = sortperm(x_in)
    x_s         = x_in[p]
    y_s         = y_in[p]

    # B. Deduplicate
    mask        = [true; diff(x_s) .> 1e-8]
    x_clean     = x_s[mask]
    y_clean     = y_s[mask]

    # C. PCHIP interpolation
    if length(x_clean) >= 2
        itp = PCHIPInterpolation(y_clean, x_clean; extrapolation = ExtrapolationType.Linear)
        vals    = itp.(eval_grid)
        derivs  = DataInterpolations.derivative.(Ref(itp), eval_grid)
        return vals, derivs
    else
        val     = isempty(y_clean) ? 0.0 : y_clean[1]
        return fill(val, length(eval_grid)), zeros(length(eval_grid))
    end
end

# 3. Solve using the repeated transition method 
function fnSolveAggregateRTM!(params, ss_endo, simu, lee; warm = true)

    # A. Unpacking business 
    @unpack P, p⃗, Nₚ, L, δᴿᵀᴹ, δᵍ, x⃗, α, β, c, Nₓ, ωᵍ, ωᴿᵀᴹ₁, ωᴿᵀᴹ₂, ωᴿᵀᴹ₃, x̅, x̲,q̲, q̅ = params
    @unpack p⃗̂, p⃗̂ᵢ, N⃗, q⃗, f⃗, S⃗, M⃗, Y⃗     = simu
    T                                   = length(p⃗̂)

    # B. Initialise predicted paths
    if isfile("results/rtm_warm_start.jld2") && warm
        println("I don't have eternity - warm start: Loading q, N, Π̃, n⃗")
        jldopen("results/rtm_warm_start.jld2", "r") do file
            lee.q⃗ .= file["warm_q"]
            lee.N⃗ .= file["warm_N"]
            for t in 1:T
                lee.n⃗[t]     = file["warm_n"][t]
                Nₙₜ          = length(lee.n⃗[t])
                lee.Π̃[t]     = file["warm_Pi"][t]
                lee.Π̃ⁱᵐᵖ[t]  = zeros(Nₓ, Nₙₜ)
                lee.R⃗[t]     = zeros(Nₙₜ)
                lee.R⃗ᵥ[t]    = zeros(Nₙₜ)
                lee.∂R⃗[t]    = zeros(Nₙₜ)
                lee.∂R⃗ᵥ[t]   = zeros(Nₙₜ)
            end
        end
        f⃗ .= fUpdatedJobFindingRate.(lee.q⃗, Ref(params))
    else
        q̂ˢˢ = zeros(Nₚ)
        N̂ˢˢ = zeros(Nₚ)
        for j in 1:Nₚ
            êⱼ      = fnSetUpEndo(params)
            q̂ˢˢ[j]  = find_zero(q -> fEqResidual(q, p⃗[j], params, êⱼ), (q̲, q̅), Bisection())
            fⱼ      = fUpdatedJobFindingRate(q̂ˢˢ[j], params)
            N̂ˢˢ[j], _, _, _, _ = fAggregation(params, êⱼ, p⃗[j], fⱼ, q̂ˢˢ[j])
        end
        lee.q⃗ .= q̂ˢˢ[p⃗̂ᵢ]
        lee.N⃗ .= N̂ˢˢ[p⃗̂ᵢ] .* (1 .+ 1e-3 .* randn(T))
        f⃗ .= fUpdatedJobFindingRate.(lee.q⃗, Ref(params))
        ℑˢˢ = [Spline1D(ss_endo.n⃗ᵛᶠⁱ, ss_endo.Π[i, :]; k=3, bc="extrapolate") for i in 1:Nₓ]
        for t in 1:T
            lee.n⃗[t]    = fn⃗(params, p⃗̂[t], f⃗[t], lee.q⃗[t])
            Nₙₜ         = length(lee.n⃗[t])
            lee.Π̃[t]    = zeros(Nₓ, Nₙₜ)
            lee.Π̃ⁱᵐᵖ[t] = zeros(Nₓ, Nₙₜ)
            lee.R⃗[t]    = zeros(Nₙₜ)
            lee.R⃗ᵥ[t]   = zeros(Nₙₜ)
            lee.∂R⃗[t]   = zeros(Nₙₜ)
            lee.∂R⃗ᵥ[t]  = zeros(Nₙₜ)
            for i in 1:Nₓ
                lee.Π̃[t][i, :] .= ℑˢˢ[i](lee.n⃗[t])
            end
        end
    end
    println("Done with setting the starting point")
    
    # C. RTM outer loop
    εᴿᵀᴹ    = Inf
    nᴿᵀᴹ    = 1
    vᶠ      = zeros(Nₓ)
    nᶠ      = zeros(Nₓ)
    vʰ      = zeros(Nₓ)
    nʰ      = zeros(Nₓ)
    ϵʰⁱˢᵗ   = Float64[]
    while εᴿᵀᴹ > δᴿᵀᴹ

        # D. Backward solution (t = T → 1)
        for t in T:-1:1
            Nₙₜ = length(lee.n⃗[t])
            Πᶠ  = zeros(Nₓ, Nₙₜ)
            Πʰ  = zeros(Nₓ, Nₙₜ)

            # E. Retrieve current period's exogenous state and nth-iteration conjectured price
            pₜ  = p⃗̂[t]
            pᵢₜ = p⃗̂ᵢ[t]
            qₜ  = lee.q⃗[t]
            
            # F. Construct 𝔼ₜ[Π̃ₜ₊₁] via RTM matching
            Π̃ₜ₊₁        = Vector{Matrix{Float64}}(undef, Nₚ)
            # Realised state handled directly (remember about grid resizing)
            τʳ          = t == T ? 1 : t+1
            Π̃ᵣ          = zeros(Nₓ, length(lee.n⃗[t]))
            for i in 1:Nₓ
                Π̃ᵣ[i,:] .= Spline1D(lee.n⃗[τʳ], lee.Π̃[τʳ][i,:]; k=3, bc="extrapolate")(lee.n⃗[t])
            end
            Π̃ₜ₊₁[pᵢₜ]   = Π̃ᵣ
            # For the realised p′: use Π̃ʳ[t+1] directly from this backward pass
            # For each counterfactual p′ ≠ pₜ₊₁: find τ s.t. N⃗ᵖ[τ] ≈ N⃗ᵖ[t+1] and p⃗̂ᵢ[τ] == p′
            for pⱼ in filter(!=(pᵢₜ), 1:Nₚ)

                # i. Find matching period τ via interpolation on N
                τ⃗ᶜ              = findall(p⃗̂ᵢ .== pⱼ)
                N⃗ᶜ              = lee.N⃗[τ⃗ᶜ]
                Nᵉᶠᶠ            = t == T ? lee.N⃗[1] : lee.N⃗[t+1]
                idx             = sortperm(abs.(N⃗ᶜ .- Nᵉᶠᶠ))
                τ̲               = τ⃗ᶜ[idx[1]]
                τ̅               = τ⃗ᶜ[idx[2]]
                ωτ              = (lee.N⃗[τ̅] - Nᵉᶠᶠ) / (lee.N⃗[τ̅] - lee.N⃗[τ̲])
                ωτ              = clamp(ωτ, 0.0, 1.0)

                # ii. Borrow and interpolate Π̃ onto lee.n⃗[t]
                Π̃τ̲  = zeros(Nₓ, length(lee.n⃗[t]))
                Π̃τ̅  = zeros(Nₓ, length(lee.n⃗[t]))
                for i in 1:Nₓ
                    Π̃τ̲[i,:] .= Spline1D(lee.n⃗[τ̲], lee.Π̃[τ̲][i,:]; k=3, bc="extrapolate")(lee.n⃗[t])
                    Π̃τ̅[i,:] .= Spline1D(lee.n⃗[τ̅], lee.Π̃[τ̅][i,:]; k=3, bc="extrapolate")(lee.n⃗[t])
                end
                Π̃ₜ₊₁[pⱼ] = ωτ .* Π̃τ̲ .+ (1 - ωτ) .* Π̃τ̅

            end

            # G. Assemble expected value: 𝔼ₜ[Π̃ₜ₊₁] = Σⱼ P[pᵢ,j] × Π̃ₜ₊₁(·; pⱼ)
            𝔼Π̃ = sum(P[pᵢₜ, pⱼ] .* Π̃ₜ₊₁[pⱼ] for pⱼ in 1:Nₚ)
            
            # H. Single Bellman update and policy extraction
            n⃗ₜ          = lee.n⃗[t]
            W           = fW(params, pₜ, f⃗[t], qₜ, n⃗ₜ)
            Πᶠˡᵒʷ       = pₜ .* x⃗ .* (n⃗ₜ' .^ α) .- W .* n⃗ₜ'
            Πᶜ          = Πᶠˡᵒʷ .+ β .* 𝔼Π̃
            for i in 1:Nₓ
                # Fire 
                valᶠ, idᶠ   = findmax(view(Πᶜ, i, :))
                vᶠ[i]       = valᶠ
                nᶠ[i]       = n⃗ₜ[idᶠ]
                # Hire 
                valʰ, idʰ   = findmax(view(Πᶜ, i, :) .- (c / qₜ) .* n⃗ₜ)
                vʰ[i]       = valʰ
                nʰ[i]       = n⃗ₜ[idʰ]
            end
            Πᶠ              .= vᶠ .* (nᶠ .< n⃗ₜ') .- 1e8 .* (nᶠ .>= n⃗ₜ')
            Πʰ              .= (vʰ .+ c / qₜ .* n⃗ₜ') .* (nʰ .> n⃗ₜ') .- 1e8 .* (nʰ .<= n⃗ₜ')
            lee.Π̃ⁱᵐᵖ[t]     = max.(Πᶠ, max.(Πʰ, Πᶜ))

            # I. Extract and store policy functions R⃗ₜ, R⃗ᵥₜ, n⃗ₜ into lee.R⃗ᵀ[t] etc.
            𝓃₁              = findlast(nᶠ .< nʰ[Nₓ])
            𝓃₂              = findfirst(nʰ .> nᶠ[1])
            isnothing(𝓃₁) && (𝓃₁ = Nₓ)
            isnothing(𝓃₂) && (𝓃₂ = 1)

            # I1. Firing threshold 
            𝕟ᴿ                  = nᶠ[1:𝓃₁]
            𝕩ᴿ                  = x⃗[1:𝓃₁]
            if length(𝕟ᴿ) < 2
                lee.R⃗[t]  .= x̲
                lee.∂R⃗[t] .= 0.0
            else
                lee.R⃗[t], lee.∂R⃗[t] = fRobustSplineAggregateUncertainty(𝕟ᴿ, 𝕩ᴿ, lee.n⃗[t])
            end
            lee.R⃗[t]            = clamp.(lee.R⃗[t], x̲, x̅)

            # I2. Hiring threshold 
            𝕟ᴿⱽ             = nʰ[𝓃₂:end]
            𝕩ᴿⱽ             = x⃗[𝓃₂:end]
            # R⃗ᵥᵗᵐᵖ, ∂R⃗ᵥᵗᵐᵖ   = fRobustSplineAggregateUncertainty(𝕟ᴿⱽ, 𝕩ᴿⱽ, lee.n⃗[t])
            # lee.R⃗ᵥ[t]       = clamp.(min.(R⃗ᵥᵗᵐᵖ, x̅), x̲, x̅)
            # lee.∂R⃗ᵥ[t]      = ∂R⃗ᵥᵗᵐᵖ
            if length(𝕟ᴿⱽ) < 2
                lee.R⃗ᵥ[t]       .= x̅
                lee.∂R⃗ᵥ[t]      .= 0.0
            else
                R⃗ᵥᵗᵐᵖ, ∂R⃗ᵥᵗᵐᵖ   = fRobustSplineAggregateUncertainty(𝕟ᴿⱽ, 𝕩ᴿⱽ, lee.n⃗[t])
                lee.R⃗ᵥ[t]       = clamp.(min.(R⃗ᵥᵗᵐᵖ, x̅), x̲, x̅)
                lee.∂R⃗ᵥ[t]      = ∂R⃗ᵥᵗᵐᵖ
            end 

            # I3. Ensure no impossible results 
            lee.R⃗ᵥ[t]       = max.(lee.R⃗ᵥ[t], lee.R⃗[t] .* (1 + 1e-3))
        end 

        #%% J. Forward simulation (t = 1 → T)
        Nₜᵖʳᵉᵛ = lee.N⃗[T]
        for t in 1:T

            # J1. Retrieve pₜ and Nₜ₋₁
            pₜ      = p⃗̂[t]
            Nₜ₋₁    = Nₜᵖʳᵉᵛ

            # J2. Solve equilibrium q loop (reuse fEqResidual logic)
            # qₜ          = lee.q⃗[t]     
            # Nₜ          = lee.N⃗[t]
            # εᵍ          = Inf
            # nᵍ          = 0
            # while εᵍ > δᵍ && nᵍ < 500

            #     # !!!. DIAGNOSTIC WINDOW (delete once life is good again)
            #     N_t, Y_t, S_t, M_t, A_t = fAggregationAggregateUncertainty(params, lee.R⃗ᵥ[t], lee.∂R⃗ᵥ[t], lee.R⃗[t], lee.∂R⃗[t], lee.n⃗[t], pₜ)
            #     if isnan(N_t) || isnan(S_t) || isnan(M_t) || minimum(lee.R⃗ᵥ[t] .- lee.R⃗[t]) < 1e-6
            #         @show t, qₜ, Nₜ
            #         @show minimum(lee.R⃗[t]), maximum(lee.R⃗[t])
            #         @show minimum(lee.R⃗ᵥ[t]), maximum(lee.R⃗ᵥ[t])
            #         @show minimum(lee.R⃗ᵥ[t] .- lee.R⃗[t])
            #         @show N_t, S_t, M_t
            #         @warn "NaN at t=$t, skipping forward update"
            #         break
            #     end

            #     # i. fAggregation → N*, S*, M*
            #     n⃗ₜ              = lee.n⃗[t]
            #     _, _, S, M, _   = fAggregationAggregateUncertainty(params, lee.R⃗ᵥ[t], lee.∂R⃗ᵥ[t], lee.R⃗[t], lee.∂R⃗[t], n⃗ₜ, pₜ)
            #     # ii. Beveridge: N¹ = Nₜ₋₁ + M - S, f¹ = M/(L - N¹)
            #     N⁺¹             = Nₜ₋₁ + M - S
            #     f⃗[t]            = M / (L - N⁺¹ + 1e-8) 

            #     # iii. Update q and check convergence
            #     q¹          = fUpdatedJobFindingRateInverse(f⃗[t], params)  
            #     εᵍ          = abs(N⁺¹ - Nₜ)
            #     Nₜ          = N⁺¹
            #     qₜ          = ωᵍ * q¹ + (1 - ωᵍ) * qₜ
            #     qₜ          = clamp(qₜ, q̲, q̅)

            #     # Diagnostic 
            #     if nᵍ % 25 == 0 && (t <= 3 || t >= T-2 || t % 200 == 0)
            #         @printf "  t=%4d | εᵍ=%.3e | qₜ=%.4f | Nₜ=%.4f | f=%.4f \n" t εᵍ qₜ Nₜ f⃗[t]
            #     end
            #     nᵍ += 1
            # end

            # J2. Single-pass forward step
            n⃗ₜ              = lee.n⃗[t]
            _, Y, S, M, A   = fAggregationAggregateUncertainty(params, lee.R⃗ᵥ[t], lee.∂R⃗ᵥ[t], lee.R⃗[t], lee.∂R⃗[t], n⃗ₜ, pₜ)
            Nₜ              = Nₜ₋₁ + M - S
            f⃗[t]            = M / (L - Nₜ + 1e-8)
            qₜ              = fUpdatedJobFindingRateInverse(f⃗[t], params)
            qₜ              = clamp(qₜ, q̲, q̅)
            if isnan(Nₜ)
                lee.q⃗ⁱᵐᵖ[t] = lee.q⃗[t]
                lee.N⃗ⁱᵐᵖ[t] = Nₜ₋₁
                Nₜᵖʳᵉᵛ      = Nₜ₋₁
                continue
            end
            Nₜᵖʳᵉᵛ          = Nₜ

            # J3. Store realised paths
            lee.q⃗ⁱᵐᵖ[t]     = qₜ
            lee.N⃗ⁱᵐᵖ[t]     = Nₜ
            simu.S⃗[t]       = S
            simu.M⃗[t]       = M
            simu.Y⃗[t]       = Y
            simu.A⃗[t]       = A


            # # (iv) Correct the grid error 
            # n⃗ₒₗₐ                    = lee.n⃗[t]
            # lee.n⃗[t]                = fn⃗(params, pₜ, f⃗[t], qₜ)
            # lee.R⃗[t], lee.∂R⃗[t]     = fRobustSplineAggregateUncertainty(n⃗ₒₗₐ, lee.R⃗[t], lee.n⃗[t])
            # lee.R⃗ᵥ[t], lee.∂R⃗ᵥ[t]   = fRobustSplineAggregateUncertainty(n⃗ₒₗₐ, lee.R⃗ᵥ[t], lee.n⃗[t])
            # lee.R⃗[t]  = clamp.(lee.R⃗[t], x̲, x̅)
            # lee.R⃗ᵥ[t] = clamp.(lee.R⃗ᵥ[t], x̲, x̅)
            # lee.R⃗ᵥ[t] = max.(lee.R⃗ᵥ[t], lee.R⃗[t] .* (1 + 1e-3))   
            # Π̃ₒₗₐ                    = lee.Π̃[t]
            # lee.Π̃[t]                = zeros(Nₓ, length(lee.n⃗[t]))
            # for i in 1:Nₓ
            #     lee.Π̃[t][i,:]       .= Spline1D(n⃗ₒₗₐ, Π̃ₒₗₐ[i,:]; k=3, bc="extrapolate")(lee.n⃗[t])
            # end
            # lee.Π̃ⁱᵐᵖ[t]             = zeros(Nₓ, length(lee.n⃗[t]))

            # J3. Store realised paths
            # _, Y, S, M, A   = fAggregationAggregateUncertainty(params,lee.R⃗ᵥ[t],lee.∂R⃗ᵥ[t],lee.R⃗[t],lee.∂R⃗[t],lee.n⃗[t],pₜ)
            # lee.q⃗ⁱᵐᵖ[t]     = qₜ
            # lee.N⃗ⁱᵐᵖ[t]     = Nₜ
            # simu.S⃗[t]       = S 
            # simu.M⃗[t]       = M 
            # simu.Y⃗[t]       = Y
            # simu.A⃗[t]       = A
        end

        #%% J4. Convergence check and damped update 
        εᴿᵀᴹ            = maximum(abs.(lee.q⃗ .- lee.q⃗ⁱᵐᵖ))
        push!(ϵʰⁱˢᵗ, εᴿᵀᴹ)
        fnPlotConvergenceRTM(ϵʰⁱˢᵗ, lee, ss_endo, params)
        lee.q⃗           .= ωᴿᵀᴹ₁ .* lee.q⃗ⁱᵐᵖ .+ (1 - ωᴿᵀᴹ₁) .* lee.q⃗
        lee.N⃗           .= ωᴿᵀᴹ₂ .* lee.N⃗ⁱᵐᵖ .+ (1 - ωᴿᵀᴹ₂) .* lee.N⃗
        for t in 1:T
            lee.Π̃[t]    = ωᴿᵀᴹ₃ .* lee.Π̃ⁱᵐᵖ[t] .+ (1 - ωᴿᵀᴹ₃) .* lee.Π̃[t]
        end
        lee.εᴿᵀᴹ        = εᴿᵀᴹ

        # Print progress
        @printf "RTM Iteration: %4d | ε = %.6f \n" nᴿᵀᴹ εᴿᵀᴹ
        if nᴿᵀᴹ % 5 == 0
            @save "results/rtm_warm_start.jld2" warm_q=lee.q⃗ warm_N=lee.N⃗ warm_Pi=lee.Π̃ warm_n=lee.n⃗
            println("  💾 Results saved for insurance at iteration $nᴿᵀᴹ")
        end
        nᴿᵀᴹ           += 1
    end

    #%% K. Collect results into the simu structure 
    simu.N⃗ .= lee.N⃗
    simu.q⃗ .= lee.q⃗
    simu.f⃗ .= fUpdatedJobFindingRate.(lee.q⃗, Ref(params))
end

# 4. Live plotting function for RTM convergence
function fnPlotConvergenceRTM(error_history, lee, ss_endo, params)
    
    T       = length(lee.q⃗)
    t_grid  = 1:T

    p1 = plot(error_history, title="RTM error history", xlabel="Iteration", 
            ylabel="sup |q - qⁱᵐᵖ|", yscale=:log10, lw=2, color=:black, 
            legend=false, grid=true)

    p2 = plot(t_grid, lee.q⃗ .- lee.q⃗ⁱᵐᵖ, title="q: predicted − implied", xlabel="t", 
            lw=1, color=:navy, legend=false, grid=true)
    hline!(p2, [0.0], color=:gray, ls=:dash)

    p3 = plot(t_grid, lee.N⃗ .- lee.N⃗ⁱᵐᵖ, title="N: predicted − implied", xlabel="t", 
            lw=1, color=:maroon, legend=false, grid=true)
    hline!(p3, [0.0], color=:gray, ls=:dash)

    plt = plot(p1, p2, p3, layout=(1, 3), size=(1500, 500), margin=5Plots.mm)
    display(plt)
    return plt
end

# Old Krussel and Smith (1998) draft → May continue it if I find time 
# Content: 
# 1A. Collapsed forecasting parameters 
# 1B. Prepare aggregate N grid 

# 1A. Collapsed forecasting parameters 
# function fnCollapsedForecastingParameters(ν₀,νₙ,νₚ,θ₀,θₙ,θₚ)

#     # A. Compute the simplification 
#     𝔼θₙ         = θₙ*νₙ
#     𝔼θₚ         = θₚ+θₙ*νₚ
#     𝔼θ₀         = θ₀+θₙ*ν₀

#     # B. Return values
#     return 𝔼θₙ, 𝔼θₚ, 𝔼θ₀
# end 

# # 1B. Prepare aggregate grids
# function fnAggregateGridsStateSpace!(params,KS,ν₀,νₙ,νₚ,θ₀,θₙ,θₚ)

#     # A. Unpacking business 
#     @unpack p⃗,Ñₙ,Ñₜ = params 

#     # B. Prepare the "collapsed" forecasting parameters
#     KS.𝔼θₙ, KS.𝔼θₚ, KS.𝔼θ₀   = fnCollapsedForecastingParameters(ν₀,νₙ,νₚ,θ₀,θₙ,θₚ)

#     # C. Aggregate employment state space 
#     𝒩⃗̄               = max(3.333,(ν₀+νₚ*maximum(p⃗))/(1-νₙ))
#     𝒩̲⃗               = min(3.220,(ν₀+νₚ*minimum(p⃗))/(1-νₙ))
#     Δ𝒩              = (𝒩⃗̄-𝒩̲⃗)/(Ñₙ-1)
#     KS.𝒩⃗            = collect(𝒩̲⃗:Δ𝒩:𝒩⃗̄)
# end 
