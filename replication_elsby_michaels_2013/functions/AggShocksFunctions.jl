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
    𝔼x      = (1 / p̄ₓ) * x̲^ξ * (ξ /(ξ - 1)) * (R⃗.^(-ξ+1)-R⃗ᵥ.^(-ξ+1)) ./ (𝐆R⃗ᵥ .- 𝐆R⃗)

    # C. Compute aggregate values 
    N       = fSimpsonRule(n⃗ .*𝐡n⃗, n⃗)                   # Employed
    Y       = fSimpsonRule(𝔼x .* n⃗.^α .*𝐡n⃗,n⃗)           # Production
    S       = λ * fSimpsonRule((1 .- 𝐇n⃗).*𝐆R⃗,n⃗)         # Separations 
    M       = λ * fSimpsonRule(𝐇n⃗.*(1 .- 𝐆R⃗ᵥ),n⃗)        # Matches 
    A       = fSimpsonRule(p .* 𝔼x .* n⃗.^(α-1).*𝐡n⃗,n⃗)   # Total marginal product of labour  
    return N, Y, S, M, A 
end 

# 2. Solve using the repeated transition method 
function fnSolveAggregateRTM!(params, ss_endo, simu, lee)

    #%% 1. Unpacking business 
    @unpack P, p⃗, Nₚ, L, δᴿᵀᴹ, δᵍ, x⃗, α, β, c, Nₓ, ωᵍ, ωᴿᵀᴹ₁, ωᴿᵀᴹ₂, ωᴿᵀᴹ₃ = params
    @unpack p⃗̂, p⃗̂ᵢ, N⃗, q⃗, f⃗, S⃗, M⃗, Y⃗     = simu
    T                                   = length(p⃗̂)

    #%% 2. Initialise predicted paths
    # A. Use SS values as starting conjecture for {qₜ⁽⁰⁾, Nₜ⁽⁰⁾}ᵀ
    lee.q⃗   .= ss_endo.q̂
    lee.N⃗   .= ss_endo.N .* (1 .+ 1e-2 *randn(T))
    # Initialise Π̃ with SS value functions via fVFI!
    ℑˢˢ = [Spline1D(ss_endo.n⃗, ss_endo.Π[i,:]; k=3, bc="extrapolate") for i in 1:Nₓ]
    @. f⃗ = fUpdatedJobFindingRate(lee.q⃗, Ref(params))
    for t in 1:T
        lee.n⃗[t]            = fn⃗(params, p⃗̂[t], f⃗[t], lee.q⃗[t])
        for i in 1:Nₓ
            lee.Π̃[t][i,:]  .= ℑˢˢ[i](lee.n⃗[t])
        end
    end

    #%% 3. RTM outer loop
    εᴿᵀᴹ    = Inf
    nᴿᵀᴹ    = 1
    vᶠ      = zeros(Nₓ)
    nᶠ      = zeros(Nₓ)
    vʰ      = zeros(Nₓ)
    nʰ      = zeros(Nₓ)
    while εᴿᵀᴹ > δᴿᵀᴹ

        #%% 3.1. Backward solution (t = T → 1)
        for t in T:-1:1
            Nₙₜ = length(lee.n⃗[t])
            Πᶠ  = zeros(Nₓ, Nₙₜ)
            Πʰ  = zeros(Nₓ, Nₙₜ)

            # A. Retrieve current period's exogenous state and nth-iteration conjectured price
            pₜ  = p⃗̂[t]
            pᵢₜ = p⃗̂ᵢ[t]
            qₜ  = lee.q⃗[t]
            
            # B. Construct 𝔼ₜ[Π̃ₜ₊₁] via RTM matching
            Π̃ₜ₊₁        = Vector{Matrix{Float64}}(undef, Nₚ)
            # Realised state handled directly
            Π̃ₜ₊₁[pᵢₜ] = t == T ? lee.Π̃[1] : lee.Π̃[t+1]
            # For the realised p′: use Π̃ʳ[t+1] directly from this backward pass
            # For each counterfactual p′ ≠ pₜ₊₁: find τ s.t. N⃗ᵖ[τ] ≈ N⃗ᵖ[t+1] and p⃗̂ᵢ[τ] == p′
            Threads.@threads for pⱼ in filter(!=(pᵢₜ), 1:Nₚ)

                # i. Find matching period τ via interpolation on N
                τ⃗ᶜ              = findall(p⃗̂ᵢ .== pⱼ)
                N⃗ᶜ              = lee.N⃗[τ⃗ᶜ]
                Nᵉᶠᶠ            = t == T ? lee.N⃗[1] : lee.N⃗[t+1]
                idx             = sortperm(abs.(N⃗ᶜ .- Nᵉᶠᶠ))
                τ̲               = τ⃗ᶜ[idx[1]]
                τ̅               = τ⃗ᶜ[idx[2]]
                ωτ              = (lee.N⃗[τ̅] - Nᵉᶠᶠ) / (lee.N⃗[τ̅] - lee.N⃗[τ̲])
                ωτ              = clamp(ωτ, 0.0, 1.0)

                # ii. Borrow and interpolate Π̃ for this counterfactual state
                Π̃ₜ₊₁[pⱼ]        = ωτ .* lee.Π̃[τ̲] .+ (1 - ωτ) .* lee.Π̃[τ̅]

            end

            # C. Assemble expected value: 𝔼ₜ[Π̃ₜ₊₁] = Σⱼ P[pᵢ,j] × Π̃ₜ₊₁(·; pⱼ)
            𝔼Π̃ = sum(P[pᵢₜ, pⱼ] .* Π̃ₜ₊₁[pⱼ] for pⱼ in 1:Nₚ)
            
            # D. Single Bellman update and policy extraction
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
            lee.Π̃ⁱᵐᵖ[t]        .= max.(Πᶠ, max.(Πʰ, Πᶜ))

            # E. Extract and store policy functions R⃗ₜ, R⃗ᵥₜ, n⃗ₜ into lee.R⃗ᵀ[t] etc.
            ℑⁿˡ         = Spline1D(x⃗, nᶠ; k=3, bc="extrapolate")
            ℑⁿʰ         = Spline1D(x⃗, nʰ; k=3, bc="extrapolate")
            𝓃₁          = findlast(nᶠ .< nʰ[Nₓ])
            𝓃₂          = findfirst(nʰ .> nᶠ[1])

            # E1. Firing threshold 
            𝕟ᴿ          = nᶠ[1:𝓃₁]
            𝒾ᴿ          = unique(i -> 𝕟ᴿ[i], 1:length(𝕟ᴿ))
            𝕩ᴿ          = x⃗[1:𝓃₁]
            ℑᴿ          = Spline1D(𝕟ᴿ[𝒾ᴿ], 𝕩ᴿ[𝒾ᴿ]; k=3, bc="extrapolate")
            lee.R⃗[t]    = ℑᴿ(lee.n⃗[t])
            lee.∂R⃗[t]   = Dierckx.derivative(ℑᴿ, lee.n⃗[t])

            # E2. Hiring threshold 
            𝕟ᴿⱽ         = nʰ[𝓃₂:end]
            𝒾ᴿⱽ         = unique(i -> 𝕟ᴿⱽ[i], 1:length(𝕟ᴿⱽ))
            𝕩ᴿⱽ         = x⃗[𝓃₂:end]
            ℑᴿⱽ         = Spline1D(𝕟ᴿⱽ[𝒾ᴿⱽ], 𝕩ᴿⱽ[𝒾ᴿⱽ]; k=3, bc="extrapolate")
            lee.R⃗ᵥ[t]   = min.(ℑᴿⱽ(n⃗ₜ), x̅)
            lee.∂R⃗ᵥ[t]  = Dierckx.derivative(ℑᴿⱽ, n⃗ₜ)
        end 

        #%% 3.2. Forward simulation (t = 1 → T)
        for t in 1:T

            # A. Retrieve pₜ and Nₜ₋₁
            pₜ      = p⃗̂[t]
            Nₜ₋₁    = t == 1 ? lee.N⃗[T] : lee.N⃗[t-1]

            # B. Solve equilibrium q loop (reuse fEqResidual logic)
            # Inner loop: given q guess, compute f, run fAggregation,
            # update N from Beveridge curve N = Nₜ₋₁ + M - S, iterate until convergence
            qₜ          = lee.q⃗[t]     
            Nₜ          = lee.N⃗[t]
            lee.n⃗[t]    = fn⃗(params, pₜ, f⃗[t], qₜ)
            εᵍ          = Inf
            while εᵍ > δᵍ

                # i. fAggregation → N*, S*, M*
                n⃗ₜ              = lee.n⃗[t]
                _, _, S, M, _   = fAggregationAggregateUncertainty(params,view(lee.R⃗ᵥ[t][:]),view(lee.∂R⃗ᵥ[t][:]),view(lee.R⃗[t][:]),view(lee.∂R⃗[t]),n⃗ₜ,pₜ)

                # ii. Beveridge: N¹ = Nₜ₋₁ + M - S, f¹ = M/(L - N¹)
                N⁺¹             = Nₜ₋₁ + M - S
                f⃗[t]            = M / (L - N⁺¹ + 1e-8) 

                # iii. Update q and check convergence
                q¹          = fUpdatedJobFindingRateInverse(f⃗[t], params)  
                εᵍ          = abs(N⁺¹ - Nₜ)
                Nₜ          = N⁺¹
                qₜ          = ωᵍ * q¹ + (1 - ωᵍ) * qₜ
                lee.n⃗[t]    = fn⃗(params, pₜ, f⃗[t], qₜ)

            end  

            # C. Store realised paths
            _, Y, S, M, A   = fAggregationAggregateUncertainty(params,view(lee.R⃗ᵥ[t][:]),view(lee.∂R⃗ᵥ[t][:]),view(lee.R⃗[t][:]),view(lee.∂R⃗[t]),view(lee.n⃗[t][:]),pₜ)
            lee.q⃗ⁱᵐᵖ[t]     = qₜ
            lee.N⃗ⁱᵐᵖ[t]     = Nₜ
            simu.S⃗[t]       = S 
            simu.M⃗[t]       = M 
            simu.Y⃗[t]       = Y
            simu.A⃗[t]       = A
        end

        #%% 3.3. Convergence check and damped update 
        εᴿᵀᴹ            = maximum(abs.(lee.q⃗ .- lee.q⃗ⁱᵐᵖ))
        lee.q⃗           .= ωᴿᵀᴹ₁ .* lee.q⃗ⁱᵐᵖ .+ (1 - ωᴿᵀᴹ₁) .* lee.q⃗
        lee.N⃗           .= ωᴿᵀᴹ₂ .* lee.N⃗ⁱᵐᵖ .+ (1 - ωᴿᵀᴹ₂) .* lee.N⃗
        for t in 1:T
            lee.Π̃[t]    .= ωᴿᵀᴹ₃ .* lee.Π̃ⁱᵐᵖ[t] .+ (1 -  ωᴿᵀᴹ₃) .* lee.Π̃[t]
        end
        lee.εᴿᵀᴹ        = εᴿᵀᴹ

        # Print progress
        @printf "RTM Iteration: %4d | ε = %.6f \n" nᴿᵀᴹ εᴿᵀᴹ
        nᴿᵀᴹ           += 1
    end

    #%% 4. Collect results into the simu structure 
    simu.N⃗ .= lee.N⃗
    simu.q⃗ .= lee.q⃗
    simu.f⃗ .= fUpdatedJobFindingRate.(lee.q⃗, Ref(params))
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
