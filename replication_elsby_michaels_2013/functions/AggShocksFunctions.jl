# Important note on notation with aggregate uncertainty 
# 𝒩⃗: caligraphic letters:           State space
# N⃗: standard letters:              Simulation  
# Π̃: greek letters with a tilde:    Value functions

function fnSolveAggregateRTM!(params, ss_endo, simu, lee)

    #%% 1. Unpacking business 
    @unpack P, p⃗, Nₚ, L, δᴿᵀᴹ, δᵍ, x⃗, α, β, c = params
    @unpack p⃗̂, p⃗̂ᵢ, N⃗, q⃗, f⃗, S⃗, M⃗, Y⃗     = simu
    T                                   = length(p⃗̂)

    #%% 2. Initialise predicted paths
    # A. Use SS values as starting conjecture for {qₜ⁽⁰⁾, Nₜ⁽⁰⁾}ᵀ
    lee.q⃗   .= ss_endo.q̂
    lee.N⃗   .= ss_endo.N .* (1 .+ 1e-2 *randn(T))
    # Initialise Π̃ with SS value functions via fVFI!
    ℑˢˢ = [Spline1D(ss_endo.n⃗, ss_endo.Π[i,:]; k=3, bc="extrapolate") for i in 1:Nₓ]
    @. f⃗ = fUpdatedJobFindingRate(lee.q⃗ᵖ, Ref(params))
    for t in 1:T
        lee.n⃗ᵀ[t]           = fn⃗(params, p⃗̂[t], f⃗[t], lee.q⃗ᵖ[t])
        for i in 1:Nₓ
            lee.Π̃[t][i,:]  .= ℑˢˢ[i](lee.n⃗ᵀ[t])
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
            qₜ  = lee.q⃗ᵖ[t]
            
            # B. Construct 𝔼ₜ[Π̃ₜ₊₁] via RTM matching
            Π̃ₜ₊₁        = Vector{Matrix{Float64}}(undef, Nₚ)
            # Realised state handled directly
            Π̃ₜ₊₁[pᵢₜ] = t == T ? lee.Π̃[1] : lee.Π̃[t+1]
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
                ωτ              = (lee.N⃗[τ̅] - lee.N⃗[t+1]) / (lee.N⃗[τ̅] - lee.N⃗[τ̲])
                ωτ              = clamp(ωτ, 0.0, 1.0)

                # ii. Borrow and interpolate Π̃ for this counterfactual state
                Π̃ₜ₊₁[pⱼ]        = ωτ .* lee.Π̃[τ_low] .+ (1 - ωτ) .* lee.Π̃[τ_high]

            end

            # C. Assemble expected value: 𝔼ₜ[Π̃ₜ₊₁] = Σⱼ P[pᵢ,j] × Π̃ₜ₊₁(·; pⱼ)
            𝔼Π̃ = sum(P[pᵢₜ, pⱼ] .* Π̃ₜ₊₁[pⱼ] for pⱼ in 1:Nₚ)
            
            # D. Single Bellman update and policy extraction
            W           = fW(params, pₜ, f⃗[t], qₜ, lee.n⃗[t])
            Πᶠˡᵒʷ       = pₜ .* x⃗ .* (lee.n⃗[t]' .^ α) .- W .* lee.n⃗[t]'
            Πᶜ          = Πᶠˡᵒʷ .+ β .* 𝔼Π̃
            for i in 1:Nₓ
                # Fire 
                valᶠ, idᶠ   = findmax(view(Πᶜ, i, :))
                vᶠ[i]       = valᶠ
                nᶠ[i]       = lee.n⃗[t][idᶠ]
                # Hire 
                valʰ, idʰ   = findmax(view(Πᶜ, i, :) .- (c / qₜ) .* lee.n⃗[t])
                vʰ[i]       = valʰ
                nʰ[i]       = lee.n⃗[t][idʰ]
            end
            Πᶠ              .= vᶠ .* (nᶠ .< lee.n⃗[t]') .- 1e8 .* (nᶠ .>= lee.n⃗[t]')
            Πʰ              .= (vʰ .+ c / qₜ .* lee.n⃗[t]') .* (nʰ .> lee.n⃗[t]') .- 1e8 .* (nʰ .<= lee.n⃗[t]')
            lee.Π̃[t]        .= max.(Πᶠ, max.(Πʰ, Πᶜ))

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
            lee.R⃗ᵥ[t]   = min.(ℑᴿⱽ(lee.n⃗[t]), x̅)
            lee.∂R⃗ᵥ[t]  = Dierckx.derivative(ℑᴿⱽ, lee.n⃗[t])
        end 

        #%% 3.2. Forward simulation (t = 1 → T)
        for t in 1:T

            # A. Retrieve pₜ and Nₜ₋₁

            # B. Solve equilibrium q loop (reuse fEqResidual logic)
            # Inner loop: given q guess, compute f, run fAggregation,
            # update N from Beveridge curve N = Nₜ₋₁ + M - S, iterate until convergence
            qₜ      = lee.q⃗ᵖ[t]     # initialise from predicted path
            Nₜ      = lee.N⃗ᵖ[t]
            εᵍ      = Inf
            while εᵍ > δᵍ

                # i. Compute f from q

                # ii. fAggregation → N*, S*, M*

                # iii. Beveridge: N¹ = Nₜ₋₁ + M - S, f¹ = M/(L - N¹)

                # iv. Update q and check convergence

            end  

            # C. Store realised paths
            lee.q⃗ʳ[t]  = qₜ
            lee.N⃗ʳ[t]  = Nₜ
            simu.S⃗[t]  = x # separations
            simu.M⃗[t]  = x # matches
            simu.Y⃗[t]  = x # output

        end

        #%% 3.3. Convergence check and damped update
        εᴿᵀᴹ            = maximum(abs.(lee.q⃗ʳ .- lee.q⃗ᵖ))
        lee.q⃗ᵖ         .= params.ω .* lee.q⃗ʳ .+ (1 - params.ω) .* lee.q⃗ᵖ
        lee.N⃗ᵖ         .= params.ω .* lee.N⃗ʳ .+ (1 - params.ω) .* lee.N⃗ᵖ
        for t in 1:T
            lee.Π̃ᵖ[t]  .= params.ω .* lee.Π̃ʳ[t] .+ (1 - params.ω) .* lee.Π̃ᵖ[t]
        end
        lee.εᴿᵀᴹ        = εᴿᵀᴹ

        # Print progress
        @printf "RTM Iteration: %4d | ε = %.6f \n" nᴿᵀᴹ εᴿᵀᴹ
        nᴿᵀᴹ           += 1
    end

    #%% 4. Collect results into simu
    simu.N⃗ .= lee.N⃗ʳ
    simu.q⃗ .= lee.q⃗ʳ
    simu.f⃗ .= fUpdatedJobFindingRate.(lee.q⃗ʳ, Ref(params))
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
