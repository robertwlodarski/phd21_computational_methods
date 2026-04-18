# Important note on notation with aggregate uncertainty 
# 𝒩⃗: caligraphic letters:           State space
# N⃗: standard letters:              simulation  
# Π̃: greek letters with a tilde:    Value functions


function fnSolveAggregateRTM!(params, ss_endo, simu, lee)

    #%% 1. Unpacking
    @unpack P, p⃗, Nₚ, L, δᴿᵀᴹ,  δᵍ      = params
    @unpack p⃗̂, p⃗̂ᵢ, N⃗, q⃗, f⃗, S⃗, M⃗, Y⃗     = simu
    T                                   = length(p⃗̂)

    #%% 2. Initialise predicted paths
    # Use SS values as starting conjecture for {qₜ⁽⁰⁾, Nₜ⁽⁰⁾}ᵀ
    # Initialise Π̃ᵖ with SS value functions via fVFI!

    #%% 3. RTM outer loop
    εᴿᵀᴹ    = Inf
    nᴿᵀᴹ    = 1
    while εᴿᵀᴹ > δᴿᵀᴹ

        #%% 3.1. Backward solution (t = T → 1)
        for t in T:-1:1

            # A. Retrieve current period's exogenous state and nth-iteration conjectured price
            pₜ  = p⃗̂[t]
            pᵢₜ = p⃗̂ᵢ[t]
            qₜ  = lee.q⃗ᵖ[t]
            
            # B. Construct 𝔼ₜ[Π̃ₜ₊₁] via RTM matching
            # For the realised p′: use Π̃ʳ[t+1] directly from this backward pass
            # For each counterfactual p′ ≠ pₜ₊₁: find τ s.t. N⃗ᵖ[τ] ≈ N⃗ᵖ[t+1] and p⃗̂ᵢ[τ] == p′
            for pⱼ in 1:Nₚ

                # i. Find matching period τ (the RTM "repeated transition" step)

                # ii. Borrow Π̃ᵖ[τ] for this counterfactual state

            end

            # C. Assemble expected value: 𝔼ₜ[Π̃ₜ₊₁] = Σⱼ P[pᵢ,j] × Π̃ₜ₊₁(·; pⱼ)

            # D. Single Bellman update and policy extraction

            # E. Extract and store policy functions R⃗ₜ, R⃗ᵥₜ, n⃗ₜ into lee.R⃗ᵀ[t] etc.

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
