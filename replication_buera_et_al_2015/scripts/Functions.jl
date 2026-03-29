# Content:

# A. VFI components
# 1. Utility function 
# 2. Static policy functions 
# 3. Initial VFI guesses 
# 4. Solve for assets 
# 5. VFI

# B. Aggregating 
# 1. Job destruction
# 2. Update labour market (unemployment and matching)
# 3. Compute flows 
# 4. Forward iteration
# 5. Aggregate states 

# C. Solving  
# 1A. Logistic function for bounding 
# 1B. Inverse logistic for bounding
# 2. Combined steady state sesiduals
# 3. Solve for the steady state 

# -----------------------------------------------
# A. VFI components 
# -----------------------------------------------

# 1. Utility function 
function fnUtility(C, params)

    # A. Unpacking business 
    @unpack σ, c̲ = params

    # B. Define a pure scalar function with the short-circuiting ternary operator
    u(c) = c <= 0.0 ? (c̲^(1 - σ)) / (1 - σ) : (c^(1 - σ)) / (1 - σ)

    # C. Broadcast the custom function over C (works for scalars and arrays!)
    return u.(C)
end

# 2. Static policy functions  
function fnStaticPolicies!(params, endo)

    # A. Unpacking business 
    @unpack α, δ, θ, A, z⃗, λ, a⃗ = params

    # B. Simplifying notation 
    A1 = α * A * z⃗ ./ (endo.rₜ + δ)
    A2 = θ / α * (endo.rₜ + δ) / endo.wₜ
    A3 = A2^(θ)

    # C. Getting the capital 
    kˣ = (A1 .* A3) .^ (1 / (1 - α - θ))
    endo.𝐤 .= min.(kˣ, λ .* a⃗')
    endo.𝕀ᶜ .= (kˣ .>= λ .* a⃗')

    # D. Getting the labour 
    endo.𝐥 .= (θ .* A .* z⃗ .* endo.𝐤 .^ α ./ endo.wₜ) .^ (1 / (1 - θ))

    # E. Compute the associated profit 
    endo.Π .= A .* z⃗ .* (endo.𝐤) .^ α .* (endo.𝐥) .^ θ .- endo.wₜ .* endo.𝐥 .- (endo.rₜ + δ) .* endo.𝐤
end

# 3. Initial VFI guesses 
function fnInitialVFIGuess!(params, endo)

    # A. Unpacking business 
    @unpack β, a⃗ = params

    if endo.𝐕[1, 1] == 0.0
        # B. Workers VF 
        endo.𝐕ᵂ .= 0.75 * (1 - β)^(-1) * fnUtility(endo.rₜ .* a⃗' .- endo.τₜ .+ endo.wₜ, params)

        # C. Entrepreneurs VF 
        endo.𝐕ᴱ .= 0.75 * (1 - β)^(-1) * fnUtility(endo.rₜ .* a⃗' .- endo.τₜ .+ endo.Π, params)

        # D. Overall VF 
        endo.𝐕 .= max.(endo.𝐕ᵂ, endo.𝐕ᴱ)
        endo.𝐨 .= (endo.𝐕ᵂ .<= endo.𝐕ᴱ)
    else 
    end 

end

# 4. Solve for assets 
function fnFindAssets(iz, ia, exp_spline, params, endo)

    # A. Unpack parameters 
    @unpack c̲, a⃗, β = params

    # B. Cash variables
    Cash_w = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ
    Cash_e = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz, ia] - endo.τₜ

    # C. Bounds 
    Lower = a⃗[1]
    Upper_w = min(Cash_w - c̲, a⃗[end])
    Upper_e = min(Cash_e - c̲, a⃗[end])

    # D. Computation of assets: Worker 
    if Upper_w <= Lower
        A_w         = a⃗[1]
        𝐕ᵂ          = fnUtility(Cash_w - a⃗[1], params) + exp_spline(a⃗[1])
        𝔼𝐕ᵂ         = exp_spline(A_w)
    else
        ℜʷ          = optimize(ap -> -(fnUtility(Cash_w - ap, params) + exp_spline(ap)),a⃗[1],Upper_w)
        A_w         = Optim.minimizer(ℜʷ)
        𝐕ᵂ          = - Optim.minimum(ℜʷ)
        𝔼𝐕ᵂ         = exp_spline(A_w)
    end

    # E. Computation of assets: Entrepreneurs 
    if Upper_e <= Lower
        A_e         = a⃗[1]
        𝐕ᴱ          = fnUtility(Cash_e - a⃗[1], params) + exp_spline(a⃗[1])
        𝔼𝐕ᴱ         = exp_spline(A_e)
    else
        ℜᴱ          = optimize(ap -> -(fnUtility(Cash_e - ap, params) + exp_spline(ap)),a⃗[1],Upper_e)
        A_e         = Optim.minimizer(ℜᴱ)
        𝐕ᴱ          = - Optim.minimum(ℜᴱ)
        𝔼𝐕ᴱ         = exp_spline(A_e)
    end
    return A_w, A_e, 𝐕ᵂ, 𝐕ᴱ, 𝔼𝐕ᵂ, 𝔼𝐕ᴱ
end

# 5. Value function iteration 
function fnVFI!(params, endo)

    # A. Unpacking business
    @unpack ψ, μ⃗, δᵛᶠⁱ, 𝒾̄ᵛᶠⁱ, z⃗, c̲, a⃗, β, λᵛᶠⁱ,Nᶻ,Nᵃ = params

    # B. Static policies
    fnStaticPolicies!(params, endo)

    # C. Initial guess 
    fnInitialVFIGuess!(params, endo)

    # D. Prepare the VFI loop 
    εᵛᶠⁱ    = 1.0
    𝓃ᵛᶠⁱ    = 1
    𝔼𝐕ᵂ     = zeros(Nᶻ,Nᵃ) 
    𝔼𝐕ᴱ     = zeros(Nᶻ,Nᵃ)

    while (εᵛᶠⁱ > δᵛᶠⁱ && 𝓃ᵛᶠⁱ < 𝒾̄ᵛᶠⁱ)

        # D1. Interpolate the expected value properly 
        ℑᴱ                      = [Spline1D(a⃗, view(endo.𝐕ᴱ, iz, :); k=1, bc="nearest") for iz in eachindex(z⃗)]
        ℑᵂ                      = [Spline1D(a⃗, view(endo.𝐕ᵂ, iz, :); k=1, bc="nearest") for iz in eachindex(z⃗)]
        ℑ̄(iz, a)                = max(ℑᴱ[iz](a), ℑᵂ[iz](a))
        ℑ(iz, a)                = β * (ψ * ℑ̄(iz, a) + (1 - ψ) * sum(μ⃗[jz] * ℑ̄(jz, a) for jz in eachindex(z⃗)))

        # D2. Update the expected VF and open the loop for productivity 
        endo.𝔼𝐕 .= ψ .* endo.𝐕 .+ (1 - ψ) .* (μ⃗' * endo.𝐕)
        𝐕ᵖʳᵉᵛ = copy(endo.𝐕)
        Threads.@threads for ia in eachindex(a⃗)
            @inbounds for iz in eachindex(z⃗)
                # B1. Solve for each a and z 
                # I. Find assets
                ℑᶻ(a) = ℑ(iz, a)
                endo.𝐚ʷ[iz, ia], endo.𝐚ᵉ[iz, ia], endo.𝐕ᵂ[iz, ia], endo.𝐕ᴱ[iz, ia], 𝔼𝐕ᵂ[iz, ia], 𝔼𝐕ᴱ[iz, ia] = fnFindAssets(iz, ia, ℑᶻ, params, endo)

                # II. VFs for different occupations
                endo.𝐜ʷ[iz, ia] = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ - endo.𝐚ʷ[iz, ia]
                endo.𝐜ᵉ[iz, ia] = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz, ia] - endo.τₜ - endo.𝐚ᵉ[iz, ia]

                # III. Occupational decision and updated VF
                endo.𝐨[iz, ia] = (endo.𝐕ᴱ[iz, ia] >= endo.𝐕ᵂ[iz, ia])
                endo.𝐚[iz, ia] = endo.𝐨[iz, ia] ? endo.𝐚ᵉ[iz, ia] : endo.𝐚ʷ[iz, ia]
                endo.𝐜[iz, ia] = max(endo.𝐨[iz, ia] ? endo.𝐜ᵉ[iz, ia] : endo.𝐜ʷ[iz, ia], c̲)
                endo.𝐕[iz, ia] = max(endo.𝐕ᴱ[iz, ia], endo.𝐕ᵂ[iz, ia])
            end
        end

        # D2. Error and update 
        εᵛᶠⁱ = maximum(abs.(𝐕ᵖʳᵉᵛ .- endo.𝐕))
        𝓃ᵛᶠⁱ += 1
        # Add this heartbeat print statement!
        # if 𝓃ᵛᶠⁱ % 50 == 0
        #     println("→→→→ VFI Iteration: $𝓃ᵛᶠⁱ, εᵛᶠⁱ = $(round(εᵛᶠⁱ, digits=5))")
        # end
    end
    # println("→→→→ VFI done: 𝓃 = $𝓃ᵛᶠⁱ, ε = $εᵛᶠⁱ")
end

# -----------------------------------------------
# B. Aggregating 
# -----------------------------------------------

# 1. Job destruction
function fnJobDestruction!(params, endo)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, z⃗, a⃗, l⃗ = params

    # B. Initialise 
    endo.S = 0.0
    endo.D = 0.0

    # C. Loop 
    for iz in eachindex(z⃗)
        for ia in eachindex(a⃗)
            for il in eachindex(l⃗)

                # C1. Compute key elements 
                Mass = endo.g[iz, ia, il, 1]
                LabDemand = endo.𝐨[iz, ia] ? endo.𝐥[iz, ia] : 0.0
                JobsDestroyed = max(0.0, l⃗[il] - LabDemand)
                Switchers = (l⃗[il] > 0) && (endo.𝐨[iz, ia] == false)

                # C2. Start adding up 
                endo.S += Mass * Switchers
                endo.D += Mass * JobsDestroyed
            end
        end
    end
    endo.JD = endo.S + endo.D
end

# 2. Update labour market (unemployment and matching)
function fnUpdateLabourMarket!(params, endo)

    # A. Unpacking business 
    @unpack γ = params

    # B. Compute the key values
    endo.U = sum(@views(endo.g[:, :, :, 2]))
    endo.M = γ * (endo.U + endo.JD)
    endo.W = sum(@views(endo.g[:, :, 1, 1]))
    endo.E = 1.0 - endo.U - endo.W
end

# 3. Compute flows 
function fnComputeFlows(iz, ia, il, iu, params, endo)

    # A. Unpacking business
    @unpack l⃗ = params
    # Legend: 
    # iu    = 1 -> employed, entrepreneur 
    # iu    = 2 -> unemployed 
    # E     = Entrepreneurs 
    # W     = Employed workers 
    # U     = Unemployed workers

    # B. Flow indicators 
    JFR = min(1.0, max(endo.M / (endo.U + endo.JD), 0.0))
    JDR = min(1.0, max(0.0, endo.JD / endo.W))
    WU = 0.0 + (l⃗[il] == l⃗[1]) * JDR * (1 - JFR) * (iu == 1) * (endo.𝐨[iz, ia] == false)
    WE = 0.0 + (l⃗[il] == l⃗[1]) * (endo.𝐨[iz, ia] == true) * (iu == 1)
    EU = 0.0 + (l⃗[il] > l⃗[1]) * (endo.𝐨[iz, ia] == false) * (iu == 1) * (1 - JFR)
    EW = 0.0 + (l⃗[il] > l⃗[1]) * (endo.𝐨[iz, ia] == false) * (iu == 1) * JFR
    UE = 0.0 + (endo.𝐨[iz, ia] == true) * (iu == 2)
    UW = 0.0 + (endo.𝐨[iz, ia] == false) * (iu == 2) * JFR
    EE = 0.0 + (l⃗[il] > l⃗[1]) * (endo.𝐨[iz, ia] == true) * (iu == 1)
    WW = 0.0 + (l⃗[il] == l⃗[1]) * (endo.𝐨[iz, ia] == false) * (iu == 1) * ((1 - JDR) + (JDR * JFR))
    UU = 0.0 + (endo.𝐨[iz, ia] == false) * (iu == 2) * (1 - JFR)

    # C. Combine indicators to account for flows that matter
    F1t2 = WU + EU           # To unemployment
    F2t1 = UW + UE           # From unemployment 
    F2t2 = UU
    F1t1 = EE + WW + WE + EW
    return F1t2, F2t1, F2t2, F1t1
end

# 4. Forward iteration 
function fnForwardIteration!(params, endo)

    # A. Unpacking business 
    @unpack z⃗, a⃗, l⃗, δᵈⁱˢᵗ, ψ, μ⃗, Nᶻ, Nᵃ, Nˡ, Nᵘ = params

    # B. Starting the loop for a PDF
    εᵈⁱˢᵗ   = 1.0
    gⁿᵉˣᵗ   = zeros(Nᶻ, Nᵃ, Nˡ, Nᵘ)
    𝓃ᵈⁱˢᵗ   = 1
    B       = zeros(Nᵃ, Nˡ, 2)

    # C. Precompute the invariant elements 
    # I. Initialisation 
    ibᵃ     = zeros(Int, Nᶻ, Nᵃ)
    iuᵃ     = zeros(Int, Nᶻ, Nᵃ)
    wᵃ      = zeros(Float64, Nᶻ, Nᵃ)
    ibˡ     = zeros(Int, Nᶻ, Nᵃ)
    iuˡ     = zeros(Int, Nᶻ, Nᵃ)
    wˡ      = zeros(Float64, Nᶻ, Nᵃ)

    # II. Loop 
    @inbounds for ia in eachindex(a⃗)
        for iz in eachindex(z⃗)

            # (a) Next asset mass 
            Aₜ          = endo.𝐚[iz, ia]
            ibᵃ[iz,ia]  = clamp(searchsortedlast(a⃗, Aₜ), 1, length(a⃗) - 1)
            iuᵃ[iz,ia]  = ibᵃ[iz,ia] + 1
            wᵃ[iz,ia]   = clamp((a⃗[iuᵃ[iz,ia]] - Aₜ) / (a⃗[iuᵃ[iz,ia]] - a⃗[ibᵃ[iz,ia]]), 0.0, 1.0)
            # (b) Labour mass 
            Lₜ          = endo.𝐨[iz, ia] ? endo.𝐥[iz, ia] : 0.0
            ibˡ[iz,ia]  = clamp(searchsortedlast(l⃗, Lₜ), 1, length(l⃗) - 1)
            iuˡ[iz,ia]  = ibˡ[iz,ia] + 1
            wˡ[iz,ia]   = clamp((l⃗[iuˡ[iz,ia]] - Lₜ) / (l⃗[iuˡ[iz,ia]] - l⃗[ibˡ[iz,ia]]), 0.0, 1.0)
        end 
    end 

    # D. Begin the loop 
    while (εᵈⁱˢᵗ > δᵈⁱˢᵗ)

        # E. Aggregate states 
        fnJobDestruction!(params, endo)         # Getting JD, S 
        fnUpdateLabourMarket!(params, endo)     # Getting U, M, W, E  

        @inbounds for iu in 1:2 
            for il in eachindex(l⃗) 
                for ia in eachindex(a⃗)
                    for iz in eachindex(z⃗)
                        # I. Common terms
                        Mass = endo.g[iz, ia, il, iu]
                        if Mass == 0.0
                            continue
                        end

                        # II. Update the flows 
                        # (c) Get flows that matter
                        f¹², f²¹, f²², f¹¹ = fnComputeFlows(iz, ia, il, iu, params, endo)
                        FlowU = f¹² + f²²
                        FlowO = f²¹ + f¹¹

                        # III. With switching productivity
                        # To unemployment (iu=2)
                        B[ibᵃ[iz,ia], ibˡ[iz,ia], 2]    += Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowU
                        B[ibᵃ[iz,ia], iuˡ[iz,ia], 2]    += Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowU
                        B[iuᵃ[iz,ia], ibˡ[iz,ia], 2]    += Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowU
                        B[iuᵃ[iz,ia], iuˡ[iz,ia], 2]    += Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowU
                        # To employment/entrepreneurship (iu=1)
                        B[ibᵃ[iz,ia], ibˡ[iz,ia], 1]    += Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowO
                        B[ibᵃ[iz,ia], iuˡ[iz,ia], 1]    += Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowO
                        B[iuᵃ[iz,ia], ibˡ[iz,ia], 1]    += Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowO
                        B[iuᵃ[iz,ia], iuˡ[iz,ia], 1]    += Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowO

                        # IV. The same productivity
                        # (c) Update (asset, labour) = (BB, BU, UB, UU)
                        # To unemployment 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], ibˡ[iz,ia], 2]    += ψ * Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowU
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], iuˡ[iz,ia], 2]    += ψ * Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], ibˡ[iz,ia], 2]    += ψ * Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], iuˡ[iz,ia], 2]    += ψ * Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowU
                        # To others 
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], ibˡ[iz,ia], 1]    += ψ * Mass * wᵃ[iz,ia] * wˡ[iz,ia] * FlowO
                        gⁿᵉˣᵗ[iz, ibᵃ[iz,ia], iuˡ[iz,ia], 1]    += ψ * Mass * wᵃ[iz,ia] * (1 - wˡ[iz,ia]) * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], ibˡ[iz,ia], 1]    += ψ * Mass * (1 - wᵃ[iz,ia]) * wˡ[iz,ia] * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ[iz,ia], iuˡ[iz,ia], 1]    += ψ * Mass * (1 - wᵃ[iz,ia]) * (1 - wˡ[iz,ia]) * FlowO
                    end
                end
            end
        end

        # IV. Add the mass to the updated productivity states
        @inbounds for iu in 1:2, il in 1:Nˡ ,ia in 1:Nᵃ, iz in 1:Nᶻ
            gⁿᵉˣᵗ[iz, ia, il, iu] += (1 - ψ) * μ⃗[iz] * B[ia, il, iu]
        end
        
        # V. Updating business 
        εᵈⁱˢᵗ = maximum(abs(gⁿᵉˣᵗ[i] - endo.g[i]) for i in eachindex(gⁿᵉˣᵗ))
        endo.g .= gⁿᵉˣᵗ
        fill!(gⁿᵉˣᵗ, 0.0)
        fill!(B, 0.0)
        if 𝓃ᵈⁱˢᵗ % 100 == 0
            println("→→→ Dist. iteration: $𝓃ᵈⁱˢᵗ, εᵈⁱˢᵗ: $(round(εᵈⁱˢᵗ, digits=8))")
        end
        𝓃ᵈⁱˢᵗ += 1
    end
end

# 5. Aggregate states 
function fnAggregateStates!(params, endo)

    # A. Unpacking business 
    @unpack Nᶻ, Nᵃ, Nˡ, δᵈⁱˢᵗ, a⃗ = params

    # B. Iterate to find the distribution and its marginal 
    fnForwardIteration!(params, endo)
    endo.g̃ .= dropdims(sum(endo.g, dims=4), dims=4)

    # D. Update the labour market 
    fnJobDestruction!(params, endo)
    fnUpdateLabourMarket!(params, endo)

    # E. Weight for the marginal z × a distribution
    ωᵃ = dropdims(sum(endo.g, dims=(3, 4)), dims=(3, 4))

    # F. Capital demand and supply 
    endo.Kˢ = sum(ωᵃ .* a⃗')
    endo.Kᵈ = sum(ωᵃ .* endo.𝐤 .* endo.𝐨)

    # G. Labour demand and supply 
    endo.Lˢ = endo.W
    endo.Lᵈ = sum(ωᵃ .* endo.𝐥 .* endo.𝐨)
end

# -----------------------------------------------
# C. Solving  
# -----------------------------------------------

# 1. Labour residual.
function fnLabourResidual!(w, params, endo, r, τ)
    # A. Update the state
    endo.wₜ = w
    endo.rₜ = r
    endo.τₜ = τ

    # B. Run the parallelized engine
    fnVFI!(params, endo)
    fnAggregateStates!(params, endo)
    
    # C. Return the Labor market error
    εᴸ      = (endo.Lᵈ / max(endo.Lˢ, 1e-4)) - 1.0
    println("→→→ Wage search        | w = $(round(w, digits=4)), εᴸ = $(round(εᴸ, digits=4))")
    return εᴸ
end

# 2. Labour market 
function fnBudgetResidual!(τ, params, endo, r)

    # A. Unpacking business 
    @unpack w̲, w̅, δᴸ = params

    # B. Find the wage that clears labor for THIS tax and interest rate
    wˣ      = find_zero(w -> fnLabourResidual!(w, params, endo, r, τ), (w̲, w̅), Bisection(),xatol = δᴸ)
    endo.wₜ = wˣ 

    # C. Return the Budget error
    εᵗ      = (endo.U * endo.wₜ / max(endo.τₜ, 1e-4)) - 1.0
    println("→→ Tax loop        | τ = $(round(τ, digits=4)), εᵗ = $(round(εᵗ, digits=4)) [Cleared w = $(round(endo.wₜ, digits=4))]")
    return εᵗ
end

# 3. Government budget 
function fnCapitalResidual!(params, endo, r)
        
        # A. Unpacking business 
        @unpack τ̲, τ̅, δᵗ = params

        # B. Solve
        τˣ      = find_zero(τ -> fnBudgetResidual!(τ, params, endo, r), (τ̲, τ̅), Bisection(), xatol = δᵗ)
        endo.τₜ = τˣ

        # C. Return error 
        εᴷ = (endo.Kᵈ / max(endo.Kˢ, 1e-4)) - 1.0
        println("→ Capital loop     | r = $(round(r, digits=4)), εᴷ = $(round(εᴷ, digits=4)) [Cleared τ = $(round(endo.τₜ, digits=4))]")
        return εᴷ
    end

# 4. Steady state 
function fnSolveSteadyState!(params, endo)
    
    # A. Unpacking business 
    @unpack r̲, r̅, δʳ = params

    # B. The Final Solve
    rˣ = find_zero(r -> fnCapitalResidual!(params, endo, r), (r̲, r̅), Bisection(), xatol = δʳ)
    endo.rₜ = rˣ

    # C. Lock in the results
    fnVFI!(params, endo)
    fnAggregateStates!(params, endo)
    println("\n--- Steady state ---")
    println("Wage (w):   $(round(endo.wₜ, digits=6))")
    println("Interest(r):$(round(endo.rₜ, digits=6))")
    println("Tax (τ):    $(round(endo.τₜ, digits=6))")
end

# # 1. Labour market clearing
# function fnLabourMarketClearing!(params, endo)
    
#     # A. Unpacking business 
#     @unpack w̲, w̅, δᴸ, λᵗ, τ̲, τ̅, κᵗ, δᵗ = params

#     # B. Initial setup 
#     endo.wₜ     = 1.6       
#     endo.τₜ     = 0.05 * endo.wₜ
#     εᵗ          = 1.0
#     w_low       = w̲         
#     w_high      = w̅    
    
#     # C. Tax loop 
#     while abs(εᵗ) > δᵗ
        
#         # I. Reset bisection bounds for the new tax environment     
#         εᴸ          = 1.0

#         # D. Wage bisection 
#         while abs(εᴸ) > δᴸ

#             # I. Run the model 
#             if (w_high - w_low) < 1e-4
#                 break
#             end
#             fnVFI!(params, endo)
#             fnAggregateStates!(params, endo)

#             # II. Calculate labour error 
#             εᴸ = (endo.Lᵈ / max(endo.Lˢ, 1e-4)) - 1.0
#             if abs(εᴸ) <= δᴸ
#                 break
#             end

#             # III. The bisection bounds update
#             if εᴸ > 0.0
#                 w_low = endo.wₜ
#             else
#                 w_high = endo.wₜ
#             end
#             endo.wₜ = 0.5 * w_low + 0.5 * w_high
#             println("→→→ Wage search        | w = $(round(endo.wₜ,digits=4)) [Floor: $(round(w_low,digits=3)), Ceil: $(round(w_high,digits=3))], εᴸ = $(round(εᴸ,digits=4))")
#         end

#         # E. Check the budget
#         εᵗ          = (endo.U * endo.wₜ / max(endo.τₜ, 1e-4)) - 1.0
#         τ_target    = clamp(endo.τₜ * (1.0 + κᵗ * εᵗ), τ̲, τ̅)
#         endo.τₜ     = λᵗ * endo.τₜ + (1 - λᵗ) * τ_target
#         println("→→ Tax loop            | τ = $(round(endo.τₜ,digits=4)), εᵗ = $(round(εᵗ,digits=4)) [Cleared w = $(round(endo.wₜ,digits=4))]")
    
#         # F. Update for performance 
#         w_low       = max(w̲, endo.wₜ * 0.90)
#         w_high      = min(w̅, endo.wₜ * 1.10)
#     end
# end

# # 2. Capital markets clearing 
# function fnCapitalMarketsClearing!(params, endo)

#     # A. Unpacking business 
#     @unpack r̲, r̅ = params

#     # B. Prepare closure 
#     function ℊ(r, params, endo)
#         endo.rₜ = r
#         fnLabourMarketClearing!(params, endo)
#         return (endo.Kᵈ / max(endo.Kˢ, 1e-4)) - 1.0
#     end

#     # C. Solve for interest 
#     rˣ = find_zero(r -> ℊ(r, params, endo), (r̲, r̅), Bisection())
#     endo.rₜ = rˣ
#     fnLabourMarketClearing!(params, endo)
#     println("→ Capital loop             | w = $(round(endo.wₜ,digits=3)), r = $(round(endo.rₜ,digits=3)), τ = $(round(endo.τₜ,digits=3))")
# end

# # 3. Solve the model 
# function fnSolveSteadyState!(params, endo)
#     fnCapitalMarketsClearing!(params, endo)
#     println("Steady state ready found")
#     println("Equilibrium values: w = $(endo.wₜ), r = $(endo.rₜ), τ = $(endo.τₜ)")
# end


# # 1A. Logistic function for bounding 
# function fnLogistic(x,a,b)
#     # a: minimum 
#     # b: maximum 
#     return a + (b - a) / (1 + exp(-x))
# end 

# # 1B. Inverse logistic for bounding 
# function fnInverseLogistic(y,a,b)
#     # a: minimum 
#     # b: maximum 
#     return log( (y - a) / (b - y))
# end 
# # 2. Combined steady state residuals
# function fnSteadyStateResiduals!(F, x, params, endo)

#     # A. Unpacking business 
#     @unpack w̲, w̅, r̲, r̅, τ̲, τ̅ = params 

#     # B. Unpack the simultaneous guesses
#     # x = [w_guess, r_guess, τ_guess]
#     endo.wₜ = fnLogistic(x[1],w̲,w̅)
#     endo.rₜ = fnLogistic(x[2],r̲,r̅)
#     endo.τₜ = fnLogistic(x[3],τ̲,τ̅)

#     # C. Run the heavy lifting 
#     fnVFI!(params, endo)
#     fnAggregateStates!(params, endo)

#     # D. Populate the residual vector F 
#     F[1] = (endo.Lᵈ / max(endo.Lˢ, 1e-4)) - 1.0                
#     F[2] = (endo.Kᵈ / max(endo.Kˢ, 1e-4)) - 1.0                
#     F[3] = (endo.τₜ / max(endo.wₜ * endo.U, 1e-4)) - 1.0

#     # E. Print the trace
#     println("→ Wage       | w = $(round(endo.wₜ,digits=4)), εᴸ = $(round(F[1],digits=4))")
#     println("→ Capital    | r = $(round(endo.rₜ,digits=4)), εᴷ = $(round(F[2],digits=4))")
#     println("→ Tax        | τ = $(round(endo.τₜ,digits=4)), εᵗ = $(round(F[3],digits=4))\n")

#     # E. The NaN Bouncer
#     if any(isnan, F) || any(isinf, F)
#         println("→→→ [WARNING] The world collapsed! Penalty.")
#         F .= 1e6 
#         return
#     end
# end 

# # 3. Solve for the steady state 
# function fnSolveSteadyState!(params, endo)

#     # A. Unpacking business 
#     @unpack w̲, w̅, r̲, r̅, τ̲, τ̅ = params 

#     # B. Set your initial guesses: [w_initial, r_initial, τ_initial]
#     x⁰          = [
#                 fnInverseLogistic(1.5, w̲,w̅),
#                 fnInverseLogistic(0.02,r̲,r̅),
#                 fnInverseLogistic(0.075,τ̲,τ̅)
#                 ] 

#     # C. Run the model
#     sol         = nlsolve((F, x) -> fnSteadyStateResiduals!(F, x, params, endo), x⁰, 
#                     method        = :broyden, 
#                     show_trace    = true, 
#                     ftol          = 1e-4)

#     # D. Run the multivariate solver
#     endo.wₜ = fnLogistic(sol.zero[1],w̲,w̅)
#     endo.rₜ = fnLogistic(sol.zero[2],r̲,r̅)
#     endo.τₜ = fnLogistic(sol.zero[3],τ̲,τ̅)

#     # E. Run the model one final time 
#     fnVFI!(params, endo)
#     fnAggregateStates!(params, endo)
#     println("Equilibrium values: w = $(endo.wₜ), r = $(endo.rₜ), τ = $(endo.τₜ)")
# end
