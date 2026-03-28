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

    # B. Workers VF 
    endo.𝐕ᵂ .= 0.1 * (1 - β)^(-1) * fnUtility(endo.rₜ .* a⃗' .- endo.τₜ .+ endo.wₜ, params)

    # C. Entrepreneurs VF 
    endo.𝐕ᴱ .= 0.1 * (1 - β)^(-1) * fnUtility(endo.rₜ .* a⃗' .- endo.τₜ .+ endo.Π, params)

    # D. Overall VF 
    endo.𝐕 .= max.(endo.𝐕ᵂ, endo.𝐕ᴱ)
    endo.𝐨 .= (endo.𝐕ᵂ .<= endo.𝐕ᴱ)

end

# 4. Solve for assets 
function fnFindAssets(iz, ia, params, endo)

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
        A_w = a⃗[1]
        iʷ = 1
        𝐕ᵂ = fnUtility(Cash_w - a⃗[1], params) + β * @views(endo.𝔼𝐕[iz, 1])
    else
        Objʷ = fnUtility(Cash_w .- a⃗, params) .+ β .* @views(endo.𝔼𝐕[iz, :])
        𝐕ᵂ, iʷ = findmax(Objʷ)
        A_w = a⃗[iʷ]
    end

    # E. Computation of assets: Entrepreneurs 
    if Upper_e <= Lower
        A_e = a⃗[1]
        iᵉ = 1
        𝐕ᴱ = fnUtility(Cash_e - a⃗[1], params) + β * @views(endo.𝔼𝐕[iz, 1])
    else
        Objᵉ = fnUtility(Cash_e .- a⃗, params) .+ β .* @views(endo.𝔼𝐕[iz, :])
        𝐕ᴱ, iᵉ = findmax(Objᵉ)
        A_e = a⃗[iᵉ]
    end
    return A_w, A_e, 𝐕ᵂ, 𝐕ᴱ, iʷ, iᵉ
end

# 5. Value function iteration 
function fnVFI!(params, endo)

    # A. Unpacking business
    @unpack ψ, μ⃗, δᵛᶠⁱ, 𝒾̄ᵛᶠⁱ, z⃗, c̲, a⃗, β, λᵛᶠⁱ = params

    # B. Static policies
    fnStaticPolicies!(params, endo)

    # C. Initial guess 
    fnInitialVFIGuess!(params, endo)

    # D. Prepare the VFI loop 
    εᵛᶠⁱ    = 1.0
    𝓃ᵛᶠⁱ    = 1
    𝐢ᵃʷ     = ones(Int, length(z⃗), length(a⃗))
    𝐢ᵃᵉ     = ones(Int, length(z⃗), length(a⃗))

    while (εᵛᶠⁱ > δᵛᶠⁱ && 𝓃ᵛᶠⁱ < 𝒾̄ᵛᶠⁱ)

        # D1. Update the expected VF and open the loop for productivity 
        endo.𝔼𝐕 .= ψ .* endo.𝐕 .+ (1 - ψ) .* (μ⃗' * endo.𝐕)
        𝐕ᵖʳᵉᵛ = copy(endo.𝐕)
        @inbounds for iz in eachindex(z⃗)

            # B1. Solve for each a and z 
            for ia in eachindex(a⃗)
                
                # Fully solve
                if (𝓃ᵛᶠⁱ % 25 == 0 | 𝓃ᵛᶠⁱ < 25)

                    # I. Find assets
                    endo.𝐚ʷ[iz, ia], endo.𝐚ᵉ[iz, ia], endo.𝐕ᵂ[iz, ia], endo.𝐕ᴱ[iz, ia],𝐢ᵃʷ[iz,ia],𝐢ᵃᵉ[iz,ia] = fnFindAssets(iz, ia, params, endo)

                    # II. VFs for different occupations
                    endo.𝐜ʷ[iz, ia] = a⃗[ia] * (1 + endo.rₜ) + endo.wₜ - endo.τₜ - endo.𝐚ʷ[iz, ia]
                    endo.𝐜ᵉ[iz, ia] = a⃗[ia] * (1 + endo.rₜ) + endo.Π[iz, ia] - endo.τₜ - endo.𝐚ᵉ[iz, ia]

                    # III. Occupational decision and updated VF
                    endo.𝐨[iz, ia] = (endo.𝐕ᴱ[iz, ia] >= endo.𝐕ᵂ[iz, ia])
                    endo.𝐚[iz, ia] = endo.𝐨[iz, ia] ? endo.𝐚ᵉ[iz, ia] : endo.𝐚ʷ[iz, ia]
                    endo.𝐜[iz, ia] = max(endo.𝐨[iz, ia] ? endo.𝐜ᵉ[iz, ia] :  endo.𝐜ʷ[iz, ia], c̲)
                    endo.𝐕[iz, ia] = max(endo.𝐕ᴱ[iz, ia], endo.𝐕ᵂ[iz, ia])
                else 
                    # Howard improvement
                    endo.𝐕ᴱ[iz, ia] = fnUtility(endo.𝐜ᵉ[iz, ia], params) + β * endo.𝔼𝐕[iz, 𝐢ᵃᵉ[iz,ia]]
                    endo.𝐕ᵂ[iz, ia] = fnUtility(endo.𝐜ʷ[iz, ia], params) + β * endo.𝔼𝐕[iz, 𝐢ᵃʷ[iz,ia]]
                    endo.𝐕[iz, ia]  = max(endo.𝐕ᴱ[iz, ia], endo.𝐕ᵂ[iz, ia])
                end 
            end
        end

        # D2. Error and update 
        εᵛᶠⁱ = maximum(abs.(𝐕ᵖʳᵉᵛ .- endo.𝐕))
        𝓃ᵛᶠⁱ += 1
        endo.𝐕 .= λᵛᶠⁱ .* 𝐕ᵖʳᵉᵛ .+ (1.0 - λᵛᶠⁱ) .* endo.𝐕
        # Add this heartbeat print statement!
        # if 𝓃ᵛᶠⁱ % 25 == 0
        #     println("→→→→ VFI Iteration: $𝓃ᵛᶠⁱ, εᵛᶠⁱ = $(round(εᵛᶠⁱ, digits=5))")
        # end
    end
    #println("→→→→ VFI done: 𝓃 = $𝓃ᵛᶠⁱ, ε = $εᵛᶠⁱ")
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
    while (εᵈⁱˢᵗ > δᵈⁱˢᵗ)

        # C. Aggregate states 
        fnJobDestruction!(params, endo)         # Getting JD, S 
        fnUpdateLabourMarket!(params, endo)     # Getting U, M, W, E  

        @inbounds for ia in eachindex(a⃗)
            for il in eachindex(l⃗)
                for iz in eachindex(z⃗)
                    for iu in 1:2

                        # I. Common terms
                        Mass = endo.g[iz, ia, il, iu]
                        if Mass == 0.0
                            continue
                        end
                        # (a) Next asset mass 
                        Aₜ = endo.𝐚[iz, ia]
                        ibᵃ = clamp(searchsortedlast(a⃗, Aₜ), 1, length(a⃗) - 1)
                        iuᵃ = ibᵃ + 1
                        wᵃ = clamp((a⃗[iuᵃ] - Aₜ) / (a⃗[iuᵃ] - a⃗[ibᵃ]), 0.0, 1.0)
                        # (b) Labour mass 
                        Lₜ = endo.𝐨[iz, ia] ? endo.𝐥[iz, ia] : 0.0
                        ibˡ = clamp(searchsortedlast(l⃗, Lₜ), 1, length(l⃗) - 1)
                        iuˡ = ibˡ + 1
                        wˡ = clamp((l⃗[iuˡ] - Lₜ) / (l⃗[iuˡ] - l⃗[ibˡ]), 0.0, 1.0)
                        # (c) Get flows that matter
                        f¹², f²¹, f²², f¹¹ = fnComputeFlows(iz, ia, il, iu, params, endo)
                        FlowU = f¹² + f²²
                        FlowO = f²¹ + f¹¹

                        # II. With switching productivity
                        # To unemployment (iu=2)
                        B[ibᵃ, ibˡ, 2] += Mass * wᵃ * wˡ * FlowU
                        B[ibᵃ, iuˡ, 2] += Mass * wᵃ * (1 - wˡ) * FlowU
                        B[iuᵃ, ibˡ, 2] += Mass * (1 - wᵃ) * wˡ * FlowU
                        B[iuᵃ, iuˡ, 2] += Mass * (1 - wᵃ) * (1 - wˡ) * FlowU
                        # To employment/entrepreneurship (iu=1)
                        B[ibᵃ, ibˡ, 1] += Mass * wᵃ * wˡ * FlowO
                        B[ibᵃ, iuˡ, 1] += Mass * wᵃ * (1 - wˡ) * FlowO
                        B[iuᵃ, ibˡ, 1] += Mass * (1 - wᵃ) * wˡ * FlowO
                        B[iuᵃ, iuˡ, 1] += Mass * (1 - wᵃ) * (1 - wˡ) * FlowO

                        # III. The same productivity
                        # (c) Update (asset, labour) = (BB, BU, UB, UU)
                        # To unemployment 
                        gⁿᵉˣᵗ[iz, ibᵃ, ibˡ, 2] += ψ * Mass * wᵃ * wˡ * FlowU
                        gⁿᵉˣᵗ[iz, ibᵃ, iuˡ, 2] += ψ * Mass * wᵃ * (1 - wˡ) * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ, ibˡ, 2] += ψ * Mass * (1 - wᵃ) * wˡ * FlowU
                        gⁿᵉˣᵗ[iz, iuᵃ, iuˡ, 2] += ψ * Mass * (1 - wᵃ) * (1 - wˡ) * FlowU
                        # To others 
                        gⁿᵉˣᵗ[iz, ibᵃ, ibˡ, 1] += ψ * Mass * wᵃ * wˡ * FlowO
                        gⁿᵉˣᵗ[iz, ibᵃ, iuˡ, 1] += ψ * Mass * wᵃ * (1 - wˡ) * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ, ibˡ, 1] += ψ * Mass * (1 - wᵃ) * wˡ * FlowO
                        gⁿᵉˣᵗ[iz, iuᵃ, iuˡ, 1] += ψ * Mass * (1 - wᵃ) * (1 - wˡ) * FlowO
                    end
                end
            end
        end

        # IV. Add the mass to the updated productivity states
        @inbounds for iu in 1:2, il in 1:Nˡ, ia in 1:Nᵃ, iz in 1:Nᶻ
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

# 1. Labour market clearing
function fnLabourMarketClearing!(params, endo)
    
    # A. Unpacking business 
    @unpack w̲, w̅, δᴸ, κᴸ, λᵗ, τ̲, τ̅, κᵗ, δᵗ = params

    # B. Initial guesses
    endo.wₜ     = 1.1
    endo.τₜ     = 0.05 * endo.wₜ
    εᵗ          = 1.0

    # C. Budget loop 
    while abs(εᵗ) > δᵗ
        
        κᴸˣ     = κᴸ
        εᴸ      = 1.0
        εᴸ²     = 1.0

        # D.Labour market loop 
        while abs(εᴸ) > δᴸ

            # I. Run the model 
            fnVFI!(params, endo)
            fnAggregateStates!(params, endo)

            # II. Wage error 
            εᴸ = (endo.Lᵈ / max(endo.Lˢ, 1e-4)) - 1.0
            if εᴸ * εᴸ² < 0.0
                κᴸˣ = max(κᴸˣ * 0.5, 0.01)
            end

            # III. Update Wage (Tax remains fixed here)
            endo.wₜ = clamp(endo.wₜ * (1.0 + κᴸˣ * εᴸ), w̲, w̅)
            εᴸ²     = εᴸ
            
            println("→→→ Wage Search    | w = $(round(endo.wₜ,digits=3)), εᴸ = $(round(εᴸ,digits=3))")
        end

        # E. OUTER UPDATE: Once wage has settled, check the budget
        # Calculate budget error (Tax Gap)
        εᵗ = endo.U * endo.wₜ - endo.τₜ
        
        # Proportional update for the Tax
        τ_next  = clamp(endo.τₜ * (1.0 + κᵗ * (εᵗ/max(endo.τₜ, 1e-4))), τ̲, τ̅)
        endo.τₜ = λᵗ * endo.τₜ + (1 - λᵗ) * τ_next
        println("→→ Tax loop            | w = $(round(endo.wₜ,digits=3)), εᵗ = $(round(εᵗ,digits=5)), τ = $(round(endo.τₜ,digits=3))")
    end
end


# 2. Capital markets clearing 
function fnCapitalMarketsClearing!(params, endo)

    # A. Unpacking business 
    @unpack r̲, r̅ = params

    # B. Prepare closure 
    function ℊ(r, params, endo)
        endo.rₜ = r
        fnLabourMarketClearing!(params, endo)
        return (endo.Kᵈ / max(endo.Kˢ, 1e-4)) - 1.0
    end

    # C. Solve for interest 
    rˣ = find_zero(r -> ℊ(r, params, endo), (r̲, r̅), Bisection())
    endo.rₜ = rˣ
    fnLabourMarketClearing!(params, endo)
    println("→ Capital loop             | w = $(round(endo.wₜ,digits=3)), r = $(round(endo.rₜ,digits=3)), τ = $(round(endo.τₜ,digits=3))")
end

# 3. Solve the model 
function fnSolveSteadyState!(params, endo)
    fnCapitalMarketsClearing!(params, endo)
    println("Steady state ready found")
    println("Equilibrium values: w = $(endo.wₜ), r = $(endo.rₜ), τ = $(endo.τₜ)")
end


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
#     endo.τₜ = min(endo.τₜ, 0.9*endo.wₜ)

#     # C. Run the heavy lifting 
#     fnVFI!(params, endo)
#     fnAggregateStates!(params, endo)

#     # D. Populate the residual vector F 
#     F[1] = (endo.Lᵈ / max(endo.Lˢ, 1e-4)) - 1.0                
#     F[2] = (endo.Kᵈ / max(endo.Kˢ, 1e-4)) - 1.0                
#     F[3] = (endo.τₜ / max(endo.wₜ * endo.U, 1e-4)) - 1.0
# end 

# # 3. Solve for the steady state 
# function fnSolveSteadyState!(params, endo)

#     # A. Unpacking business 
#     @unpack w̲, w̅, r̲, r̅, τ̲, τ̅ = params 

#     # B. Set your initial guesses: [w_initial, r_initial, τ_initial]
#     x⁰          = [
#                 fnInverseLogistic(1.0, w̲,w̅),
#                 fnInverseLogistic(0.01,r̲,r̅),
#                 fnInverseLogistic(0.25,τ̲,τ̅)
#                 ] 

#     # C. Define the closure 
#     function 𝒻_optim(x)
#         F = zeros(3)
#         fnSteadyStateResiduals!(F, x, params, endo)
#         return F[1]^2 + F[2]^2 + F[3]^2
#     end

#     # D. Run the derivative-free Nelder-Mead simplex
#     sol = optimize(𝒻_optim, x⁰, NelderMead(), Optim.Options(show_trace=true, iterations=500))

#     # D. Run the multivariate solver (defaults to a highly efficient Trust-Region method)
#     endo.wₜ = fnLogistic(sol.minimizer[1],w̲,w̅)
#     endo.rₜ = fnLogistic(sol.minimizer[2],r̲,r̅)
#     endo.τₜ = fnLogistic(sol.minimizer[3],τ̲,τ̅)

#     # E. Run the model one final time 
#     fnVFI!(params, endo)
#     fnAggregateStates!(params, endo)
#     println("Equilibrium values: w = $(endo.wₜ), r = $(endo.rₜ), τ = $(endo.τₜ)")
# end
